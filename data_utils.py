#!/usr/bin/env uv run

import glob
import io
import re
from typing import Iterable, Literal, Self, Sequence

import polars as pl

type VariableDate = Sequence[int, int | None, int | None]


class ZantiksFile:
    def __init__(self, path: str):
        self.path = path
        parts = path.split("/")
        self.filename = parts.pop()
        self.groups = tuple(parts.pop().upper())
        self.assay_type = parts.pop()
        self.genotypes_path = "/".join(parts) + "/genotypes.csv"
        self.day = parts.pop()
        self.month = parts.pop()
        self.year = parts.pop()
        self.prefix = "".join(parts)

        self.is_xy = "xy" in self.filename or "XY" in self.filename

    def __repr__(self):
        return f"ZantiksFile({self.prefix}, {self.year}, {self.month}, {self.day}, {self.assay_type}, {self.groups}, {self.filename}) - is_xy: {self.is_xy}"

    def __str__(self):
        return self.path


class ZantiksData:
    def __init__(self, zantiks_file: ZantiksFile):
        self.info = zantiks_file

        with open(zantiks_file.path, "r") as f:
            if not zantiks_file.is_xy:
                lines = f.readlines()
                f = io.StringIO("".join(lines[3:-1]))
            data = pl.read_csv(f)
        if zantiks_file.is_xy:
            self.data = self._expand_arenas(data)
        else:
            self.data = self._expand_zones_and_arenas(data)

        with open(zantiks_file.genotypes_path, "r") as f:
            self.genotypes: pl.DataFrame = pl.read_csv(f).drop_nulls("Cluster")

        self._attach_genotypes()

    def __repr__(self):
        return "\n" + str(self.genotypes) + "\n" + str(self.data)

    def _expand_zones_and_arenas(self, df: pl.DataFrame) -> pl.DataFrame:
        pattern = re.compile(r"T\.A(\d+)\.Z(\d+)")

        long_df = df.unpivot(
            index=["TIME", "BIN_NUMBER"],
            on=[col for col in df.columns if re.match(pattern, col)],
            variable_name="ZONE_CODE",
            value_name="TIME_IN_ZONE",
        )
        long_df = long_df.with_columns(
            [
                pl.col("ZONE_CODE")
                .str.extract(r"T\.A(\d+)", 1)
                .cast(pl.Int32)
                .alias("ARENA"),
                pl.col("ZONE_CODE")
                .str.extract(r"Z(\d+)", 1)
                .cast(pl.Int32)
                .alias("ZONE"),
            ]
        )
        long_df = long_df.drop("ZONE_CODE")

        zero_arenas = (
            long_df.group_by("ARENA")
            .agg(pl.sum("TIME_IN_ZONE").alias("TOTAL_TIME"))
            .filter(pl.col("TOTAL_TIME") == 0)
            .select("ARENA")
        )
        if not zero_arenas.is_empty():
            zero_arena_list = zero_arenas["ARENA"].to_list()
            result = long_df.filter(~pl.col("ARENA").is_in(zero_arena_list))
        else:
            result = long_df

        result = result.select(["TIME", "BIN_NUMBER", "ARENA", "ZONE", "TIME_IN_ZONE"])
        return result

    def _expand_arenas(self, df: pl.DataFrame) -> pl.DataFrame:
        pattern = re.compile(r"[XY]_A(\d+)")

        long_df = df.unpivot(
            index="RUNTIME", on=[col for col in df.columns if re.match(pattern, col)]
        )
        long_df = long_df.with_columns(
            [
                pl.col("variable")
                .str.extract(r"(\d+)", 1)
                .cast(pl.Int32)
                .alias("ARENA"),
                pl.col("variable").str.extract(r"([XY])", 1).alias("direction"),
            ]
        )
        long_df = long_df.filter(pl.col("direction") == "X").join(
            long_df.filter(pl.col("direction") == "Y"),
            on=("ARENA", "RUNTIME"),
        )
        long_df = long_df.select(
            "RUNTIME",
            "ARENA",
            X=pl.col("value").cast(pl.Float64),
            Y=pl.col("value_right").cast(pl.Float64),
        )

        zero_arenas = (
            long_df.group_by("ARENA")
            .sum()
            .filter(pl.col("X") <= 0.0)
            .filter(pl.col("Y") <= 0.0)
            .select("ARENA")
        )
        if not zero_arenas.is_empty():
            zero_arena_list = zero_arenas["ARENA"].to_list()
            result = long_df.filter(~pl.col("ARENA").is_in(zero_arena_list))
        else:
            result = long_df

        return result

    def _attach_genotypes(self) -> None:
        if "A" in self.info.groups:
            genotype_mapping = self.genotypes.select("Cluster").with_row_index(offset=1)
        elif self.info.groups == ("B",):
            genotype_mapping = self.genotypes.select(
                (pl.int_range(pl.len(), dtype=pl.Int32) - 2).alias("index"), "Cluster"
            )
        elif self.info.groups == ("C", "D"):
            genotype_length = len(self.genotypes["Cluster"])
            genotype_mapping = self.genotypes.select(
                (
                    pl.int_range(pl.len(), dtype=pl.Int32) - (genotype_length // 2) + 1
                ).alias("index"),
                "Cluster",
            )
        else:
            raise ValueError(
                f"Zantiks data groups not understood. File: {self.info.path}, Groups: {self.info.groups}."
            )
        self.data = self.data.join(
            genotype_mapping, left_on="ARENA", right_on="index", maintain_order="left"
        )

    def get_genotype(self, genotype: str | Iterable[str]) -> pl.DataFrame:
        if genotype is str:
            genotype = (genotype,)

        return self.data.filter(pl.col("Cluster").is_in(genotype))


class DataLoader:
    def __init__(self):
        self.pathnames: list[str] = []
        self.zantiks_files = []

    def add_by_filter(
        self,
        *,
        dates: Iterable[VariableDate] | None = None,
        assay_types: Iterable[str] | None = None,
        groups: Iterable[str] | None = None,
        data_type: Literal["zantiks"] | Literal["position"] | None = None,
    ) -> Self:
        pathglob = "data/"

        if dates:
            date_patterns = []
            for date in dates:
                year = str(date[0])
                month = f"{date[1]:02d}" if date[1] is not None else "*"
                day = f"{date[2]:02d}" if date[2] is not None else "*"
                date_patterns.append(f"{year}/{month}/{day}")

            pathglob += f"{self._iterable_glob(date_patterns)}/"
        else:
            pathglob += "*/*/*/"

        if assay_types:
            pathglob += f"{self._iterable_glob(assay_types)}/"
        else:
            pathglob += "*/"

        if groups:
            pathglob += f"{self._iterable_glob(groups)}/"
        else:
            pathglob += "*/"

        # can't parse data_type yet
        pathglob += "*.csv"

        # collect files
        matching_paths = glob.glob(pathglob)

        if data_type == "zantiks":
            matching_paths = [path for path in matching_paths if "xy" not in path]
        elif data_type == "position":
            matching_paths = [path for path in matching_paths if "xy" in path]

        self.add_by_name(matching_paths)

        return self

    def _iterable_glob(self, strings: Iterable[str]) -> str:
        if len(tuple(strings)) == 1:
            return strings[0]
        else:
            return f"{{{','.join(strings)}}}"

    def add_by_glob(self, pathglob: str) -> Self:
        self.add_by_name(glob.glob(pathglob))

    def add_by_name(self, paths: str | Iterable[str]) -> Self:
        if paths is str:
            paths = (paths,)

        self.pathnames += paths
        self.zantiks_files += [ZantiksFile(path) for path in paths]

        return self

    def load_all(self) -> list[ZantiksData]:
        data = []
        for zantiks_file in self.zantiks_files:
            data.append(ZantiksData(zantiks_file))

        return data


if __name__ == "__main__":
    dl = DataLoader().add_by_filter(
        assay_types=("social_preference",), data_type="position"
    )
    dfs = dl.load_all()
    for df in dfs:
        print(df.info)
        print(df.get_genotype(("WT", "HOM")))
