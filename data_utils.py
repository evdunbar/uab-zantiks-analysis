#!/usr/bin/env uv run

import glob
import io
import re
from abc import ABC, abstractmethod
from typing import Iterable, Literal, Self, Sequence

import polars as pl

type VariableDate = Sequence[int, int | None, int | None]


class Assay(ABC):
    @property
    @abstractmethod
    def total_arenas(self) -> int:
        pass

    @property
    @abstractmethod
    def arenas_per_group(self) -> int | None:
        pass


class LightDarkPreference3wpf(Assay):
    def total_arenas(self) -> int:
        return 12

    def arenas_per_group(self) -> int | None:
        return 4


class LightDarkPreference6dpf(Assay):
    def total_arenas(self) -> int:
        return 12

    def arenas_per_group(self) -> int | None:
        return 12


class LightDarkTransition(Assay):
    def total_arenas(self) -> int:
        return 48

    def arenas_per_group(self) -> int | None:
        return None


class MirrorBiting(Assay):
    def total_arenas(self) -> int:
        return 20

    def arenas_per_group(self) -> int | None:
        return 4


class SocialPreference(Assay):
    def total_arenas(self) -> int:
        return 10

    def arenas_per_group(self) -> int | None:
        return 4


class StartleResponse(Assay):
    def total_arenas(self) -> int:
        return 48

    def arenas_per_group(self) -> int | None:
        return None


class Ymaze15(Assay):
    def total_arenas(self) -> int:
        return 15

    def arenas_per_group(self) -> int | None:
        return 12


class Ymaze4(Assay):
    def total_arenas(self) -> int:
        return 4

    def arenas_per_group(self) -> int | None:
        return 4


class ZantiksFile:
    assay_types = {
        "light_dark_preference": LightDarkPreference3wpf,
        "light_dark_transition": LightDarkTransition,
        "mirror_biting": MirrorBiting,
        "social_preference": SocialPreference,
        "startle_response": StartleResponse,
        "ymaze_15": Ymaze15,
        "ymaze_4": Ymaze4,
    }

    def __init__(self, path: str):
        self.path = path
        parts = path.split("/")
        self.filename = parts[-1]
        filename_parts = self.filename.split("-")
        self.groups = self._parse_groups(parts[-2])
        self.assay_type = self.assay_types[filename_parts[0]]()
        self.genotypes_path = self._find_run_day_path() + "/genotypes.csv"
        self.year = filename_parts[1][:4]
        self.month = filename_parts[1][4:6]
        self.day = filename_parts[1][6:8]

        self.is_xy = "xy" in self.filename or "XY" in self.filename

    def __repr__(self):
        return f"ZantiksFile({self.year}, {self.month}, {self.day}, {self.assay_type}, {self.groups}, {self.filename}) - is_xy: {self.is_xy}"

    def _parse_groups(self, groups: str) -> tuple[str]:
        if groups.startswith("zantiks"):
            return (groups[-1].upper(),)
        else:
            return tuple(groups.upper())

    def _find_run_day_path(self):
        date_folders_expression = re.compile(r".*\d{4}/\d{2}/\d{2}/")
        return date_folders_expression.search(self.path)[0]


class ZantiksData:
    def __init__(self, zantiks_file: ZantiksFile, use_genotypes: bool = True):
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

        if use_genotypes:
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
        pathglob = "**/"

        if groups:
            pathglob += f"{self._iterable_glob(groups)}/"

        if assay_types:
            pathglob += self._iterable_glob(assay_types)
        elif pathglob[-1] != "*":
            pathglob += "*"

        if dates:
            date_patterns = []
            for date in dates:
                year = str(date[0])
                month = f"{date[1]:02d}" if date[1] is not None else ""
                day = f"{date[2]:02d}" if date[2] is not None else ""
                date_patterns.append(f"{year}{month}{day}")

            pathglob += f"-{self._iterable_glob(date_patterns)}*"
        elif pathglob[-1] != "*":
            pathglob += "*"

        if pathglob[-1] != "*":
            pathglob += "*"
        pathglob += ".csv"

        matching_paths = glob.glob(pathglob, recursive=True)
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

    def load_all(self, *args, **kwargs) -> list[ZantiksData]:
        data = []
        for zantiks_file in self.zantiks_files:
            data.append(ZantiksData(zantiks_file, *args, **kwargs))

        return data


if __name__ == "__main__":
    dl = DataLoader().add_by_filter(
        assay_types=("light_dark_transition",), data_type="position"
    )
    dfs = dl.load_all()
    for df in dfs:
        print(df.info)
        # print(df.get_genotype(("WT", "HOM")))
