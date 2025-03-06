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
            self.data = self.expand_zones_and_arenas(pl.read_csv(f))

        with open(zantiks_file.genotypes_path, "r") as f:
            self.genotypes = pl.read_csv(f)

        self.genotype_filter = None

    def __repr__(self):
        return "\n" + str(self.genotypes) + "\n" + str(self.data)

    def set_genotype_filter(self, *args: str) -> None:
        if not args:
            self.genotype_filter = None
            return

        # else
        new_filter = []
        for arg in args:
            upper_arg = arg.upper()
            if upper_arg == "WT":
                new_filter.append("WT")
            elif upper_arg == "HOM":
                new_filter.append("HOM")
            elif upper_arg == "HET":
                new_filter.append("HET")
            else:
                print(f"{arg} is not a recognized genotype")
        if not new_filter:
            return

        self.genotype_filter = new_filter

    def expand_zones_and_arenas(self, df: pl.DataFrame):
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
                .cast(pl.Int64)
                .alias("ARENA"),
                pl.col("ZONE_CODE")
                .str.extract(r"Z(\d+)", 1)
                .cast(pl.Int64)
                .alias("ZONE"),
            ]
        )

        result = long_df.drop("ZONE_CODE").select(
            ["TIME", "BIN_NUMBER", "ARENA", "ZONE", "TIME_IN_ZONE"]
        )
        return result

    def filtered_data(self) -> pl.DataFrame:
        pass


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
        assay_types=("social_preference",), data_type="zantiks"
    )
    dfs = dl.load_all()
    print(dfs)
