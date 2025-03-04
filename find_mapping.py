#!/usr/bin/env uv run
from typing import Literal

import numpy as np
import polars as pl
from numpy.typing import NDArray
from PIL import Image
from polars import selectors as cs


class MapFinder:
    def __init__(self, filenames: list[str]):
        # parameters
        self.pos_dfs = [self._load_dataframe(filename) for filename in filenames]

        # constant
        self.arena_map = (
            self._load_image("data/social_preference_arenas.bmp").sum(-1) > 10
        )

        # computed
        self.min_dfs = [df.min() for df in self.pos_dfs]
        self.max_dfs = [df.max() for df in self.pos_dfs]
        self.x_mins = self.direction_to_numpy(self.min_dfs, "X")
        self.x_maxs = self.direction_to_numpy(self.max_dfs, "X")
        self.y_mins = self.direction_to_numpy(self.min_dfs, "Y")
        self.y_maxs = self.direction_to_numpy(self.max_dfs, "Y")

    def _load_dataframe(self, filepath: str) -> pl.DataFrame:
        with open(filepath, "r") as f:
            df = pl.read_csv(f)
        return df.drop("RUNTIME")

    def _load_image(self, filepath: str) -> NDArray:
        with open(filepath, "rb") as f:
            image = np.asarray(Image.open(f))
        return image

    def direction_to_numpy(
        self, dfs: list[pl.DataFrame], direction: Literal["X", "Y"]
    ) -> NDArray:
        arrays = [
            df.select(cs.starts_with(direction)).to_numpy().astype(float) for df in dfs
        ]
        return np.ma.masked_array(arrays, np.isnan(arrays))

    def arena_mins_maxs(self) -> tuple[NDArray, NDArray, NDArray, NDArray]:
        x_mins = []
        x_maxs = []
        y_mins = []
        y_maxs = []
        was_in_arena = False
        for i in range(self.arena_map.shape[1]):
            col = self.arena_map[:, i]
            in_arena = np.sum(col) > 0
            if in_arena and not was_in_arena:
                was_in_arena = True
                x_mins.append(i)
            elif was_in_arena and not in_arena:
                was_in_arena = False
                x_maxs.append(i)

        for x_min, x_max in zip(x_mins, x_maxs):
            was_in_arena = False
            for i, row in enumerate(self.arena_map[:, x_min:x_max]):
                in_arena = np.sum(row) > 0
                if in_arena and not was_in_arena:
                    was_in_arena = True
                    y_mins.append(i)
                elif was_in_arena and not in_arena:
                    was_in_arena = False
                    y_maxs.append(i)

        x_mins = np.asarray(x_mins)
        x_maxs = np.asarray(x_maxs)
        y_mins = np.asarray(y_mins)
        y_maxs = np.asarray(y_maxs)
        return x_mins, x_maxs, y_mins, y_maxs


if __name__ == "__main__":
    import glob

    filenames = glob.glob("data/*/*/*/social_preference/*/*xy_position.csv")
    mf = MapFinder(filenames)
    print(mf.arena_mins_maxs())
