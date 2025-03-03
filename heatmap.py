#!/usr/bin/env uv run
import math

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from numpy.typing import NDArray
from scipy import signal


class Heatmap:
    def __init__(self, path: str):
        self.xy_data = self._load_xy(path)
        self.num_arenas = int((self.xy_data.width - 1) / 2)

        width = 12776
        height = 8540
        self.map = np.zeros((height, width))

    def _load_xy(self, filepath: str) -> pl.DataFrame:
        """
        - Returns columns of the form [XY]_A[0-9]+
        - Represents X and Y tracking for Arenas
        - Data in millimeters
        - Top left of plate is (0, 0)
        - Bottom right is (127.76, 85.4)
        """
        with open(filepath, "r") as f:
            xy_data = pl.read_csv(f)

        return xy_data

    def make_map(self) -> NDArray:
        for i in range(1, self.num_arenas + 1):
            arena = np.nan_to_num(
                self.xy_data.select(f"X_A{i}", f"Y_A{i}").to_numpy().astype(float)
            )
            self.map += self._process_arena(arena)

        return self.map

    def _process_arena(self, arena_data: NDArray) -> NDArray:
        for x, y in arena_data:
            if x <= 0.0 and y <= 0.0:
                continue
            scaled_x = math.floor(x * 100)
            scaled_y = math.floor(y * 100)
            self.map[scaled_y, scaled_x] += 1

        return self.map

    def show_map(self) -> None:
        self.map = signal.decimate(self.map, q=10, axis=0, n=3)
        self.map = signal.decimate(self.map, q=10, axis=1, n=3)
        self.map = signal.decimate(self.map, q=10, axis=0, n=3)
        self.map = signal.decimate(self.map, q=10, axis=1, n=3)

        plt.figure(dpi=100)
        plt.imshow(self.map, cmap="magma", vmin=0.0)
        plt.show()


if __name__ == "__main__":
    import glob

    filenames = glob.glob("*/*/*/social_preference/*/*xy_position.csv")
    for filename in filenames:
        heatmap = Heatmap(filename)
        heatmap.make_map()
        heatmap.show_map()
