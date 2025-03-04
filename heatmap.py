#!/usr/bin/env uv run
import math
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from numpy.typing import NDArray
from PIL import Image
from scipy import ndimage

mpl.rc("text", usetex=True)


class Heatmap:
    def __init__(self, data_path: str, arena_map_path: str):
        self.xy_data = self._load_xy(data_path)
        arena_map = self._load_image(arena_map_path)
        self.arena_mask = arena_map.sum(-1) < 10

        self.num_arenas = int((self.xy_data.width - 1) / 2)
        self.map = np.zeros_like(self.arena_mask, dtype=float)

        self.scale = 10.02272713
        self.bias = 76.05767812
        self.date = self._value_from_path(data_path, r".*(\d{4})(\d{2})(\d{2})T")
        self.group = self._value_from_path(data_path, r".*/([abcd]+)/.*")

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

    def _load_image(self, filepath: str) -> NDArray:
        with open(filepath, "rb") as f:
            image = np.asarray(Image.open(f))

        return image

    def _value_from_path(self, path: str, query: str) -> tuple | None:
        regexp = re.compile(query)
        if (maybe_matches := regexp.search(path)) is not None:
            return maybe_matches.groups()
        else:
            return None

    def make_map(self, *, show: bool = False):
        for i in range(1, self.num_arenas + 1):
            self.show_map()
            arena = np.nan_to_num(
                self.xy_data.select(f"X_A{i}", f"Y_A{i}").to_numpy().astype(float)
            )
            self._process_arena(arena)

        if show:
            self.show_map()

    def _process_arena(self, arena_data: NDArray) -> NDArray:
        for x, y in arena_data:
            if x <= 0.0 and y <= 0.0:
                continue
            scaled_x = math.floor(x * self.scale + self.bias)
            scaled_y = math.floor(y * self.scale + self.bias)
            self.map[scaled_y, scaled_x] += 1

        return self.map

    def show_map(self, sum_radius: float = 10) -> None:
        masked_plot = np.ma.masked_array(
            self.map, mask=np.bitwise_invert(self.arena_mask)
        )
        masked_plot[self.arena_mask] = 0.0
        plot_buffer = self.map.copy()
        plot_buffer = ndimage.gaussian_filter(plot_buffer, sum_radius)
        # plot_buffer[self.arena_mask] = 0.0

        title = r"\textbf{Social Preference}"
        if self.group is not None:
            title += f" \\textit{{Group(s) {self.group[0].upper()}}}"
        if self.date is not None:
            title += f" {self.date[1]}-{self.date[2]}-{self.date[0]}"

        plt.figure(dpi=100)
        plt.title(title)
        plt.imshow(plot_buffer, cmap="magma")
        plt.imshow(masked_plot, cmap="Greens_r")
        plt.axis("off")
        plt.show()


if __name__ == "__main__":
    import glob

    filenames = glob.glob("data/*/*/*/social_preference/*/*xy_position.csv")
    # for filename in filenames:
    for filename in (filenames[2],):
        Heatmap(filename, "data/social_preference_arenas.bmp").make_map(show=True)
