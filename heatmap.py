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
    def __init__(
        self, data_path: str, arena_map_path: str, scale_factors: tuple[NDArray]
    ):
        self.xy_data = self._load_xy(data_path)
        arena_map = self._load_image(arena_map_path)
        self.arena_mask = arena_map.sum(-1) < 10
        self.x_bias = scale_factors[0][0]
        self.x_scale = scale_factors[0][1]
        self.y_bias = scale_factors[1][0]
        self.y_scale = scale_factors[1][1]

        self.num_arenas = int((self.xy_data.width - 1) / 2)
        self.map = np.zeros_like(self.arena_mask, dtype=float)

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
            scaled_x = math.floor(x * self.x_scale + self.x_bias)
            scaled_y = math.floor(y * self.y_scale + self.y_bias)
            self.map[scaled_y, scaled_x] += 1

        return self.map

    def show_map(self, sum_radius: float = 10) -> None:
        masked_plot = np.ma.masked_array(
            self.map, mask=np.bitwise_invert(self.arena_mask)
        )
        masked_plot[self.arena_mask] = 0.0
        plot_buffer = self.map.copy()
        plot_buffer = ndimage.gaussian_filter(plot_buffer, sum_radius)

        plot_title = r"\textbf{Social Preference}"
        save_title = "social_preference"
        if self.group is not None:
            plot_title += f" \\textit{{Group(s) {self.group[0].upper()}}}"
            save_title += f"_{self.group[0].upper()}"
        if self.date is not None:
            plot_title += f" {self.date[1]}-{self.date[2]}-{self.date[0]}"
            save_title += f"_{self.date[1]}-{self.date[2]}-{self.date[0]}"

        plt.figure(dpi=100)
        plt.title(plot_title)
        im = plt.imshow(plot_buffer, cmap="viridis")
        plt.imshow(masked_plot, cmap="binary")
        plt.colorbar(
            im, label="Frequency", orientation="vertical", shrink=0.805, aspect=13
        )
        plt.xticks([])
        plt.yticks([])
        plt.gca().set_xticklabels([])
        plt.gca().set_yticklabels([])
        plt.box(on=True)
        plt.savefig(
            f"data/figures/{save_title}.png",
            dpi=500,
            bbox_inches="tight",
        )
        plt.show()


if __name__ == "__main__":
    import glob

    import find_mapping

    filenames = glob.glob("data/*/*/*/social_preference/*/*xy_position.csv")
    mapping_finder = find_mapping.MapFinder(filenames)
    mapping_finder.find_mappings()
    scale_factors = mapping_finder.mappings

    for filename in filenames:
        Heatmap(filename, "data/social_preference_arenas.bmp", scale_factors).make_map(
            show=True
        )
