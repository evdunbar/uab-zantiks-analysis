#!/usr/bin/env uv run
import math
from typing import Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from numpy.typing import NDArray
from PIL import Image
from scipy import ndimage

import data_utils

mpl.rc("text", usetex=True)


class Heatmap:
    def __init__(
        self,
        data: data_utils.ZantiksData,
        genotypes: str | Iterable[str],
        arena_map_path: str,
        scale_factors: tuple[NDArray],
    ):
        self.info = data.info
        self.xy_data = data.get_genotype(genotypes)
        arena_map = self._load_image(arena_map_path)
        self.arena_mask = arena_map.sum(-1) < 10
        self.x_bias = scale_factors[0][0]
        self.x_scale = scale_factors[0][1]
        self.y_bias = scale_factors[1][0]
        self.y_scale = scale_factors[1][1]

        self.map = np.zeros_like(self.arena_mask, dtype=float)

    def _load_image(self, filepath: str) -> NDArray:
        with open(filepath, "rb") as f:
            image = np.asarray(Image.open(f))

        return image

    def make_map(self, *, by_arena: bool = False, show: bool = False):
        for arena_id in self.xy_data["ARENA"].unique().sort():
            arena = np.nan_to_num(
                self.xy_data.filter(pl.col("ARENA") == arena_id)
                .select("X", "Y")
                .to_numpy()
                .astype(float)
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
        if self.info.groups is not None:
            plot_title += f" \\textit{{Group(s) {self.info.groups}}}"
            save_title += f"_{self.info.groups}"
        if self.info.day is not None:
            plot_title += f" {self.info.month}-{self.info.day}-{self.info.year}"
            save_title += f"_{self.info.month}-{self.info.day}-{self.info.year}"

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
    import find_mapping

    dataloader = data_utils.DataLoader().add_by_filter(
        assay_types=("social_preference",), data_type="position"
    )
    all_zantiks_data = dataloader.load_all()
    mapping_finder = find_mapping.MapFinder(
        [zantiks_data.info.path for zantiks_data in all_zantiks_data]
    )
    mapping_finder.find_mappings()
    scale_factors = mapping_finder.mappings

    for zantiks_data in all_zantiks_data:
        Heatmap(
            zantiks_data,
            (
                "WT",
                "HOM",
            ),
            "data/social_preference_arenas.bmp",
            scale_factors,
        ).make_map(by_arena=True, show=True)
