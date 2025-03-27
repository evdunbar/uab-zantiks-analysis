#!/usr/bin/env uv run
import math
from typing import Iterable

import cv2
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
    ):
        self.info = data.info
        self.xy_data = data.get_genotype(genotypes)
        self.arena_ids = self.xy_data["ARENA"].unique().sort()

        self.arena_map = self._load_image(data.info.assay_type.arena_bmp_path)
        self.arena_mask = np.zeros(
            (self.arena_map.shape[0], self.arena_map.shape[1])
        ).astype(bool)
        for arena in self.arena_ids:
            self.arena_mask = self.arena_mask | self.get_arena_coords(arena - 1)

        self.camera_params = self._load_camera_params(
            data.info.assay_type.camera_parameters_path
        )

        self.map = np.zeros_like(self.arena_mask).astype(float)

    def _load_image(self, filepath: str) -> NDArray:
        with open(filepath, "rb") as f:
            image = np.asarray(Image.open(f))

        cleaned_image = np.where(image < 8, 0, image)
        return cleaned_image

    def _load_camera_params(self, filepath: str):
        camera_params = {}
        with open(filepath, "rb") as f:
            camera_params_file = np.load(f)
            camera_params["rvec"] = camera_params_file["rvec"]
            camera_params["tvec"] = camera_params_file["tvec"]
            camera_params["mtx"] = camera_params_file["mtx"]
            camera_params["dist"] = camera_params_file["dist"]

        return camera_params

    def _arena_color(self, arena: int) -> NDArray:
        # fmt: off
        color_cycle = np.asarray(((120, 120, 248),  # blue
                                  (120, 248, 120),  # green
                                  (120, 248, 248),  # cyan
                                  (248, 120, 120),  # red
                                  (248, 120, 248),  # pink
                                  (248, 248, 120))) # yellow
        # fmt: on
        num_colors = 6
        darken_by = 8

        times_to_darken, color = divmod(arena, num_colors)
        return color_cycle[color] - times_to_darken * darken_by

    def get_arena_coords(self, arena: int) -> NDArray[bool]:
        return np.all(self.arena_map == self._arena_color(arena), axis=-1)

    def make_map(self, *, by_arena: bool = False, show: bool = False):
        for arena_id in self.arena_ids:
            dirty_arena_points = (
                self.xy_data.filter(pl.col("ARENA") == arena_id)
                .select("X", "Y")
                .to_numpy()
                .astype(float)
            )
            dirty_3d_points = np.column_stack(
                (dirty_arena_points, np.zeros(dirty_arena_points.shape[0]))
            ).astype(np.float32)
            arena_points = np.nan_to_num(
                np.squeeze(
                    cv2.projectPoints(
                        dirty_3d_points,
                        self.camera_params["rvec"],
                        self.camera_params["tvec"],
                        self.camera_params["mtx"],
                        self.camera_params["dist"],
                    )[0]  # first value is the points
                )
            )
            self._process_arena(arena_points)

        if show:
            self.show_map()

    def _process_arena(self, arena_data: NDArray) -> NDArray:
        for x, y in arena_data:
            if x <= 0.0 and y <= 0.0 or x > self.map.shape[1] or y > self.map.shape[0]:
                continue
            self.map[math.floor(y), math.floor(x)] += 1

        return self.map

    def show_map(self, sum_radius: float = 10) -> None:
        masked_plot = np.ma.masked_array(self.map, mask=self.arena_mask)
        masked_plot[np.logical_not(self.arena_mask)] = 0.0
        plot_buffer = self.map.copy()
        plot_buffer = ndimage.gaussian_filter(plot_buffer, sum_radius)

        plot_title = f"\\textbf{{{self.info.assay_type.pretty_name}}}"
        save_title = self.info.assay_type.name
        if self.info.groups is not None:
            plot_title += f" \\textit{{{self._format_groups()}}}"
            save_title += f"_{''.join(self.info.groups)}"
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

    def _format_groups(self) -> str | None:
        if not self.info.groups:
            return None
        else:
            if len(self.info.groups) == 1:
                return f"Group {self.info.groups[0]}"
            elif len(self.info.groups) == 2:
                return f"Groups {r' \& '.join(self.info.groups)}"
            else:
                return f"Groups {', '.join(self.info.groups)}"


if __name__ == "__main__":
    # fmt: off
    dataloader = data_utils.DataLoader().add_by_filter(
        # assay_types=("light_dark_preference_3wpf",), data_type="position"
        # assay_types=("light_dark_preference_6dpf",), data_type="position"
        # assay_types=("light_dark_transition",), data_type="position"
        # assay_types=("mirror_biting",), data_type="position"
        # assay_types=("social_preference",), data_type="position"
        # assay_types=("startle_response",), data_type="position"
        # assay_types=("ymaze_15",), data_type="position"
        assay_types=("ymaze_4",), data_type="position"
    )
    # fmt: on
    all_zantiks_data = dataloader.load_all()
    for zantiks_data in all_zantiks_data:
        Heatmap(zantiks_data, ("WT", "HOM", "HET")).make_map(by_arena=True, show=True)
