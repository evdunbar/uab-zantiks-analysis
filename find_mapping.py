#!/usr/bin/env uv run
from typing import Iterable

import alphashape
import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
from PIL import Image

import data_utils


class MapFinder:
    def __init__(self, data: Iterable[data_utils.ZantiksData]):
        # parameters
        self.data: Iterable[data_utils.ZantiksData] = []
        for datum in data:
            if datum.info.is_xy:
                self.data.append(datum)
            else:
                print(f"{datum} is not position data")

        # constant
        self.arena_map = (
            self._load_image("data/social_preference_arenas.bmp").sum(-1) > 10
        )

        # computed
        # self.min_dfs = [df.min() for df in self.pos_dfs]
        # self.max_dfs = [df.max() for df in self.pos_dfs]
        # self.x_mins = self.direction_to_numpy(self.min_dfs, "X").min(axis=0)
        # self.x_maxs = self.direction_to_numpy(self.max_dfs, "X").max(axis=0)
        # self.y_mins = self.direction_to_numpy(self.min_dfs, "Y").min(axis=0)
        # self.y_maxs = self.direction_to_numpy(self.max_dfs, "Y").max(axis=0)
        # self.xs = np.ma.concatenate((self.x_mins, self.x_maxs))
        # self.ys = np.ma.concatenate((self.y_mins, self.y_maxs))
        # arena_x_mins, arena_x_maxs, arena_y_mins, arena_y_maxs = (
        #     self.find_arena_mins_maxs()
        # )
        # self.arena_xs = np.concatenate((arena_x_mins, arena_x_maxs))
        # self.arena_ys = np.concatenate((arena_y_mins, arena_y_maxs))

    def _load_image(self, filepath: str) -> NDArray:
        with open(filepath, "rb") as f:
            image = np.asarray(Image.open(f))
        return image

    def scatter(self, separate: bool):
        if separate:
            for datum in self.data:
                plt.figure(dpi=200)
                plt.title(datum.info.path)
                plt.scatter(datum.data["X"], datum.data["Y"])
                plt.show()
        else:
            plt.figure(dpi=200)
            for datum in self.data:
                plt.scatter(datum.data["X"], datum.data["Y"])

            plt.show()

    def make_shape(self):
        # shape = alphashape.alphashape()

    # def direction_to_numpy(
    #     self, dfs: list[pl.DataFrame], direction: Literal["X", "Y"]
    # ) -> NDArray:
    #     arrays = [
    #         df.select(cs.starts_with(direction)).to_numpy().astype(float) for df in dfs
    #     ]
    #     return np.ma.masked_array(arrays, np.isnan(arrays)).squeeze()
    #
    # def find_arena_mins_maxs(self) -> tuple[NDArray, NDArray, NDArray, NDArray]:
    #     x_mins = []
    #     x_maxs = []
    #     y_mins = []
    #     y_maxs = []
    #     was_in_arena = False
    #     for i in range(self.arena_map.shape[1]):
    #         col = self.arena_map[:, i]
    #         in_arena = np.sum(col) > 0
    #         if in_arena and not was_in_arena:
    #             was_in_arena = True
    #             x_mins.append(i)
    #         elif was_in_arena and not in_arena:
    #             was_in_arena = False
    #             x_maxs.append(i)
    #
    #     for x_min, x_max in zip(x_mins, x_maxs):
    #         was_in_arena = False
    #         for i, row in enumerate(self.arena_map[:, x_min:x_max]):
    #             in_arena = np.sum(row) > 0
    #             if in_arena and not was_in_arena:
    #                 was_in_arena = True
    #                 y_mins.append(i)
    #             elif was_in_arena and not in_arena:
    #                 was_in_arena = False
    #                 y_maxs.append(i)
    #
    #     x_mins = np.asarray(x_mins)
    #     x_maxs = np.asarray(x_maxs)
    #     y_mins = np.asarray(y_mins)
    #     y_maxs = np.asarray(y_maxs)
    #     return x_mins, x_maxs, y_mins, y_maxs
    #
    # def find_mappings(self):
    #     self.mappings = []
    #     for positions, arena_positions in (
    #         (self.xs, self.arena_xs),
    #         (self.ys, self.arena_ys),
    #     ):
    #         mask = np.bitwise_not(positions.mask)
    #         positions = positions.data[mask]
    #         arena_positions = arena_positions[mask]
    #
    #         # a + b(positions) = arena_positions
    #         m = positions[:, np.newaxis] ** (0, 1)
    #         least_squares = linalg.lstsq(m, arena_positions)
    #         coefficients = least_squares[0]
    #         self.mappings.append(coefficients)


if __name__ == "__main__":
    dl = data_utils.DataLoader().add_by_filter(
        assay_types=("ymaze_4",), data_type="position"
    )
    data = dl.load_all(use_genotypes=False)
    # for datum in data:
    #     print(datum.info.path)
    mf = MapFinder(data)
    mf.make_shape()
