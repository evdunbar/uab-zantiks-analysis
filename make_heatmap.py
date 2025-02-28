#!/usr/bin/env uv run
import math

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from numpy.typing import NDArray
from scipy import ndimage


def load_xy(filepath: str) -> pl.DataFrame:
    with open(filepath, "r") as f:
        xy_data = pl.read_csv(f)

    return xy_data


def make_map(xy_data: NDArray) -> NDArray:
    heatmap = np.zeros((1080, 1440))

    for x, y in xy_data:
        if x <= 0.0 and y <= 0.0:
            continue
        scaled_x = math.floor(x * 1440 / 100)
        scaled_y = math.floor(y * 1080 / 100)
        heatmap[scaled_y, scaled_x] += 1

    return heatmap


def show_map(map: NDArray) -> None:
    # mask = map > 0.0
    # mean = map[mask].mean()
    # std = map[mask].std()
    map = ndimage.gaussian_filter(map, 10)

    plt.figure(dpi=300)
    # plt.imshow(map, cmap="magma", vmin=0.0, vmax=mean)
    plt.matshow(map, cmap="magma")
    plt.show()


if __name__ == "__main__":
    data = load_xy(
        "./2025/02/18/social_preference/ab/social_preference-20250218T142029-xy_position.csv"
    )
    map = np.zeros((1080, 1440))
    for i in range(1, 9):
        arena = np.nan_to_num(
            data.select(f"X_A{i}", f"Y_A{i}").to_numpy().astype(float)
        )
        map += make_map(arena)
        print(map.min(), map.max(), map.mean(), map.std())

    show_map(map)
