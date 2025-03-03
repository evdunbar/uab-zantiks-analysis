#!/usr/bin/env uv run
import math

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from numpy.typing import NDArray
from scipy import signal


def load_xy(filepath: str) -> pl.DataFrame:
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


def make_map(xy_data: NDArray) -> NDArray:
    width = 12776
    height = 8540
    heatmap = np.zeros((height, width))

    for x, y in xy_data:
        if x <= 0.0 and y <= 0.0:
            continue
        scaled_x = math.floor(x * 100)
        scaled_y = math.floor(y * 100)
        heatmap[scaled_y, scaled_x] += 1

    return heatmap


def show_map(heatmap: NDArray) -> None:
    # mask = heatmap > 0.0
    # mean = heatmap[mask].mean()
    # std = heatmap[mask].std()
    # heatmap = ndimage.gaussian_filter(heatmap, 50)
    heatmap = signal.decimate(heatmap, q=10, axis=0, n=3)
    heatmap = signal.decimate(heatmap, q=10, axis=1, n=3)
    heatmap = signal.decimate(heatmap, q=10, axis=0, n=3)
    heatmap = signal.decimate(heatmap, q=10, axis=1, n=3)

    plt.figure(dpi=100)
    # plt.imshow(heatmap, cmap="magma", vmin=0.0, vmax=mean)
    plt.matshow(heatmap, cmap="magma")
    plt.show()


if __name__ == "__main__":
    import glob

    filenames = glob.glob("*/*/*/social_preference/*/*xy_position.csv")
    for filename in filenames:
        data = load_xy(filename)
        heatmap = np.zeros((8540, 12776))
        for i in range(1, 11):
            arena = np.nan_to_num(
                data.select(f"X_A{i}", f"Y_A{i}").to_numpy().astype(float)
            )
            heatmap += make_map(arena)
            print(heatmap.min(), heatmap.max(), heatmap.mean(), heatmap.std())

        show_map(heatmap)
