import marimo

__generated_with = "0.11.22"
app = marimo.App(width="medium", css_file="data/marimo.css")


@app.cell
def _():
    import math
    import os
    from typing import Iterable

    import cv2
    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import polars as pl
    from numpy.typing import NDArray
    from PIL import Image
    from scipy import ndimage

    import data_utils
    return (
        Image,
        Iterable,
        NDArray,
        cv2,
        data_utils,
        math,
        mo,
        ndimage,
        np,
        os,
        pl,
        plt,
    )


@app.cell
def _(mo):
    assay_types = ("light_dark_preference", "light_dark_transition", "mirror_biting", "social_preference", "startle_response", "ymaze_4", "ymaze_15")
    assay_picker = mo.ui.radio(assay_types)
    assay_picker
    return assay_picker, assay_types


@app.cell
def _(Iterable, NDArray, alphashape, data_utils, pl, plt):
    class MapFinder:
        def __init__(self, data: Iterable[data_utils.ZantiksData]):
            # parameters
            self.data: Iterable[data_utils.ZantiksData] = []
            for datum in data:
                if datum.info.is_xy:
                    self.data.append(datum)
                else:
                    print(f"{datum} is not position data")
            self.assay_type = self.data[0].info.assay_type
            self.arena_data = self._convert_to_arenas()

            # constant
            # self.arena_map = (
            #     self._load_image("data/social_preference_arenas.bmp").sum(-1) > 10
            # )

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

        # def _load_image(self, filepath: str) -> NDArray:
        #     with open(filepath, "rb") as f:
        #         image = np.asarray(Image.open(f))
        #     return image

        def _convert_to_arenas(self) -> list[NDArray]:
            arenas = []
            for arena_idx in range(self.assay_type.total_arenas):
                # combine all data
                for i, datum in enumerate(self.data):
                    new_df = (
                        datum.data.filter(pl.col("ARENA") == arena_idx + 1)
                        .select("X", "Y")
                        .drop_nulls()
                    )
                    if i == 0:
                        df = new_df
                    else:
                        df = df.vstack(new_df)

                # keep only unique points to speed up calculations
                arenas.append(df.unique().to_numpy())

            return arenas

        def scatter(self, separate: bool):
            if separate:
                for datum in self.data:
                    plt.figure(dpi=200)
                    plt.title(datum.info.path)
                    plt.scatter(datum.data["X"], datum.data["Y"], s=0.5)
                    plt.show()
            else:
                plt.figure(dpi=200)
                for datum in self.data:
                    plt.scatter(datum.data["X"], datum.data["Y"], s=0.5)

                plt.show()

        def make_shapes(self):
            self.shapes = []

            for points in self.arena_data:
                shape = alphashape.alphashape(points, 1.0)
                self.shapes.append(shape)

        def show_shapes(self):
            for shape, points in zip(self.shapes, self.arena_data):
                fig, ax = plt.subplots()
                # ax.scatter(points[:, 0], points[:, 1])
                x, y = shape.exterior.xy
                ax.fill(x, y, alpha=0.2)

                ax.set_aspect("equal")
                plt.show()
    return (MapFinder,)


@app.cell
def _(Image, NDArray, np, os):
    def load_arena(filepath: str) -> NDArray:
        with open(filepath, "rb") as f:
            image = np.asarray(Image.open(f))

        # images have weird signature in top left
        # no color channels have intensity less than 10
        cleaned_image = np.where(image < 8, 0, image)

        return cleaned_image

    def arena_color(arena: int) -> NDArray:
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

        arena -= 1 # 1-indexed to 0-indexed
        times_to_darken, color = divmod(arena, num_colors)

        return color_cycle[color] - times_to_darken * darken_by

    def get_arena(asset: NDArray, arena: int) -> NDArray:
        return np.all(asset == arena_color(arena + 1), axis=-1)

    arena_files = sorted(["data/arenas/" + filename for filename in filter(lambda s: s.endswith(".bmp"), os.listdir("data/arenas"))])
    arena_images = {filename: load_arena(filename) for filename in arena_files}
    return arena_color, arena_files, arena_images, get_arena, load_arena


@app.cell
def _(MapFinder, assay_picker, data_utils):
    dl = data_utils.DataLoader().add_by_filter(
        assay_types=(assay_picker.value,), data_type="position"
    )
    data = dl.load_all(use_genotypes=False)
    mf = MapFinder(data)
    mf.scatter(separate=False)

    if assay_picker.value == "light_dark_preference":
        bad_directions = None
    elif assay_picker.value == "light_dark_transition":
        bad_directions = None
    elif assay_picker.value == "mirror_biting":
        bad_directions = ("N", "E", "S", "W")
    elif assay_picker.value == "social_preference":
        bad_directions = ("N", "E", "S", "W")
    elif assay_picker.value == "startle_response":
        bad_directions = None
    elif assay_picker.value == "ymaze_4":
        bad_directions = ("N", "S")
    elif assay_picker.value == "ymaze_15":
        bad_directions = ("N", "S")
    return bad_directions, data, dl, mf


@app.cell
def _(NDArray, math, mf, ndimage, np, plt):
    # size of plate in mm is 127.76 x 85.4 
    def bin_data(data: list[NDArray], scale: float, pad: bool = False) -> NDArray:
        if pad:
            height_padding = int(17 * scale)
            width_padding = int(11 * scale)
            bins = np.zeros((math.ceil(85.4 * scale) + height_padding, math.ceil(127.76 * scale) + width_padding))
            for arena_data in data:
                for x, y in arena_data:
                    bins[math.floor(y * scale) + height_padding // 2, math.floor(x * scale) + width_padding // 2] += 1
        else:
            bins = np.zeros((math.ceil(85.4 * scale), math.ceil(127.76 * scale)))
            for arena_data in data:
                for x, y in arena_data:
                    bins[math.floor(y * scale), math.floor(x * scale)] += 1

        bins = ndimage.gaussian_filter(bins, scale)
        bins = np.clip(bins, a_min=0.75*scale, a_max=scale)
        bins -= np.min(bins)
        return bins

    size = 2.5
    binned = bin_data(mf.arena_data, size)
    binned_padded = bin_data(mf.arena_data, size, pad=True)

    plt.figure(dpi=300)
    plt.imshow(binned, cmap="binary")
    plt.colorbar()
    plt.show()
    return bin_data, binned, binned_padded, size


@app.cell
def _(
    arena_color,
    arena_images,
    assay_picker,
    binned_padded,
    mf,
    ndimage,
    np,
    plt,
):
    # make array of just the ground truth arenas that we have data for
    tracking_arenas = arena_images[f"data/arenas/{assay_picker.value}_arenas.bmp"]
    used_arena_idx = []
    used_arenas = np.zeros((tracking_arenas.shape[0], tracking_arenas.shape[1]))
    for arena_idx in range(len(mf.arena_data)):
        tracking_arena = np.all(tracking_arenas == arena_color(arena_idx + 1), axis=-1)
        zoom_factors = (binned_padded.shape[0] / tracking_arena.shape[0], binned_padded.shape[1] / tracking_arena.shape[1])
        resized_tracking_arena = ndimage.zoom(tracking_arena.astype(float), zoom_factors, order=0).astype(bool)
        if np.sum(binned_padded[resized_tracking_arena]) > 0.0:
            used_arenas += tracking_arena
            used_arena_idx.append(arena_idx)

    # change data to new size
    zoom_factors = (used_arenas.shape[0] / binned_padded.shape[0], used_arenas.shape[1] / binned_padded.shape[1])
    upsized_binned = ndimage.zoom(binned_padded, zoom_factors, order=0)

    plt.figure(dpi=300)
    plt.imshow(used_arenas, alpha=0.5, cmap="binary")
    plt.imshow(upsized_binned, alpha=0.5)
    plt.show()
    return (
        arena_idx,
        resized_tracking_arena,
        tracking_arena,
        tracking_arenas,
        upsized_binned,
        used_arena_idx,
        used_arenas,
        zoom_factors,
    )


@app.cell
def _(
    NDArray,
    assay_picker,
    bad_directions,
    cv2,
    get_arena,
    mf,
    np,
    tracking_arenas,
    used_arena_idx,
):
    # get 8 cardinal points for each arena in binned data and zantiks assets
    def get_cardinals(data: NDArray, leave_out: list[str] | None = None) -> NDArray:
        directions = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]
        if leave_out:
            for direction in leave_out:
                directions.remove(direction)

        cardinals = np.empty((len(directions), 2))
        for i, direction in enumerate(directions):
            if direction == "N":
                idx = np.argmin(data[:, 1])
            elif direction == "NE":
                idx = np.argmin(data[:, 1] - data[:, 0])
            elif direction == "E":
                idx = np.argmin(data[:, 0])
            elif direction == "SE":
                idx = np.argmax(data[:, 1] + data[:, 0])
            elif direction == "S":
                idx = np.argmax(data[:, 1])
            elif direction == "SW":
                idx = np.argmax(data[:, 1] - data[:, 0])
            elif direction == "W":
                idx = np.argmax(data[:, 0])
            elif direction == "NW":
                idx = np.argmin(data[:, 1] + data[:, 0])
            else:
                pass
            cardinals[i] = data[idx]
 
        return cardinals

    for i, j in enumerate(used_arena_idx):
        # points from real data
        arena_points = mf.arena_data[j]
        collected_cardinals = get_cardinals(arena_points, leave_out=bad_directions)
        # points from zantiks assets
        asset_arena_points = np.flip(np.argwhere(get_arena(tracking_arenas, j)))
        asset_cardinals = get_cardinals(asset_arena_points, leave_out=bad_directions)

        if i == 0:
            all_real_cardinals = collected_cardinals.copy()
            all_asset_cardinals = asset_cardinals.copy()
        else:
            all_real_cardinals = np.vstack((all_real_cardinals, collected_cardinals))
            all_asset_cardinals = np.vstack((all_asset_cardinals, asset_cardinals))

    obj_points = [np.column_stack((all_real_cardinals, np.zeros(all_real_cardinals.shape[0]))).astype(np.float32)]
    image_points = [all_asset_cardinals.astype(np.float32)]
    ret, mtx, dist, rvecs, tvecs = cv2.calibrateCamera(obj_points, image_points, tracking_arenas.shape[:2][::-1], None, None)
    np.savez_compressed(f"data/camera_params/{assay_picker.value}", mtx=mtx, dist=dist, rvec=rvecs[0], tvec=tvecs[0])
    print(ret, mtx, dist, rvecs, tvecs)
    return (
        all_asset_cardinals,
        all_real_cardinals,
        arena_points,
        asset_arena_points,
        asset_cardinals,
        collected_cardinals,
        dist,
        get_cardinals,
        i,
        image_points,
        j,
        mtx,
        obj_points,
        ret,
        rvecs,
        tvecs,
    )


@app.cell
def _(all_asset_cardinals, cv2, dist, mtx, np, obj_points, plt, rvecs, tvecs):
    out_points, _ = cv2.projectPoints(obj_points[0], rvecs[0], tvecs[0], mtx, dist)
    trans_points = np.squeeze(out_points)
    plt.figure(dpi=300)
    plt.scatter(all_asset_cardinals[:, 0], all_asset_cardinals[:, 1])
    plt.scatter(trans_points[:, 0], trans_points[:, 1])
    plt.show()
    return out_points, trans_points


if __name__ == "__main__":
    app.run()
