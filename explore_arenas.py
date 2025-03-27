import marimo

__generated_with = "0.11.22"
app = marimo.App(width="medium", css_file="data/marimo.css")


@app.cell
def _():
    import data_utils
    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    from PIL import Image
    from numpy.typing import NDArray
    return Image, NDArray, data_utils, mo, np, os, plt


@app.cell
def _(Image, NDArray, np):
    def load_arena(filepath: str) -> NDArray:
        with open(filepath, "rb") as f:
            image = np.asarray(Image.open(f))

        # images have weird signature in top left
        # no color channels have intensity less than 10
        cleaned_image = np.where(image < 8, 0, image)

        return cleaned_image
    return (load_arena,)


@app.cell
def _(load_arena, os):
    arena_files = sorted(["data/arenas/" + filename for filename in filter(lambda s: s.endswith(".bmp"), os.listdir("data/arenas"))])
    arena_images = {filename: load_arena(filename) for filename in arena_files}
    return arena_files, arena_images


@app.cell
def _(arena_files, mo):
    selected_file = mo.ui.radio(arena_files, value=arena_files[0], label="Arena file to examine")
    selected_file
    return (selected_file,)


@app.cell
def _(arena_images, selected_file):
    selected_filename = selected_file.value
    selected_arena_image = arena_images[selected_filename]
    return selected_arena_image, selected_filename


@app.cell
def _(plt, selected_arena_image, selected_filename):
    plt.imshow(selected_arena_image)
    plt.title(selected_filename)
    plt.tight_layout()
    plt.axis('off')
    plt.gca()
    return


@app.cell
def _(NDArray, np):
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
    return (arena_color,)


@app.cell
def _(mo, np, selected_arena_image):
    all_unique_values = np.unique(selected_arena_image.reshape((-1, 3)), axis=0)
    # print(all_unique_values)
    num_arenas = len(all_unique_values) - 1 # remove black background
    unique_val = mo.ui.slider(start=1, stop=num_arenas, label="Arena to view")
    unique_val
    return all_unique_values, num_arenas, unique_val


@app.cell
def _(arena_color, np, plt, selected_arena_image, unique_val):
    picked_value = arena_color(unique_val.value)
    bool_image = np.all(selected_arena_image == picked_value, axis=-1)
    plt.imshow(bool_image, cmap="Pastel1")
    plt.title(f"Arena {unique_val.value} {picked_value}")
    plt.tight_layout()
    plt.axis("off")
    plt.gca()
    return bool_image, picked_value


if __name__ == "__main__":
    app.run()
