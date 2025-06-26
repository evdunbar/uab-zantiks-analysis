import marimo

__generated_with = "0.11.22"
app = marimo.App(width="medium", css_file="data/marimo.css")


@app.cell
def _():
    import marimo as mo
    import matplotlib.pyplot as plt
    import polars as pl
    import seaborn as sns

    import data_utils
    return data_utils, mo, pl, plt, sns


@app.cell
def _(data_utils):
    sleep_config = data_utils.ConfigParser().parse("configs/jip3_test/6dpf/sleep.toml")
    sleep_files = [data_utils.ZantiksFile(sleep_config.data.files_prefix + path) for path in sleep_config.data.files]
    sleep_data = [data_utils.ZantiksData(sleep_file) for sleep_file in sleep_files]
    return sleep_config, sleep_data, sleep_files


@app.cell
def _(pl, plt, sns):
    def make_plot(df: pl.DataFrame, title: str) -> None:
        df = df.with_columns(((pl.col("TIME") - df[0, "TIME"]).cast(pl.Int32) // 3600).alias("Hour"))
        df = df.group_by("ARENA", "Hour").agg(
            pl.col("CONDITION").get(0).eq("BRIGHT").alias("Light"),
            pl.sum("DISTANCE").alias("1-hour Distance"),
            pl.col("Cluster").get(0).alias("Genotype"),
        )
        df = df.filter(pl.col("Genotype").is_in(("WT", "HET", "HOM")))
        df = df.group_by("Genotype", "Hour").agg(
            pl.col("Light").get(0),
            pl.col("1-hour Distance").mean().alias("Average Distance (mm)"),
        )
        df = df.sort("Genotype", "Hour")
        plt.figure(dpi=200)
        sns.barplot(df, x="Hour", y="Average Distance (mm)", hue="Genotype", hue_order=("WT", "HET", "HOM"))
        plt.title(title)
        plt.show()
    return (make_plot,)


@app.cell
def _(make_plot, sleep_data):
    for sleep_datum in sleep_data:
        make_plot(sleep_datum.data, sleep_datum.info.filename)
    return (sleep_datum,)


app._unparsable_cell(
    r"""
    def make_combined_plot(data: list[data_utils.ZantiksData]) -> None:
    
    """,
    name="_"
)


if __name__ == "__main__":
    app.run()
