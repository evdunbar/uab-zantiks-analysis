import polars as pl

def columns_by_value(df: pl.DataFrame, values: tuple = ("HET", "HOM", "WT")) -> dict:
    transposed_df = df.transpose(include_header=True)
    columns_by_value = {}
    for value in values:
        columns = transposed_df.filter(pl.col("column_0") == value)["column"].to_list()
        columns_by_value[value] = columns

    return columns_by_value

def genotypes_to_columns(df: pl.DataFrame) -> pl.DataFrame:
    transposed_df = df.transpose(include_header=True).drop_nulls()
    transposed_df.columns = ["arena", "genotype"]
    return transposed_df
