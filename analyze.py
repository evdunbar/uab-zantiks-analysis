import io
import itertools
from collections.abc import Callable, Iterable, Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
import polars as pl
import polars.selectors as cs
import seaborn as sns
from scipy import stats


class Analyzer:
    def __init__(
        self,
        analysis_func: Callable[
            [pl.DataFrame, pl.DataFrame, Iterable[str]], pl.DataFrame
        ],
        trial_paths: Sequence[str] | str,
        trial_names: Iterable[str] | str,
        num_zones: int,
        genotypes: Iterable[str] = ("WT", "HET", "HOM"),
        markers: Iterable[str] = ("$\u25cf$", "$\u25d0$", "$\u25cb$"),
    ):
        self.analysis_func = analysis_func
        self.trial_paths = (trial_paths,) if trial_paths is str else trial_paths
        self.trial_names = (trial_names,) if trial_names is str else trial_names
        self.num_zones = num_zones
        self.zone_names = [f"Z{i}" for i in range(1, self.num_zones + 1)]
        self.genotypes = genotypes
        # https://github.com/stipub/stixfonts/blob/master/docs/STIXTwoMath-Regular.pdf
        self.markers = markers

    def _genotypes_to_columns(self, df: pl.DataFrame) -> pl.DataFrame:
        transposed_df = df.transpose(include_header=True).drop_nulls()
        transposed_df.columns = ["Arena", "Genotype"]
        return transposed_df

    def _genotypes_wells_to_arenas(self, df: pl.DataFrame) -> pl.DataFrame:
        wells = df["Well"]
        row = wells.str.head(1).str.replace_many(["A", "B", "C", "D", "E", "F", "G", "H"], ["0", "1", "2", "3", "4", "5", "6", "7"])
        column = wells.str.strip_chars_start("ABCDEFGH")
        arena_nums = row.str.to_integer() * 12 + column.str.to_integer()
        arenas = "A" + arena_nums.cast(str)

        return df.with_columns(arenas).rename({"Well": "Arena"})

    def analyze(self) -> pl.DataFrame:
        analyzed_data = None
        for path, name in zip(self.trial_paths, self.trial_names):
            with open(f"{path}/{name}.csv", "r") as data_csv:
                data_string = data_csv.read()
                if data_string[0] != '"':
                    # get rid of superfluous zantiks metadata
                    data_string = "\n".join(data_string.splitlines()[3:-1])
                data = pl.read_csv(io.StringIO(data_string), infer_schema_length=None)
            with open(f"{path}/{name}_genotypes.csv", "r") as genotypes_csv:
                genotypes = pl.read_csv(genotypes_csv, infer_schema_length=None)

            if genotypes.shape[1] != 2:  # old genotype format
                # flip dataframe
                genotypes = self._genotypes_to_columns(genotypes)
            else:  # hrm output format
                genotypes = self._genotypes_wells_to_arenas(genotypes)

            trial_data = self.analysis_func(data, genotypes, self.zone_names)

            if analyzed_data is None:
                analyzed_data = trial_data.clone()
            else:
                analyzed_data.vstack(trial_data, in_place=True)

        return analyzed_data

    def graph(
        self, df: pl.DataFrame, show_significance: bool = False, save: bool = False
    ):
        data_types = df.columns[2:]
        for i, data_type in enumerate(data_types):
            sns.catplot(
                data=df,
                kind="bar",  # plotting bar data
                x="Genotype",
                y=data_type,
                width=0.7,  # magic number to match prism
                linewidth=1.5,  # magic number to match prism
                edgecolor=(0, 0, 0, 1),  # magic number to match prism
                facecolor=(0, 0, 0, 0.2),  # magic number to match prism
                errorbar=None,
                order=self.genotypes,
                height=5,  # inches
                aspect=1,  # aspect = width / height
            )
            for genotype, marker in zip(self.genotypes, self.markers):
                sns.swarmplot(
                    data=df.filter(pl.col("Genotype") == genotype),
                    x="Genotype",
                    y=data_type,
                    s=7,  # marker size
                    c="black",  # marker color
                    marker=marker,
                )

            # put one minor tick in between each major y-axis tick
            # n = 1 + number of minor ticks
            plt.gca().yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(n=2))

            # make axis lines and tick marks thicker
            for spine in plt.gca().spines.values():
                spine.set_linewidth(1.5)
            plt.gca().tick_params(axis="both", which="both", width=1.5)

            if save:
                plt.savefig(
                    "figures/" + self.trial_paths[-1].translate({47: 95}) + f"_{i}.png",
                    dpi=200,
                    format="png",
                )
            plt.show()

    def graph_all(
        self, df: pl.DataFrame, show_significance: bool = False, save: bool = False
    ):
        tetragrams = [
            "".join(tetragram) for tetragram in itertools.product(("L", "R"), repeat=4)
        ]
        data_no_arenas = df[1:]
        df = data_no_arenas.unpivot(
            index="Genotype",
            on=tetragrams,
            variable_name="Combination",
            value_name="Value",
        )
        print(df)

        for genotype, marker in zip(self.genotypes, self.markers):
            sns.catplot(
                data=df.filter(pl.col("Genotype") == genotype),
                kind="bar",  # plotting bar data
                x="Combination",
                y="Value",
                width=0.7,  # magic number to match prism
                linewidth=1.5,  # magic number to match prism
                edgecolor=(0, 0, 0, 1),  # magic number to match prism
                facecolor=(0, 0, 0, 0.2),  # magic number to match prism
                errorbar=None,
                height=5,  # inches
                aspect=4,  # aspect = width / height
            )
            sns.swarmplot(
                data=df.filter(pl.col("Genotype") == genotype),
                x="Combination",
                y="Value",
                s=7,  # marker size
                c="black",  # marker color
                marker=marker,
            )

            # put one minor tick in between each major y-axis tick
            # n = 1 + number of minor ticks
            plt.gca().yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(n=2))

            # make axis lines and tick marks thicker
            for spine in plt.gca().spines.values():
                spine.set_linewidth(1.5)
            plt.gca().tick_params(axis="both", which="both", width=1.5)

        if save:
            plt.savefig(
                "figures/" + self.trial_paths[-1].translate({47: 95}) + f"_{i}.png",
                dpi=200,
                format="png",
            )
        plt.show()

    def stat(self, df: pl.DataFrame, to_compare: str) -> dict | None:
        all_data = [
            df.filter(pl.col("Genotype") == genotype)[to_compare]
            for genotype in self.genotypes
        ]
        anova = stats.f_oneway(*all_data)

        # return nothing if not significant
        if anova.pvalue >= 0.05:
            return None

        # implied else
        combinations = itertools.combinations(self.genotypes, 2)
        tukey = {}
        for combination, idx in zip(combinations, itertools.combinations(range(3), 2)):
            pval = stats.tukey_hsd(all_data[idx[0]], all_data[idx[1]]).pvalue
            tukey[combination] = "*" if pval < 0.05 else "ns"
        return tukey


def mirror_biting(
    data: pl.DataFrame, genotypes: pl.DataFrame, zone_names: Iterable[str]
) -> pl.DataFrame:
    # filter data to only include columns from arenas that contained fish
    # sum across time to get totals for whole experiment
    filled_arenas = genotypes.get_column("Arena").to_list()
    # need list comprehension to not have "A10" match substring "A1"
    valid_data = data.select(
        cs.contains(*[arena_name + "." for arena_name in filled_arenas])
    ).sum()

    # add time columns to data
    for zone_name in zone_names:
        zone_data = valid_data.select(
            cs.starts_with("T") & cs.ends_with(zone_name)
        ).transpose()
        zone_data.columns = [f"Time in Zone {zone_name[1:]} (s)"]
        genotypes.hstack(zone_data, in_place=True)

    # add distance columns to data
    for zone_name in zone_names:
        zone_data = valid_data.select(
            cs.starts_with("D") & cs.ends_with(zone_name)
        ).transpose()
        zone_data.columns = [f"Distance in Zone {zone_name[1:]}"]
        genotypes.hstack(zone_data, in_place=True)
    total_distance_data = genotypes[:, -3:].sum_horizontal()
    total_distance_data = total_distance_data.alias("Distance Total")
    genotypes.hstack([total_distance_data], in_place=True)

    return genotypes


def y_maze_spontaneous_alternation(
    data: pl.DataFrame, genotypes: pl.DataFrame, _
) -> pl.DataFrame:
    # remove unecessary columns
    data = data.drop("", "_duplicated_0")

    # get alternation percentages
    alternation_percentages = []
    for arena_name in genotypes["Arena"]:
        arena_id = int(arena_name[1:])
        arena_rows = data.filter(
            (pl.col("ARENA") == arena_id)
            & (pl.col("ACTION") == "Enter_Zone")
            & (pl.col("ZONE") != 4)  # zone 4 is center of maze
        )

        num_triads = arena_rows["ACTION"].count() - 2
        spontaneous_alterations = 0
        for i in range(num_triads):
            triad = arena_rows["ZONE"].slice(offset=i, length=3)
            if triad.n_unique() == 3:
                spontaneous_alterations += 1

        alternation_percentages.append(spontaneous_alterations / num_triads * 100)

    trial_data = genotypes.with_columns(
        pl.Series(name="Alternation (%)", values=alternation_percentages)
    )
    return trial_data


def y_maze_free_movement_pattern(
    data: pl.DataFrame, genotypes: pl.DataFrame, _
) -> pl.DataFrame:
    # remove unecessary columns
    data = data.drop("", "_duplicated_0")

    # set up data structures
    possible_tetragrams = [
        "".join(tetragram) for tetragram in itertools.product(("L", "R"), repeat=4)
    ]
    tetragram_percentages = {
        tetragram: [0.0] * genotypes["Arena"].count()
        for tetragram in possible_tetragrams
    }
    arm_sequence_to_turn_direction = {
        (1, 2): "L",
        (1, 3): "R",
        (2, 1): "R",
        (2, 3): "L",
        (3, 1): "L",
        (3, 2): "R",
    }

    # for i, arena_name in enumerate(genotypes["Arena"]):
    for i, arena_id in enumerate(data["ARENA"].unique().sort()):
        # arena_id = int(arena_name[1:])
        arena_rows = data.filter(
            (pl.col("ARENA") == arena_id)
            & (pl.col("ACTION") == "Enter_Zone")
            & (pl.col("ZONE") != 4)  # zone 4 is center of maze
        )

        # create turn sequence
        previous_arms = arena_rows["ZONE"].shift(1)[1:]
        arms = arena_rows["ZONE"][1:]
        turns = []
        for arm_sequence in zip(previous_arms, arms):
            if arm_sequence[0] == arm_sequence[1]:
                continue
            turns.append(arm_sequence_to_turn_direction[arm_sequence])

        # tabulate tetragrams
        num_tetragrams = len(turns) - 3
        for j in range(num_tetragrams):
            tetragram = "".join(turns[j : j + 4])
            tetragram_percentages[tetragram][i] += 1 / num_tetragrams

    tetragram_percentages = pl.DataFrame(tetragram_percentages)
    trial_data = genotypes.with_columns(tetragram_percentages)
    print(trial_data)
    return trial_data


if __name__ == "__main__":
    # current_paths = (
    #     "2024/09/18/mirror_biting",
    #     "2024/09/19/mirror_biting",
    #     "2024/09/19/mirror_biting",
    # )
    # current_names = (
    #     "SLCb_b3m3_f2s_3wpf_201241",
    #     "SLCb_b3m3_f2s_3wpf_170229",
    #     "SLCb_b3m3_f2s_3wpf_174653",
    # )
    # current_paths = ("2024/09/18/y_maze", "2024/09/18/y_maze", "2024/09/18/y_maze")
    # current_names = (
    #     "SLCb_b3m3_f2s_6dpf_145620",
    #     "SLCb_b3m3_f2s_6dpf_204832",
    #     "SLCb_b3m3_f2s_6dpf_222400",
    # )
    # current_paths = (
    #     "2024/09/23/y_maze",
    #     "2024/09/23/y_maze",
    #     "2024/09/23/y_maze",
    #     "2024/09/23/y_maze",
    # )
    # current_names = (
    #     "SLCa_e1m1_f2s_6dpf_160629",
    #     "SLCa_e1m1_f2s_6dpf_180730",
    #     "SLCa_e1m1_f2s_6dpf_185241",
    #     "SLCa_e1m1_f2s_6dpf_203207",
    # )
    current_paths = ("2024/11/01/ymaze_15/a",)
    current_names = ("ymaze_15-20241101T164208",)

    analysis = Analyzer(y_maze_free_movement_pattern, current_paths, current_names, 4)
    analyzed_data = analysis.analyze()
    # print(analyzed_data)
    analysis.graph_all(analyzed_data, save=False)
    # analyzed_data.write_csv("y_maze_slca.csv")
