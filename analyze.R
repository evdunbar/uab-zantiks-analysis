#!/usr/bin/env Rscript

########
# GOAL #
########

# To turn microtracker/zantiks data into data you can copy/paste into GraphPad Prism


################
# INSTRUCTIONS #
################

# - Read the info about each parameter below and fill in according to your desires
# - If you need more info, check the README.md file
# - All relative paths are relative to the working directory


##############
# PARAMETERS #
##############

# config_file:
# - type:
#   - string
# - choices:
#   1. ""
#     - no config file specified - use the parameters defined in this file for analysis
#   2. "relative/path/to/file.toml"
#     - use the configuration options contained in the file located at this path
#     - all other configuration options will be ignored
CONFIG_FILE <- ""

# data_file_prefix:
# - type:
#   - string
# - choices:
#   1. ""
#     - nothing is prepended to the paths contained in data_files
#   2. "relative/path"
#     - prepend this string to each string in data_files
#     - for example, if you had files y_maze/a/file.csv and y_maze/b/file.csv, then you
#       could set this to "y_maze/" and then data_files would contain "a/file.csv" and
#       "b/file.csv"
DATA_FILE_PREFIX <- ""
# data_files:
# - type:
#   - vector of strings
# - one choice:
#   1. at least one "relative/path" separated by commas
#   - should contain the raw output of a zantiks or microtracker run
#   - pay attention to the order and use the same one for all other "vector of strings"
DATA_FILES <- c()
# assay_names:
# - type:
#   - vector of strings
# - choices:
#   1. developmental delay
#   2. light/dark preference
#   3. light/dark transition
#   4. microtracker
#   5. mirror biting
#   6. sleep
#   7. social preference
#   8. startle response/pre-pulse inhibition
#   9. total distance
#  10. y-maze 15
#  11. y-maze 4
ASSAY_NAMES <- c()
# wildtype_file:
# - type:
#   - string
# - choices:
#   1. ""
#   2. one "relative/path"
#   - should contain a config file for ALL wild type data for that assay
WILDTYPE_FILE <- ""

# genotyping_file_prefix:
# - type:
#   - string
# - choices:
#   1. ""
#   2. "relative/path"
GENOTYPING_FILE_PREFIX <- ""
# genotyping_files:
# - type:
#   - vector of strings
# - one choice:
#   1. at least one "relative/path" separated by commas
#   - should contain a "Well" and "Cluster" column
GENOTYPING_FILES <- c()
# counting_directions:
# - type:
#   - vector of strings
# - choices:
#   1. "down"
#     - the numbering of fish in my data file is like counting down, then right
#       on the genotyping plate
#     - usually used for zantiks plates
#   2. "across"
#     - the numbering of fish in my data file is like counting right, then down
#       on the genotyping plate
#     - usually only used for microtracker plates
COUNTING_DIRECTIONS <- c()

# fish_used_prefix:
# - type:
#   - string
# - choices:
#   1. ""
#   2. "relative/path"
FISH_USED_PREFIX <- ""
# fish_used:
# - type:
#   - vector of strings
# - one choice:
#   1. at least one "relative/path" separated by commas
FISH_USED <- c()

# output_file:
# - type:
#   - string
# - choices:
#   1. ""
#     - the output will be saved to "output.csv"
#   2. "relative/path"
#     - the output will be saved to this file
OUTPUT_FILE <- ""


###########################
# FUNCTIONS AND LIBRARIES #
###########################

library(tidyverse)
library(dplyr)
library(readr)
library(readxl)
if (CONFIG_FILE != "") {
  library(configr)

  USING_CONFIG <- TRUE

  # uses global CONFIG_FILE variable
  load_config <- function() {
    config <- read.config(CONFIG_FILE)

    # use <<- instead of <- to assign to the global variables
    DATA_FILE_PREFIX <<- config$data$files_prefix
    DATA_FILES <<- config$data$files
    ASSAY_NAMES <<- config$data$assay_names
    WILDTYPE_FILE <<- config$data$wildtype_file

    GENOTYPING_FILE_PREFIX <<- config$genotypes$files_prefix
    GENOTYPING_FILES <<- config$genotypes$files

    FISH_USED_PREFIX <<- config$fish_used$files_prefix
    FISH_USED <<- config$fish_used$files

    COUNTING_DIRECTIONS <<- config$genotypes$counting_directions

    maybe_output_file <- config$output$file
    if (length(maybe_output_file) != 0) OUTPUT_FILE <<- maybe_output_file else OUTPUT_FILE <<- ""
  }
} else {
  USING_CONFIG <- FALSE
}

# uses global variables, so no direct inputs
validate_inputs <- function() {
  # these values are required
  required_values <- list(DATA_FILES, GENOTYPING_FILES, ASSAY_NAMES, FISH_USED, COUNTING_DIRECTIONS)
  for (required_value in required_values) {
    if (length(required_value) < 1) {
      print(required_value)
      stop("DATA_FILES, GENOTYPING_FILES, ASSAY_NAMES, FISH_USED, and COUNTING_DIRECTIONS must not be empty")
    }
  }

  # all required values need to have the same length to match up
  if (length(unique(map(required_values, length))) != 1) {
    for (required_value in required_values) {
      print("lengths:")
      print(length(required_value))
    }
    stop("DATA_FILES, GENOTYPING_FILES, ASSAY_NAMES, FISH_USED, and COUNTING_DIRECTIONS must be the same length")
  }

  # we have to know how to analyze an assay
  valid_assay_names <- c(
    "developmental delay",
    "light/dark preference",
    "light/dark preference 6dpf",
    "light/dark preference 3wpf",
    "light/dark transition",
    "microtracker",
    "mirror biting",
    "sleep",
    "social preference",
    "startle response/pre-pulse inhibition",
    "total distance",
    "y-maze 15",
    "y-maze 4"
  )
  if (!all(ASSAY_NAMES %in% valid_assay_names)) {
    print("assay names:")
    print(ASSAY_NAMES)
    stop("ASSAY_NAMES must all match one of the options in the list")
  }

  valid_counting_directions <- c("down", "across")
  if (!all(COUNTING_DIRECTIONS %in% valid_counting_directions)) {
    print("counting directions:")
    print(COUNTING_DIRECTIONS)
    stop("COUNTING_DIRECTIONS must all match \"down\" or \"across\"")
  }
}

# turn full genotype file into a table that can be merged with the data
process_genotypes <- function(genotyping_file, fish_used_file, counting_direction) {
  # read in data about which fish were used
  fish_used_data <- read.table(fish_used_file) %>%
    as.matrix()

  # create matrix to match well labels in genotype data
  label_matrix <- matrix(
    outer(LETTERS[1:8],
      sprintf("%02d", 1:12),
      FUN = paste0
    ),
    nrow = 8,
    ncol = 12
  )

  # get which wells were used
  wells_used <- label_matrix[fish_used_data == "x" | fish_used_data == "X"]

  # read in data
  genotype_data <- read_csv(genotyping_file) %>%
    select(genotyping_well = Well, genotype = Cluster, clutch = Clutch) %>%
    mutate(row = str_extract(genotyping_well, "[A-H]"), column = as.integer(str_extract(genotyping_well, "[0-9]+"))) %>%
    mutate(genotyping_well = paste0(row, sprintf("%02d", column))) %>%
    filter(genotyping_well %in% wells_used)

  # sort by row or column
  if (counting_direction == "across") {
    genotype_data <- genotype_data %>%
      arrange(row, column)
  } else {
    genotype_data <- genotype_data %>%
      arrange(column, row)
  }

  # add row ids to join on
  genotype_data <- genotype_data %>%
    mutate(row_id = row_number())

  genotype_data
}


#####################################
# ASSAY SPECIFIC ANALYSIS FUNCTIONS #
#####################################

developmental_delay_analysis <- function(data_file, genotypes) {
  lines <- readLines(data_file)
  csv_text <- lines[4:(length(lines) - 1)]
  data <- read_csv(I(csv_text), col_types = "didddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd")

  processed_data <- data %>%
    pivot_longer(
      cols = A1:A96,
      names_to = "ARENA",
      names_pattern = "A([0-9]+)",
      values_to = "PIXEL_DIFFERENCE"
    ) %>%
    mutate(ARENA = as.integer(ARENA)) %>%
    relocate(ARENA) %>%
    arrange(ARENA) %>%
    group_by(ARENA) %>%
    summarize(
      `developmental delay: pixel difference` = sum(PIXEL_DIFFERENCE),
    )

  stupid_genotypes <- genotypes %>%
    mutate(
      row_as_num = match(row, LETTERS),
      row_id = case_when(
        row_as_num %in% c(1, 2, 3, 4) & column <= 6 ~ (row_as_num - 1) * 6 + column,
        row_as_num %in% c(1, 2, 3, 4) & column > 6 ~ (row_as_num - 1) * 6 + (column - 6) + 24,
        row_as_num %in% c(5, 6, 7, 8) & column <= 6 ~ (row_as_num - 5) * 6 + column + 48,
        row_as_num %in% c(5, 6, 7, 8) & column > 6 ~ (row_as_num - 5) * 6 + (column - 6) + 72,
      ),
    )

  finished_data <- processed_data %>%
    left_join(stupid_genotypes, by = join_by(ARENA == row_id)) %>%
    arrange(clutch, genotype) %>%
    select(clutch, genotype, `developmental delay: pixel difference`)

  finished_data
}

light_dark_preference_analysis <- function(data_file, genotypes) {
  lines <- readLines(data_file)
  csv_text <- lines[4:(length(lines) - 1)]
  data <- read_csv(I(csv_text), col_types = "dicdddddddddddddddddddddddd")

  processed_data <- data %>%
    pivot_longer(
      cols = A1_Z1:A12_Z2,
      names_to = c("ARENA", "ZONE"),
      names_pattern = "A([0-9]+)_Z([1,2])",
      values_to = "VALUE"
    ) %>%
    mutate(ARENA = as.integer(ARENA), ZONE = as.integer(ZONE)) %>%
    relocate(ARENA, ZONE) %>%
    arrange(ARENA, ZONE) %>%
    mutate(ZONE = case_when(ARENA %in% c(5, 6, 7, 8) & ZONE == 1 ~ 3, .default = ZONE)) %>% # swap zone numbers
    mutate(ZONE = case_when(ARENA %in% c(5, 6, 7, 8) & ZONE == 2 ~ 1, .default = ZONE)) %>% # now 1 is always dark
    mutate(ZONE = case_when(ARENA %in% c(5, 6, 7, 8) & ZONE == 3 ~ 2, .default = ZONE)) %>% # and 2 is always light
    group_by(ARENA) %>%
    summarize(
      `l/d preference: total light time` = sum(VALUE[ENDPOINT == "TIME_SPENT_IN_ZONE" & ZONE == 2]),
      `l/d preference: percent distance light` = sum(VALUE[ENDPOINT == "DISTANCE_IN_ZONE" & ZONE == 2]) /
        sum(VALUE[ENDPOINT == "DISTANCE_IN_ZONE"]) * 100
    )

  finished_data <- processed_data %>%
    left_join(genotypes, by = join_by(ARENA == row_id)) %>%
    arrange(clutch, genotype) %>%
    select(clutch, genotype, `l/d preference: total light time`, `l/d preference: percent distance light`)

  finished_data
}

light_dark_transition_analysis <- function(data_file, genotypes) {
  lines <- readLines(data_file)
  csv_text <- lines[4:(length(lines) - 1)]
  data <- read_csv(I(csv_text), col_types = "dcidddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd")

  processed_data <- data %>%
    pivot_longer(
      cols = A1_Z1:A48_Z2,
      names_to = c("ARENA", "ZONE"),
      names_pattern = "A([0-9]+)_Z([1,2])",
      values_to = "DISTANCE"
    ) %>%
    mutate(ARENA = as.integer(ARENA), ZONE = as.integer(ZONE)) %>%
    relocate(ARENA, ZONE) %>%
    arrange(ARENA, ZONE) %>%
    group_by(ARENA) %>%
    summarize(
      `l/d transition: light/dark ratio` = sum(DISTANCE[CONDITION == "BRIGHT"]) /
        sum(DISTANCE[CONDITION == "DARK"]),
      `l/d transition: thigmotaxis` = sum(DISTANCE[ZONE == 1]) / sum(DISTANCE) * 100
    )

  finished_data <- processed_data %>%
    left_join(genotypes, by = join_by(ARENA == row_id)) %>%
    arrange(clutch, genotype) %>%
    select(clutch, genotype, `l/d transition: light/dark ratio`, `l/d transition: thigmotaxis`)

  finished_data
}

microtracker_analysis <- function(data_file, genotypes) {
  if (str_ends(data_file, ".xlsx")) {
    data <- read_xlsx(data_file, skip = 25, n_max = 96)
  } else if (str_ends(data_file, ".csv")) {
    lines <- readLines(data_file)
    starting_line <- which(str_detect(lines, fixed("Well Activity"))) + 1
    csv_text <- lines[starting_line:(starting_line + 96)]
    data <- read_csv(I(csv_text), col_types = "cciiii")
  }

  data <- data %>%
    select(Well, `30`, `60`, `90`, `120`) %>%
    mutate(`microtracker: average locomotor activity` = rowMeans(across(c(`30`, `60`, `90`, `120`)))) %>%
    mutate(row = str_extract(Well, "[A-H]"), col = as.integer(str_extract(Well, "[0-9]+"))) %>%
    arrange(row, col) %>%
    mutate(row_id = row_number())

  finished_data <- data %>%
    full_join(genotypes, by = join_by(row_id)) %>%
    arrange(clutch, genotype) %>%
    select(clutch, genotype, `microtracker: average locomotor activity`)

  finished_data
}

mirror_biting_analysis <- function(data_file, genotypes) {
  lines <- readLines(data_file)
  csv_text <- lines[4:(length(lines) - 1)]
  data <- read_csv(I(csv_text), col_types = "diddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiidddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd")

  processed_data <- data %>%
    pivot_longer(
      cols = D.A1.Z1:T.A20.Z3,
      names_to = c("CATEGORY", "ARENA", "ZONE"),
      names_pattern = "([D,C,T])\\.A([0-9]+)\\.Z([1-3])",
      values_to = "VALUE"
    ) %>%
    mutate(ARENA = as.integer(ARENA), ZONE = as.integer(ZONE)) %>%
    relocate(ARENA, ZONE) %>%
    arrange(ARENA, ZONE) %>%
    group_by(ARENA) %>%
    summarize(
      `mirror biting: percent mirror time` = sum(VALUE[CATEGORY == "T" & ZONE == 1]) / sum(VALUE[CATEGORY == "T"]) * 100
    )

  finished_data <- processed_data %>%
    left_join(genotypes, by = join_by(ARENA == row_id)) %>%
    arrange(clutch, genotype) %>%
    select(
      clutch, genotype,
      `mirror biting: percent mirror time`
    )

  finished_data
}

sleep_analysis <- function(data_file, genotypes) {
  lines <- readLines(data_file)
  csv_text <- lines[4:(length(lines) - 1)]
  data <- read_csv(I(csv_text), col_types = "dcidddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd")

  processed_data <- data %>%
    pivot_longer(
      cols = A1_Z1:A48_Z2,
      names_to = c("ARENA", "ZONE"),
      names_pattern = "A([0-9]+)_Z([1,2])",
      values_to = "DISTANCE"
    ) %>%
    mutate(ARENA = as.integer(ARENA), ZONE = as.integer(ZONE)) %>%
    relocate(ARENA, ZONE) %>%
    arrange(ARENA, ZONE) %>%
    group_by(ARENA) %>%
    summarize(
      `sleep: light/dark ratio` = (sum(DISTANCE[CONDITION == "BRIGHT"]) / 13) /
        (sum(DISTANCE[CONDITION == "DARK"]) / 10),
      `sleep: 1hr distance` = sum(DISTANCE[TIME <= 3600]),
      `sleep: total light distance` = sum(DISTANCE[CONDITION == "BRIGHT"]),
      `sleep: total dark distance` = sum(DISTANCE[CONDITION == "DARK"]),
    )

  finished_data <- processed_data %>%
    left_join(genotypes, by = join_by(ARENA == row_id)) %>%
    arrange(clutch, genotype) %>%
    select(clutch, genotype, 
      `sleep: light/dark ratio`,
      `sleep: 1hr distance`,
      `sleep: total light distance`,
      `sleep: total dark distance`,
    )

  finished_data
}

social_preference_analysis <- function(data_file, genotypes) {
  lines <- readLines(data_file)
  csv_text <- lines[4:(length(lines) - 1)]
  data <- read_csv(I(csv_text), col_types = "didddddddddddddddddddddddddddddddddddddddddddddddddd")

  processed_data <- data %>%
    pivot_longer(
      cols = T.A1.Z1:T.A10.Z5,
      names_to = c("ARENA", "ZONE"),
      names_pattern = "T.A([0-9]+).Z([1-5])",
      values_to = "SECONDS"
    ) %>%
    mutate(ARENA = as.integer(ARENA), ZONE = as.integer(ZONE)) %>%
    relocate(ARENA, ZONE) %>%
    arrange(ARENA, ZONE) %>%
    group_by(ARENA) %>%
    summarize(
      `social preference: social preference index` = (sum(SECONDS[ZONE == 1]) +
        0.5 * sum(SECONDS[ZONE == 2]) -
        0.5 * sum(SECONDS[ZONE == 4]) -
        sum(SECONDS[ZONE == 5])) / sum(SECONDS)
    )

  finished_data <- processed_data %>%
    left_join(genotypes, by = join_by(ARENA == row_id)) %>%
    arrange(clutch, genotype) %>%
    select(
      clutch, genotype,
      `social preference: social preference index`
    )

  finished_data
}

startle_response_analysis <- function(data_file, genotypes) {
  lines <- readLines(data_file)
  csv_text <- lines[4:(length(lines) - 1)]
  data <- read_csv(I(csv_text), col_types = "ddcdddddddddddddddddddddddddddddddddddddddddddddddd")

  processed_data <- data %>%
    pivot_longer(cols = A1:A48, names_to = "ARENA", names_pattern = "A([0-9]+)", values_to = "DISTANCE") %>%
    mutate(ARENA = as.integer(ARENA)) %>%
    relocate(ARENA) %>%
    arrange(ARENA) %>%
    group_by(ARENA) %>%
    summarize(
      `startle response: pre-pulse` = sum(DISTANCE[PHASE == "PREPULSE"]),
      `startle response: startle alone` = sum(DISTANCE[PHASE == "STARTLE"])
    )

  finished_data <- processed_data %>%
    left_join(genotypes, by = join_by(ARENA == row_id)) %>%
    arrange(clutch, genotype) %>%
    select(clutch, genotype, `startle response: pre-pulse`, `startle response: startle alone`)

  finished_data
}

total_distance_analysis <- function(data_file, genotypes) {
  lines <- readLines(data_file)
  csv_text <- lines[4:(length(lines) - 1)]
  data <- read_csv(I(csv_text), col_types = "ddicdddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd")

  processed_data <- data %>%
    pivot_longer(
      cols = A1_Z1:A48_Z2,
      names_to = c("ARENA", "ZONE"),
      names_pattern = "A([0-9]+)_Z([1,2])",
      values_to = "DISTANCE"
    ) %>%
    mutate(ARENA = as.integer(ARENA), ZONE = as.integer(ZONE)) %>%
    relocate(ARENA, ZONE) %>%
    arrange(ARENA, ZONE) %>%
    group_by(ARENA) %>%
    summarize(
      `total distance: zone 1 distance` = sum(DISTANCE[ZONE == 1]),
      `total distance: zone 2 distance` = sum(DISTANCE[ZONE == 2]),
      `total distance: cumulative distance` = sum(DISTANCE),
      `total distance: average velocity` = sum(DISTANCE) / 3600
    )

  finished_data <- processed_data %>%
    left_join(genotypes, by = join_by(ARENA == row_id)) %>%
    arrange(clutch, genotype) %>%
    select(
      clutch, genotype,
      `total distance: zone 1 distance`,
      `total distance: zone 2 distance`,
      `total distance: cumulative distance`,
      `total distance: average velocity`
    )

  finished_data
}

y_maze_analysis <- function(data_file, genotypes) {
  # read the file and get rid of the auto generated zantiks lines
  # they would mess up csv parsing
  lines <- readLines(data_file)
  csv_text <- lines[4:(length(lines) - 1)]
  data <- read_csv(I(csv_text), col_types = "dccici")

  # set up function to map zone sequences to turn directions
  zone_sequence_mapping <- function(zone, next_zone) {
    result <- case_when(
      zone == 1 & next_zone == 2 ~ "L",
      zone == 1 & next_zone == 3 ~ "R",
      zone == 2 & next_zone == 1 ~ "R",
      zone == 2 & next_zone == 3 ~ "L",
      zone == 3 & next_zone == 1 ~ "L",
      zone == 3 & next_zone == 2 ~ "R"
    )

    result
  }

  # remove center zone and exit data, not needed for analysis
  cleaned_data <- data %>%
    filter(ZONE != 4 & ACTION != "Exit_Zone")

  prepared_data <- cleaned_data %>%
    group_by(ARENA) %>% # calculate everything per arena
    mutate(next_zone = lead(ZONE)) %>%
    mutate(next_next_zone = lead(next_zone)) %>%
    filter(ZONE != next_zone) %>% # only care about zone changes
    mutate(turn_direction = zone_sequence_mapping(ZONE, next_zone)) %>%
    mutate(next_turn = lead(turn_direction)) %>%
    mutate(next_next_turn = lead(next_turn)) %>%
    mutate(next_next_next_turn = lead(next_next_turn))

  triad_data <- prepared_data %>%
    filter(!is.na(next_next_zone)) %>%
    mutate(triad = paste0(ZONE, next_zone, next_next_zone)) %>%
    count(triad, name = "triad_count") %>%
    summarize(
      `y-maze: spontaneous alternation` = sum(triad_count[triad %in% c("123", "231", "312", "321", "213", "132")]) / sum(triad_count) * 100
    )

  tetragram_data <- prepared_data %>%
    filter(!is.na(next_next_next_turn)) %>%
    mutate(tetragram = paste0(
      turn_direction,
      next_turn,
      next_next_turn,
      next_next_next_turn
    )) %>% # probably a better way to do this...
    count(tetragram, name = "tetragram_count") %>%
    summarize(
      `y-maze: alternation` = sum(tetragram_count[tetragram %in% c("LRLR", "RLRL")]) / sum(tetragram_count) * 100,
      `y-maze: repetition` = sum(tetragram_count[tetragram %in% c("LLLL", "RRRR")]) / sum(tetragram_count) * 100,
      `y-maze: turns` = sum(tetragram_count) + 3,
    )

  processed_data <- full_join(triad_data, tetragram_data)

  # add genotyping
  processed_data <- processed_data %>%
    left_join(genotypes, by = join_by(ARENA == row_id)) %>%
    arrange(clutch, genotype) %>%
    select(
      clutch,
      genotype,
      `y-maze: alternation`,
      `y-maze: repetition`,
      `y-maze: turns`,
      `y-maze: spontaneous alternation`
    )

  processed_data
}


################
# MAIN PROGRAM #
################

if (USING_CONFIG) {
  load_config()
}

# validate user inputs here so we can assume they are right later
validate_inputs()

for (idx in seq_along(DATA_FILES)) {
  # get correct files by corresponding index
  data_file <- paste0(DATA_FILE_PREFIX, DATA_FILES[idx])
  genotyping_file <- paste0(GENOTYPING_FILE_PREFIX, GENOTYPING_FILES[idx])
  fish_used_file <- paste0(FISH_USED_PREFIX, FISH_USED[idx])
  assay_name <- ASSAY_NAMES[idx]
  counting_direction <- COUNTING_DIRECTIONS[idx]

  genotypes <- process_genotypes(genotyping_file, fish_used_file, counting_direction)

  # run analysis
  if (assay_name == "developmental delay") {
    data <- developmental_delay_analysis(data_file, genotypes)
  } else if (str_starts(assay_name, fixed("light/dark preference"))) {
    data <- light_dark_preference_analysis(data_file, genotypes)
  } else if (assay_name == "light/dark transition") {
    data <- light_dark_transition_analysis(data_file, genotypes)
  } else if (assay_name == "microtracker") {
    data <- microtracker_analysis(data_file, genotypes)
  } else if (assay_name == "mirror biting") {
    data <- mirror_biting_analysis(data_file, genotypes)
  } else if (assay_name == "sleep") {
    data <- sleep_analysis(data_file, genotypes)
  } else if (assay_name == "social preference") {
    data <- social_preference_analysis(data_file, genotypes)
  } else if (assay_name == "startle response/pre-pulse inhibition") {
    data <- startle_response_analysis(data_file, genotypes)
  } else if (assay_name == "total distance") {
    data <- total_distance_analysis(data_file, genotypes)
  } else if (assay_name == "y-maze 15" || assay_name == "y-maze 4") {
    data <- y_maze_analysis(data_file, genotypes)
  }

  # combine all data files
  if (idx == 1) {
    all_data <- data
  } else {
    all_data <- bind_rows(all_data, data)
  }
}

# add wildtype data if needed
if (WILDTYPE_FILE != "") {
  CONFIG_FILE <- WILDTYPE_FILE
  load_config()
  validate_inputs()

  for (idx in seq_along(DATA_FILES)) {
    # get correct files by corresponding index
    data_file <- paste0(DATA_FILE_PREFIX, DATA_FILES[idx])
    genotyping_file <- paste0(GENOTYPING_FILE_PREFIX, GENOTYPING_FILES[idx])
    fish_used_file <- paste0(FISH_USED_PREFIX, FISH_USED[idx])
    assay_name <- ASSAY_NAMES[idx]
    counting_direction <- COUNTING_DIRECTIONS[idx]

    genotypes <- process_genotypes(genotyping_file, fish_used_file, counting_direction)

    # run analysis
    if (assay_name == "developmental delay") {
      data <- developmental_delay_analysis(data_file, genotypes)
    } else if (str_starts(assay_name, fixed("light/dark preference"))) {
      data <- light_dark_preference_analysis(data_file, genotypes)
    } else if (assay_name == "light/dark transition") {
      data <- light_dark_transition_analysis(data_file, genotypes)
    } else if (assay_name == "microtracker") {
      data <- microtracker_analysis(data_file, genotypes)
    } else if (assay_name == "mirror biting") {
      data <- mirror_biting_analysis(data_file, genotypes)
    } else if (assay_name == "sleep") {
      data <- sleep_analysis(data_file, genotypes)
    } else if (assay_name == "social preference") {
      data <- social_preference_analysis(data_file, genotypes)
    } else if (assay_name == "startle response/pre-pulse inhibition") {
      data <- startle_response_analysis(data_file, genotypes)
    } else if (assay_name == "total distance") {
      data <- total_distance_analysis(data_file, genotypes)
    } else if (assay_name == "y-maze 15" || assay_name == "y-maze 4") {
      data <- y_maze_analysis(data_file, genotypes)
    }

    # combine all data files
    if (idx == 1) {
      wildtype_data <- data
    } else {
      wildtype_data <- bind_rows(wildtype_data, data)
    }
  }

  # combine all wildtype data with no repeats
  all_data_no_wildtype <- filter(all_data, genotype != "WT")
  wildtype_data <- filter(wildtype_data, genotype == "WT")
  all_data <- bind_rows(all_data_no_wildtype, wildtype_data)
}

all_data <- all_data %>%
  arrange(clutch, desc(genotype)) %>%
  filter(genotype %in% c("HET", "HOM", "WT"))
if (OUTPUT_FILE == "") {
  write_csv(all_data, "output.csv")
} else {
  write_csv(all_data, OUTPUT_FILE)
}
print(all_data)
