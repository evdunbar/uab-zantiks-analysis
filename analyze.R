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
CONFIG_FILE <- "configs/jip3_test/3wpf/social_preference.toml"

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
DATA_FILES <- c(
)
# assay_names:
# - type:
#   - vector of strings
# - choices:
#   1. light/dark preference
#   2. light/dark transition
#   3. microtracker
#   4. mirror biting
#   5. social preference
#   6. startle response/pre-pulse inhibition
#   7. total distance
#   8. y-maze 15
#   9. y-maze 4
ASSAY_NAMES <- c(
)

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
GENOTYPING_FILES <- c(
)
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
COUNTING_DIRECTIONS <- c(
)


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
FISH_USED <- c(
)

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

    GENOTYPING_FILE_PREFIX <<- config$genotypes$files_prefix
    GENOTYPING_FILES <<- config$genotypes$files

    ASSAY_NAMES <<- config$data$assay_names

    FISH_USED_PREFIX <<- config$fish_used$files_prefix
    FISH_USED <<- config$fish_used$files

    COUNTING_DIRECTIONS <<- config$genotypes$counting_directions
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
  valid_assay_names <- c("light/dark preference",
                         "light/dark transition",
                         "microtracker",
                         "mirror biting",
                         "social preference",
                         "startle response/pre-pulse inhibition",
                         "total distance",
                         "y-maze 15",
                         "y-maze 4")
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
    select(genotyping_well = Well, genotype = Cluster) %>%
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

  return(genotype_data)
}


#####################################
# ASSAY SPECIFIC ANALYSIS FUNCTIONS #
#####################################

light_dark_preference_analysis <- function(data_file, genotypes) {
  lines <- readLines(data_file)
  csv_text <- lines[4:(length(lines) - 1)]
  data <- read_csv(I(csv_text), col_types = "dicdddddddddddddddddddddddd")

  processed_data <- data %>%
    pivot_longer(cols = A1_Z1:A12_Z2,
                 names_to = c("ARENA", "ZONE"),
                 names_pattern = "A([0-9]+)_Z([1,2])",
                 values_to = "VALUE") %>%
    mutate(ARENA = as.integer(ARENA), ZONE = as.integer(ZONE)) %>%
    relocate(ARENA, ZONE) %>%
    arrange(ARENA, ZONE) %>%
    mutate(ZONE = case_when(ARENA %in% c(5, 6, 7, 8) & ZONE == 1 ~ 3, .default = ZONE)) %>% # swap zone numbers
    mutate(ZONE = case_when(ARENA %in% c(5, 6, 7, 8) & ZONE == 2 ~ 1, .default = ZONE)) %>% # now 1 is always dark
    mutate(ZONE = case_when(ARENA %in% c(5, 6, 7, 8) & ZONE == 3 ~ 2, .default = ZONE)) %>% # and 2 is always light
    group_by(ARENA) %>%
    summarize(`l/d preference: total light time` = sum(VALUE[ENDPOINT == "TIME_SPENT_IN_ZONE" & ZONE == 2]),
              `l/d preference: percent distance light` = sum(VALUE[ENDPOINT == "DISTANCE_IN_ZONE" & ZONE == 2]) /
                sum(VALUE[ENDPOINT == "DISTANCE_IN_ZONE"]))

  finished_data <- processed_data %>%
    left_join(genotypes, by = join_by(ARENA == row_id)) %>%
    relocate(genotype) %>%
    arrange(genotype) %>%
    select(genotype, `l/d preference: total light time`, `l/d preference: percent distance light`)

  return(finished_data)
}

light_dark_transition_analysis <- function(data_file, genotypes) {
  lines <- readLines(data_file)
  csv_text <- lines[4:(length(lines) - 1)]
  data <- read_csv(I(csv_text), col_types = "dcidddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd")

  processed_data <- data %>%
    pivot_longer(cols = A1_Z1:A48_Z2,
                 names_to = c("ARENA", "ZONE"),
                 names_pattern = "A([0-9]+)_Z([1,2])",
                 values_to = "DISTANCE") %>%
    mutate(ARENA = as.integer(ARENA), ZONE = as.integer(ZONE)) %>%
    relocate(ARENA, ZONE) %>%
    arrange(ARENA, ZONE) %>%
    group_by(ARENA) %>%
    summarize(`l/d transition: light/dark ratio` = sum(DISTANCE[CONDITION == "BRIGHT"]) /
                sum(DISTANCE[CONDITION == "DARK"]),
              `l/d transition: thigmotaxis` = sum(DISTANCE[ZONE == 1]) / sum(DISTANCE))

  finished_data <- processed_data %>%
    left_join(genotypes, by = join_by(ARENA == row_id)) %>%
    relocate(genotype) %>%
    arrange(genotype) %>%
    select(genotype, `l/d transition: light/dark ratio`, `l/d transition: thigmotaxis`)

  return(finished_data)
}

microtracker_analysis <- function(data_file, genotypes) {
  if (str_ends(data_file, ".xlsx")) {
    data <- read_xlsx(data_file, skip = 25, n_max = 96)
  } else if (str_ends(data_file, ".csv")) {
    lines <- readLines(data_file)
    csv_text <- lines[26:122]
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
    relocate(genotype) %>%
    select(genotype, `microtracker: average locomotor activity`) %>%
    arrange(genotype)

  return(finished_data)
}

mirror_biting_analysis <- function(data_file, genotypes) {
  lines <- readLines(data_file)
  csv_text <- lines[4:(length(lines) - 1)]
  data <- read_csv(I(csv_text), col_types = "diddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiidddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd")

  processed_data <- data %>%
    pivot_longer(cols = D.A1.Z1:T.A20.Z3,
                 names_to = c("CATEGORY", "ARENA", "ZONE"),
                 names_pattern = "([D,C,T])\\.A([0-9]+)\\.Z([1-3])",
                 values_to = "VALUE") %>%
    mutate(ARENA = as.integer(ARENA), ZONE = as.integer(ZONE)) %>%
    relocate(ARENA, ZONE) %>%
    arrange(ARENA, ZONE) %>%
    group_by(ARENA) %>%
    summarize(
      `mirror biting: percent mirror time` = sum(VALUE[CATEGORY == "T" & ZONE == 1]) / sum(VALUE[CATEGORY == "T"])
    )

  finished_data <- processed_data %>%
    left_join(genotypes, by = join_by(ARENA == row_id)) %>%
    relocate(genotype) %>%
    arrange(genotype) %>%
    select(genotype,
           `mirror biting: percent mirror time`)

  return(finished_data)
}

social_preference_analysis <- function(data_file, genotypes) {
  lines <- readLines(data_file)
  csv_text <- lines[4:(length(lines) - 1)]
  data <- read_csv(I(csv_text), col_types = "didddddddddddddddddddddddddddddddddddddddddddddddddd")

  processed_data <- data %>%
    pivot_longer(cols = T.A1.Z1:T.A10.Z5,
                 names_to = c("ARENA", "ZONE"),
                 names_pattern = "T.A([0-9]+).Z([1-5])",
                 values_to = "SECONDS") %>%
    mutate(ARENA = as.integer(ARENA), ZONE = as.integer(ZONE)) %>%
    relocate(ARENA, ZONE) %>%
    arrange(ARENA, ZONE) %>%
    group_by(ARENA) %>%
    summarize(
      `social preference: social preference index` = (sum(SECONDS[ZONE == 4]) +
                                                        0.5 * sum(SECONDS[ZONE == 3]) -
                                                        0.5 * sum(SECONDS[ZONE == 2]) -
                                                        sum(SECONDS[ZONE == 1])) / sum(SECONDS)
    )

  finished_data <- processed_data %>%
    left_join(genotypes, by = join_by(ARENA == row_id)) %>%
    relocate(genotype) %>%
    arrange(genotype) %>%
    select(genotype,
           `social preference: social preference index`)

  return(finished_data)
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
    summarize(`startle response: pre-pulse` = sum(DISTANCE[PHASE == "PREPULSE"]), `startle response: startle alone` = sum(DISTANCE[PHASE == "STARTLE"]))

  finished_data <- processed_data %>%
    left_join(genotypes, by = join_by(ARENA == row_id)) %>%
    relocate(genotype) %>%
    arrange(genotype) %>%
    select(genotype, `startle response: pre-pulse`, `startle response: startle alone`)

  return(finished_data)
}

total_distance_analysis <- function(data_file, genotypes) {
  lines <- readLines(data_file)
  csv_text <- lines[4:(length(lines) - 1)]
  data <- read_csv(I(csv_text), col_types = "ddicdddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd")

  processed_data <- data %>%
    pivot_longer(cols = A1_Z1:A48_Z2,
                 names_to = c("ARENA", "ZONE"),
                 names_pattern = "A([0-9]+)_Z([1,2])",
                 values_to = "DISTANCE") %>%
    mutate(ARENA = as.integer(ARENA), ZONE = as.integer(ZONE)) %>%
    relocate(ARENA, ZONE) %>%
    arrange(ARENA, ZONE) %>%
    group_by(ARENA) %>%
    summarize(`total distance: zone 1 distance` = sum(DISTANCE[ZONE == 1]),
              `total distance: zone 2 distance` = sum(DISTANCE[ZONE == 2]),
              `total distance: cumulative distance` = sum(DISTANCE),
              `total distance: average velocity` = sum(DISTANCE) / 3600)

  finished_data <- processed_data %>%
    left_join(genotypes, by = join_by(ARENA == row_id)) %>%
    relocate(genotype) %>%
    arrange(genotype) %>%
    select(genotype,
           `total distance: zone 1 distance`,
           `total distance: zone 2 distance`,
           `total distance: cumulative distance`,
           `total distance: average velocity`)

  return(finished_data)
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

    return(result)
  }

  # remove center zone and exit data, not needed for analysis
  cleaned_data <- data %>%
    filter(ZONE != 4 & ACTION != "Exit_Zone")

  # make tetragrams and process data
  processed_data <- cleaned_data %>%
    group_by(ARENA) %>% # calculate everything per arena
    mutate(next_zone = lead(ZONE)) %>%
    filter(ZONE != next_zone) %>% # only care about zone changes
    mutate(turn_direction = zone_sequence_mapping(ZONE, next_zone)) %>%
    mutate(next_turn = lead(turn_direction)) %>%
    mutate(next_next_turn = lead(next_turn)) %>%
    mutate(next_next_next_turn = lead(next_next_turn)) %>%
    filter(!is.na(next_next_next_turn)) %>%
    mutate(tetragram = paste0(
      turn_direction,
      next_turn,
      next_next_turn,
      next_next_next_turn
    )) %>% # probably a better way to do this...
    count(tetragram, name = "count") %>%
    mutate(percentage = count / sum(count) * 100) %>%
    summarize(
      `y-maze: alternations` = sum(count[tetragram %in% c("LRLR", "RLRL")]) / sum(count),
      `y-maze: repetitions` = sum(count[tetragram %in% c("LLLL", "RRRR")]) / sum(count),
      `y-maze: turns` = sum(count)
    )

  # add genotyping
  processed_data <- processed_data %>%
    left_join(genotypes, by = join_by(ARENA == row_id)) %>%
    relocate(genotype) %>%
    arrange(genotype) %>%
    select(genotype, `y-maze: alternations`, `y-maze: repetitions`, `y-maze: turns`)

  return(processed_data)
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

  if (assay_name == "light/dark preference") {
    data <- light_dark_preference_analysis(data_file, genotypes)
  } else if (assay_name == "light/dark transition") {
    data <- light_dark_transition_analysis(data_file, genotypes)
  } else if (assay_name == "microtracker") {
    data <- microtracker_analysis(data_file, genotypes)
  } else if (assay_name == "mirror biting") {
    data <- mirror_biting_analysis(data_file, genotypes)
  } else if (assay_name == "social preference") {
    data <- social_preference_analysis(data_file, genotypes)
  } else if (assay_name == "startle response/pre-pulse inhibition") {
    data <- startle_response_analysis(data_file, genotypes)
  } else if (assay_name == "total distance") {
    data <- total_distance_analysis(data_file, genotypes)
  } else if (assay_name == "y-maze 15" || assay_name == "y-maze 4") {
    data <- y_maze_analysis(data_file, genotypes)
  }

  if (idx == 1) {
    all_data <- data
  } else {
    all_data <- bind_rows(all_data, data)
  }
}

all_data <- all_data %>%
  arrange(genotype) %>%
  filter(genotype %in% c("HET", "HOM", "WT"))
if (OUTPUT_FILE == "") {
  write_csv(all_data, "output.csv")
} else {
  write_csv(all_data, OUTPUT_FILE)
}
print(all_data)
