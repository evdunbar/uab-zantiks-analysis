#!/usr/bin/env Rscript

########
# GOAL #
########
# This script is to process microtracker/zantiks data along with genotyping data (from HRM)
# and output a file to copy and paste into Prism


##############
# PARAMETERS #
##############

# If you would prefer to use a configuration file,
# please input that here
# CONFIG_FILE <- "configs/jip3_test/ymaze_15.toml"
CONFIG_FILE <- ""

# Input the data file paths
# DATA_FILE_PREFIX will be prepended to every string in DATA_FILES
DATA_FILE_PREFIX <- ""
DATA_FILES <- c(
)

# Input the genotyping file paths
# Each genotyping file corresponds to one data file
# GENOTYPING_FILE_PREFIX will be prepended to every string in GENOTYPING_FILES
GENOTYPING_FILE_PREFIX <- ""
GENOTYPING_FILES <- c(
)
# COUNTING_DIRECTIONS can be either across or down
COUNTING_DIRECTIONS <- c(
)

# Choose assay types from:
#   - light/dark preference
#   - light/dark transition
#   - microtracker
#   - mirror biting
#   - social preference
#   - startle response/pre-pulse inhibition
#   - y-maze 15
#   - y-maze 4
#
# Each row corresponds to a data file
ASSAY_NAMES <- c(
)

# Input the fish used paths
FISH_USED_PREFIX <- ""
FISH_USED <- c(
)

# Where should the output to be saved to?
# Name should end with .csv
# Leave empty for default of "output.csv"
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

    DATA_FILE_PREFIX <- config$data$files_prefix
    DATA_FILES <- config$data$files

    GENOTYPING_FILE_PREFIX <- config$genotypes$files_prefix
    GENOTYPING_FILES <- config$genotypes$files

    ASSAY_NAMES <- config$data$assay_names

    FISH_USED_PREFIX <- config$fish_used$files_prefix
    FISH_USED <- config$fish_used$files

    COUNTING_DIRECTIONS <- config$genotypes$counting_directions
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
      print(requried_value)
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
  wells_used <- label_matrix[fish_used_data == "x"]

  # read in data that matches valid genotypes
  genotype_data <- read_csv(genotyping_file) %>%
    select(genotyping_well = Well, genotype = Cluster) %>%
    filter(genotype %in% c("HET", "HOM", "WT")) %>% # do i want this?
    filter(genotyping_well %in% wells_used) %>%
    mutate(row = str_extract(genotyping_well, "[A-H]"), column = as.integer(str_extract(genotyping_well, "[0-9]+")))

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

microtracker_analysis <- function(data_file, genotypes) {
  data <- read_xlsx(data_file, sheet = "report", skip = 25, n_max = 96) %>%
    select(Well, `30`, `60`, `90`, `120`) %>%
    mutate(`Average Locomotor Activity` = rowMeans(across(c(`30`, `60`, `90`, `120`)))) %>%
    mutate(row = str_extract(Well, "[A-H]"), col = as.integer(str_extract(Well, "[0-9]+"))) %>%
    arrange(row, col) %>%
    mutate(row_id = row_number())

  genotypes <- process_genotypes() %>%
    select(genotype, row_id)

  finished_data <- data %>%
    full_join(genotypes, by = join_by(row_id)) %>%
    relocate(genotype) %>%
    select(genotype, `Average Locomotor Activity`) %>%
    arrange(genotype)

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
      alternations = sum(count[tetragram %in% c("LRLR", "RLRL")]) / sum(count),
      repetitions = sum(count[tetragram %in% c("LLLL", "RRRR")]) / sum(count),
      turns = sum(count)
    )

  # add genotyping
  processed_data <- processed_data %>%
    left_join(genotypes, by = join_by(ARENA == row_id)) %>%
    relocate(genotype) %>%
    arrange(genotype) %>%
    select(genotype, alternations, repetitions, turns)

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
  } else if (assay_name == "microtracker") {
    data <- microtracker_analysis(data_file, genotypes)
  } else if (assay_name == "y-maze 15") {
    data <- y_maze_analysis(data_file, genotypes)
  }

  if (idx == 1) {
    all_data <- data
  } else {
    all_data <- bind_rows(all_data, data)
  }
}

all_data <- all_data %>% arrange(genotype)
if (OUTPUT_FILE == "") {
  write_csv(all_data, "output.csv")
} else {
  write_csv(all_data, OUTPUT_FILE)
}
print(all_data)
