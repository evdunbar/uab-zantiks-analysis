#!/usr/bin/env Rscript

########
# GOAL #
########
# This script is to process microtracker/zantiks data along with genotyping data (from HRM)
# and output a file to copy and paste into Prism


##############
# PARAMETERS #
##############

# Input the data file paths
# DATA_FILE_PREFIX will be prepended to every string in DATA_FILES
DATA_FILE_PREFIX <- "2024/"
DATA_FILES <- c(
  "11/01/ymaze_15/a/ymaze_15-20241101T164208.csv",
  "11/01/ymaze_15/b/ymaze_15-20241101T175733.csv",
  "11/01/ymaze_15/c/ymaze_15-20241101T191455.csv",
  "11/01/ymaze_15/d/ymaze_15-20241101T202944.csv",
  "11/04/ymaze_15/a/ymaze_15-20241104T170903.csv",
  "11/04/ymaze_15/b/ymaze_15-20241104T182427.csv",
  "11/04/ymaze_15/c/ymaze_15-20241104T194506.csv",
  "12/09/ymaze_15/a/ymaze_15-20241209T121317.csv",
  "12/09/ymaze_15/b/ymaze_15-20241209T130440.csv",
  "12/09/ymaze_15/c/ymaze_15-20241209T142016.csv"
)

# Input the genotyping file paths
# Each genotyping file corresponds to one data file
# GENOTYPING_FILE_PREFIX will be prepended to every string in GENOTYPING_FILES
GENOTYPING_FILE_PREFIX <- "2024/"
GENOTYPING_FILES <- c(
  "11/01/genotypes.csv",
  "11/01/genotypes.csv",
  "11/01/genotypes.csv",
  "11/01/genotypes.csv",
  "11/04/genotypes.csv",
  "11/04/genotypes.csv",
  "11/04/genotypes.csv",
  "12/09/genotypes.csv",
  "12/09/genotypes.csv",
  "12/09/genotypes.csv"
)

# Where should the output to be saved to?
# Name should end with .csv
# Leave empty for default of "output.csv"
OUTPUT_FILE <- ""

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
  "y-maze 15",
  "y-maze 15",
  "y-maze 15",
  "y-maze 15",
  "y-maze 15",
  "y-maze 15",
  "y-maze 15",
  "y-maze 15",
  "y-maze 15",
  "y-maze 15"
)

# Input the fish used paths
# COUNTING_DIRECTIONS can be either across or down
FISH_USED_PREFIX <- "2024/"
FISH_USED <- c(
  "11/01/ymaze_15/a/fish_used.txt",
  "11/01/ymaze_15/b/fish_used.txt",
  "11/01/ymaze_15/c/fish_used.txt",
  "11/01/ymaze_15/d/fish_used.txt",
  "11/04/ymaze_15/a/fish_used.txt",
  "11/04/ymaze_15/b/fish_used.txt",
  "11/04/ymaze_15/c/fish_used.txt",
  "12/09/ymaze_15/a/fish_used.txt",
  "12/09/ymaze_15/b/fish_used.txt",
  "12/09/ymaze_15/c/fish_used.txt"
)
COUNTING_DIRECTIONS <- c(
  "across",
  "across",
  "across",
  "across",
  "across",
  "across",
  "across",
  "across",
  "across",
  "across"
)


############
# ANALYSIS #
############

library(tidyverse)
library(dplyr)
library(readr)
library(readxl)

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
  if (counting_direction == "across") {
    wells_used <- sort(wells_used)
  }

  # read in data that matches valid genotypes
  genotype_data <- read_csv(genotyping_file) %>%
    select(genotyping_well = Well, genotype = Cluster) %>%
    filter(genotype %in% c("HET", "HOM", "WT")) %>%
    filter(genotyping_well %in% wells_used)

  print(genotype_data)
}

microtracker_analysis <- function() {
  data <- read_xlsx(DATA_FILE, sheet = "report", skip = 25, n_max = 96) %>%
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

y_maze_analysis <- function() {
  # read the file and get rid of the auto generated zantiks lines
  # they would mess up csv parsing
  lines <- readLines(DATA_FILE)
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
  genotyping_data <- process_genotypes() %>%
    select(row_id, genotype)
  processed_data <- processed_data %>%
    left_join(genotyping_data, by = join_by(ARENA == row_id)) %>%
    relocate(genotype) %>%
    arrange(genotype) %>%
    select(genotype, alternations, repetitions, turns)

  return(processed_data)
}

# main loop!
for (idx in seq_along(DATA_FILES)) {
  # get correct files by corresponding index
  data_file <- paste0(DATA_FILE_PREFIX, DATA_FILES[idx])
  genotyping_file <- paste0(GENOTYPING_FILE_PREFIX, GENOTYPING_FILES[idx])
  fish_used_file <- paste0(FISH_USED_PREFIX, FISH_USED[idx])
  assay_name <- ASSAY_NAMES[idx]
  counting_direction <- COUNTING_DIRECTIONS[idx]

  process_genotypes(genotyping_file, fish_used_file, counting_direction)

  # write_csv(output, "output.csv")
}
