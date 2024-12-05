#!/usr/bin/env Rscript

##############
# PARAMETERS #
##############

# Input the data file path
DATA_FILE <- "2024/11/01/ymaze_15/a/ymaze_15-20241101T164208.csv"

# Input the genotyping file path
GENOTYPING_FILE <- "2024/11/04/genotypes.csv"

# What do you want the output to be saved to?
# Name should end with .csv
# Leave empty for default of "output.csv"
OUTPUT_FILE <- ""

# Choose an assay type from:
#   - light/dark preference
#   - light/dark transition
#   - microtracker
#   - mirror biting
#   - social preference
#   - startle response/pre-pulse inhibition
#   - y-maze 15
#   - y-maze 4
ASSAY_NAME <- "y-maze 15"

# How do you want to specify fish? (sorry this is a lot of questions)
#   - Do you want to answer in terms of 48 or 96-well plate labelling?
#       48-well:            96-well:
#          1 2 3 4 5 6 7 8     1  2  3  4  5  6  7  8  9 10 11 12
#        A . . . . . . . .   A .  .  .  .  .  .  .  .  .  .  .  .
#        B . . . . . . . .   B .  .  .  .  .  .  .  .  .  .  .  .
#        C . . . . . . . .   C .  .  .  .  .  .  .  .  .  .  .  .
#        D . . . . . . . .   D .  .  .  .  .  .  .  .  .  .  .  .
#        E . . . . . . . .   E .  .  .  .  .  .  .  .  .  .  .  .
#        F . . . . . . . .   F .  .  .  .  .  .  .  .  .  .  .  .
#                            G .  .  .  .  .  .  .  .  .  .  .  .
#                            H .  .  .  .  .  .  .  .  .  .  .  .
LABELLING_TYPE <- 48
#   - If 48-well, YOU MUST specify if left half, right half, or top left within the 96-well genotyping plate
LOCATION_48 <- "left half"
#   - Count down columns or across rows first?
COLUMN_OR_ROW_FIRST <- "column"

# Which fish were used in this assay?
#   - Use the same order as the assay machine
#   - You can include ranges (e.g. "A1-B3") and individual wells (e.g. "F5")
FISH_USED <- c("A1-F2")


############
# ANALYSIS #
############

library(tidyverse)
library(dplyr)
library(readr)

# turn FISH_USED into a machine readable format
parse_fish_used <- function() {
  is_range <- str_detect(FISH_USED, "-")
  ranges <- str_split(FISH_USED[is_range], "-")
  individuals <- FISH_USED[!is_range]

  return(c(ranges, individuals))
}

# turn full genotype file into a table that can be merged with the data
process_genotypes <- function() {
  # get correct data and convert to 0-95 labelling
  genotype_data <- read_csv(GENOTYPING_FILE) %>%
    select(genotyping_well = Well, genotype = Cluster) %>%
    filter(genotype %in% c("HET", "HOM", "WT")) %>%
    mutate(well_id = as.integer((match(substr(genotyping_well, 1, 1), LETTERS[1:8]) - 1) * 12 +
                                  as.integer(substr(genotyping_well, 2, 3)) - 1))

  # create matrices to help with computations
  equivalence_matrix <- matrix(0:95, nrow = 8, ncol = 12, byrow = TRUE)
  if (LABELLING_TYPE == 96) {
    labeling_matrix <- matrix(outer(LETTERS[1:8], sprintf("%d", 1:12), FUN = paste0), nrow = 8, ncol = 12)
  } else if (LABELLING_TYPE == 48) {
    labeling_matrix <- matrix(outer(LETTERS[1:6], sprintf("%d", 1:8), FUN = paste0), nrow = 6, ncol = 8)
    if (LOCATION_48 == "left half") {
      equivalence_matrix <- t(equivalence_matrix[, 1:6])
    } else if (LOCATION_48 == "right half") {
      equivalence_matrix <- t(equivalence_matrix[, 7:12])
    } else if (LOCATION_48 == "top left") {
      equivalence_matrix <- equivalence_matrix[1:6, 1:8]
    } else {
      stop("LOCATION_48 must be \"left half\", \"right half\", or \"top left\"")
    }

    # filter data to the side we just specified
    genotype_data <- genotype_data %>%
      filter(well_id %in% equivalence_matrix)
  } else {
    stop("LABELLING_TYPE must be either 48 or 96")
  }

  # get matrix positions to access well id
  used_matrix <- matrix(FALSE, nrow = nrow(labeling_matrix), ncol = ncol(labeling_matrix))
  parsed_fish <- parse_fish_used()
  for (i in 1:(length(parsed_fish[1]) / 2)) {
    begin_well <- parsed_fish[1][i * 2]
    end_well <- parsed_fish[1][i * 2 + 1]
    begin_location <- which(labeling_matrix == begin, arr.ind = TRUE)
    end_location <- which(labeling_matrix == end, arr.ind = TRUE)
    begin_row <- begin_location[row]
    begin_col <- begin_location[col]
    end_row <- end_location[row]
    end_col <- end_location[col]
  }
  if (COLUMN_OR_ROW_FIRST == "row") {
  } else if (COLUMN_OR_ROW_FIRST == "column") {
  } else {
    stop("COLUMN_OR_ROW_FIRST must be either \"row\" or \"column\"")
  }

  print(genotype_data)
}

process_genotypes()

# read the file and get rid of the auto generated zantiks lines
# they would mess up csv parsing
lines <- readLines(DATA_FILE)
csv_text <- lines[4:(length(lines) - 1)]
data <- read_csv(I(csv_text), col_types = "dccici")

y_maze_analysis <- function(df) {
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
    mutate(tetragram = paste0(turn_direction,
                              next_turn,
                              next_next_turn,
                              next_next_next_turn)) %>% # probably a better way to do this...
    count(tetragram, name = "count") %>%
    mutate(percentage = count / sum(count) * 100) %>%
    ungroup()

  return(processed_data)
}

y_maze_analysis(data)
