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

# turn full genotype file into a table that can be merged with the data
process_genotypes <- function() {
  # create matrix to help with computations
  helper_matrix <- matrix(0:95, nrow = 8, ncol = 12, byrow = TRUE)
  print(helper_matrix)

  # get correct data and convert to 0-95 labelling
  genotype_data <- read_csv(GENOTYPING_FILE) %>%
    select(genotyping_well = Well, genotype = Cluster) %>%
    filter(genotype %in% c("HET", "HOM", "WT")) %>%
    mutate(well_number = as.integer((match(substr(genotyping_well, 1, 1), LETTERS[1:8]) - 1) * 12 + as.integer(substr(genotyping_well, 2, 3)) - 1))

  # correct to matrix to right format
  if (LABELLING_TYPE == 48) {
    if (LOCATION_48 == "left half") {
      helper_matrix <- t(helper_matrix[, 1:6])
    } else if (LOCATION_48 == "right half") {
      helper_matrix <- t(helper_matrix[, 7:12])
    } else if (LOCATION_48 == "top left") {
      helper_matrix <- helper_matrix[1:6, 1:8]
    } else {
      stop("LOCATION_48 must be \"left half\", \"right half\", or \"top left\"")
    }

    # filter data to the side we just specified
    genotype_data <- genotype_data %>%
      filter(well_number %in% helper_matrix)
  } else if (LABELLING_TYPE != 96) {
    stop("LABELLING_TYPE must be either 48 or 96")
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
    mutate(tetragram = paste0(turn_direction, next_turn, next_next_turn, next_next_next_turn)) %>% # probably a better way to do this...
    count(tetragram, name = "count") %>%
    mutate(percentage = count / sum(count) * 100) %>%
    ungroup()

  return(processed_data)
}

y_maze_analysis(data)
