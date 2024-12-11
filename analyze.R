#!/usr/bin/env Rscript

##############
# PARAMETERS #
##############

# Input the data file path
DATA_FILE <- "2024/06/11/mt_plate_2_25um_amantadine_clutch_1_6_dpf_rb.xlsx"

# Input the genotyping file path
GENOTYPING_FILE <- "2024/06/11/genotypes.csv"

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
ASSAY_NAME <- "microtracker"

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
LABELLING_TYPE <- 96
#   - If 48-well, YOU MUST specify if left half, right half, or top left within the 96-well genotyping plate
LOCATION_48 <- "left half"
#   - Count down columns or across rows first?
COLUMN_OR_ROW_FIRST <- "row"

# Which fish were used in this assay?
#   - Use the same order as the assay machine
#   - You can include ranges (e.g. "A1-B3") and individual wells (e.g. "F5")
FISH_USED <- c("A1-H12")


############
# ANALYSIS #
############

library(tidyverse)
library(dplyr)
library(readr)
library(readxl)

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
    num_rows <- 8
    num_cols <- 12
  } else if (LABELLING_TYPE == 48) {
    num_rows <- 6
    num_cols <- 8
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
  labeling_matrix <- matrix(outer(LETTERS[1:num_rows], sprintf("%d", 1:num_cols), FUN = paste0),
                            nrow = num_rows,
                            ncol = num_cols)

  # get matrix positions to access well id
  used_matrix <- matrix(FALSE, nrow = nrow(labeling_matrix), ncol = ncol(labeling_matrix))
  parsed_fish <- parse_fish_used()
  ranges <- parsed_fish[[1]]
  # iterate over each range by pair
  for (i in 1:(length(ranges) / 2)) {
    begin_well <- ranges[i * 2 - 1]
    end_well <- ranges[i * 2]
    begin_location <- which(labeling_matrix == begin_well, arr.ind = TRUE)
    end_location <- which(labeling_matrix == end_well, arr.ind = TRUE)
    begin_row <- begin_location[1, 1]
    begin_col <- begin_location[1, 2]
    end_row <- end_location[1, 1]
    end_col <- end_location[1, 2]

    if (COLUMN_OR_ROW_FIRST == "row") {
      # are we filling whole rows?
      # note that LABELLING_TYPE has already been validated
      rows_full <- begin_col == 1 & end_col == num_cols
      if (rows_full) {
        used_matrix[begin_row:end_row, ] <- TRUE
      } else {
        # are we one row?
        if (begin_row == end_row) {
          used_matrix[begin_row, begin_col:end_col]
        } else {
          # first row
          used_matrix[begin_row, begin_col:num_cols] <- TRUE
          # last row
          used_matrix[end_row, 1:end_col] <- TRUE
          # middle rows
          if (end_row - begin_row > 1) {
            used_matrix[(begin_row + 1):(end_row - 1), ] <- TRUE
          }
        }
      }
    } else if (COLUMN_OR_ROW_FIRST == "column") {
      # are we filling whole columns?
      # note that LABELLING_TYPE has already been validated
      cols_full <- begin_row == 1 & end_row == num_rows
      if (cols_full) {
        used_matrix[, begin_col:end_col] <- TRUE
      } else {
        # are we one column?
        if (begin_col == end_col) {
          used_matrix[begin_row:end_row, begin_col]
        } else {
          # first column
          used_matrix[begin_row:num_rows, begin_col] <- TRUE
          # last column
          used_matrix[1:end_row, end_col] <- TRUE
          # middle columns
          if (end_col - begin_col > 1) {
            used_matrix[, (begin_col + 1):(end_col - 1)] <- TRUE
          }
        }
      }
    } else {
      stop("COLUMN_OR_ROW_FIRST must be either \"row\" or \"column\"")
    }
  }
  if (length(parsed_fish) > 1) {
    singletons <- parsed_fish[[2]]
    for (well in singletons) {
      used_matrix[which(labeling_matrix == well)] <- TRUE
    }
  }

  # filter to only the used fish!!!!
  genotype_data <- genotype_data %>%
    filter(well_id %in% equivalence_matrix[used_matrix]) %>%
    arrange(genotyping_well)

  # add row IDs
  genotype_data <- genotype_data %>%
    mutate(row_id = row_number())

  return(genotype_data)
}

microtracker_analysis <- function() {
  data <- read_xlsx(DATA_FILE, sheet = "report", skip = 25, n_max = 96) %>%
    select(Well, `30`, `60`, `90`, `120`) %>%
    mutate(average = rowMeans(across(c(`30`, `60`, `90`, `120`))))

  genotypes <- process_genotypes()
  print(genotypes)

  return(data)
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
    mutate(tetragram = paste0(turn_direction,
                              next_turn,
                              next_next_turn,
                              next_next_next_turn)) %>% # probably a better way to do this...
    count(tetragram, name = "count") %>%
    mutate(percentage = count / sum(count) * 100) %>%
    ungroup()

  # add genotyping
  genotyping_data <- process_genotypes() %>%
    select(row_id, genotype)
  processed_data <- processed_data %>%
    left_join(genotyping_data, by = join_by(ARENA == row_id)) %>%
    relocate(genotype) %>%
    arrange(genotype)

  return(processed_data)
}

microtracker_analysis()
