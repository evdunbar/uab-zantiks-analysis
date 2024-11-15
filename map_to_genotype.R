#!/usr/bin/env Rscript

# Goals:
# - take in a full(ish) genotyping file generated from hrm that's been run
#   through the clean_genotypes.R script
# - take in numbers 1-8 that represent which 12 fish we're doing
# - arranged like:
#   1 5
#   2 6
#   3 7
#   4 8
# - assuming we're split into two 48 well plates, then 1 would be equivalent
#   to plate 1 columns 1 & 2

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
num_args <- length(args)
if (num_args < 2) {
  stop("Supply at least two arguments", call. = FALSE)
}
if (num_args > 8) {
  stop("Cannot specify more than 7 plate zones", call. = FALSE)
}

file_path <- args[1]
file_dir <- dirname(file_path)
zones <- args[2:num_args]

genotypes <- read_csv(file_path)
