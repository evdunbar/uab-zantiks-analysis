#!/usr/bin/env Rscript

library(stringr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Provide exactly 1 input file argument")
}
if (!str_ends(args[1], fixed(".csv"))) {
  stop("Input file must be a csv file")
}

library(dplyr)
library(readr)

csv <- read_csv(args[1])
print(count(csv, genotype))
