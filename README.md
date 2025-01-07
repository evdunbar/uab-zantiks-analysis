# Crowder Lab MicroTracker & Zantiks Analysis

## What Is This?

- This is a project to automatically analyze the output files from Zantiks and MicroTracker runs
- The goal is to combine genotype data and run data to create an output that is easy to copy and paste into GraphPad Prism

## Okay, How Do I Use It?

1. To run an analysis, you will need:
    - The output file from a MicroTracker or Zantiks run
    - Genotyping data for the plate you ran as a .csv file (this can be exported from an HRM)
    - A .txt file that indicates which fish from your genotyping data were used in your run
    - The answers to a few questions
2. This information needs to be added to the file [analyze.R](analyze.R) in the section "PARAMETERS"
    - If you think you will run this analysis multiple times, consider using a config file
        - The information required is the same, but this allows you to save the information for later use
3. To run multiple analyses at the same time, list each of the items in 1. in the same order
4. Run [analyze.R](analyze.R) and open the output file
5. Copy and paste your results into GraphPad Prism!

## What Are _________ ?

### The Main Files and Folders in This Repository

- [README.md](README.md)
    - This file
- [analyze.R](analyze.R)
    - The brains of the project
    - Fill out the info at the top and run it to get your output!
- [configs](configs/)
    - This folder contains files that can be fed into [analyze.R](analyze.R)
    - Check this out if you think you might need to run the same analysis multiple times
- [templates](templates/)
    - This folder contains templates to reuse for hand-made genotyping files and fish used files
    - Try making copies of these to preserve the original
- [extras](extras/)
    - This folder contains R scripts that may or may not be useful to you
    - They do extra things like count the number of each genotype in an output file

### Config Files

- These files let you save the PARAMETERS used to run analysis for later
- An example to start from is located in the [templates folder](templates/config.toml)
- They are written in a file format called "toml"
    - `[` and `]` on their own line denote the name of a section
    - `"` and `"` surround a *string*
    - `[` and `]` on the right hand side of an `=` represent a list
    - The left hand side of an `=` is the name of a field/variable
    - `#` denotes the start of a comment; Everything to the right of this will not be read by the computer

### Fish Used Files

- These files denote both which fish were used in a specific assay
- Common examples to start from are located in the [templates/fish_used/](templates/fish_used/)
- "O" (capital letter o) represents that this well did not have a fish in it and "X" or "x" represents that it did
- For example, a full 96-well plate for a MicroTracker Run would look like this:
```
x x x x x x x x x x x x
x x x x x x x x x x x x
x x x x x x x x x x x x
x x x x x x x x x x x x
x x x x x x x x x x x x
x x x x x x x x x x x x
x x x x x x x x x x x x
x x x x x x x x x x x x
```
