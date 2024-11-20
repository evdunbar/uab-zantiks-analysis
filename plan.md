# How Will This Work?

## Goals

- usable once i'm gone
- possibly extendable once i'm gone
    - new assays
- understandable, not a black box
    - explanations shouldn't be forced
    - documentation should exist and be relatively easy to find
- relatively simple, cannot take forever to develop
    - could ask cammie if that is something she wants

## Ideal Workflow

1. download the experiment folder from box
2. open gui and point to that folder
    - possibly point to genotypes file as well?
3. select assay type
4. select wells for each assay?
5. get csv file that can be imported into prism for analysis

## What to Write It in?

### Python

- pros
    - analysis code already written
    - PySide6 seems like a good option
- cons
    - could be slow
    - lab doesn't know any python
    - hard to package for other platforms

### R

- pros
    - some of the lab knows R
    - fairly easy to translate analysis code
- cons
    - gui libraries seem lacking
    - packaging ease unknown

### Rust

- pros
    - definitely not slow
    - polars code is likely easy to translate to rust polars
    - easy to build for cross platform
- cons
    - even harder for lab to pick up than python if they want to add things
