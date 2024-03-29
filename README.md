# Supplementary Information for "Sensitivity Analysis of Causal Effects in Observational Studies with Multivalued Treatments"

### Author of the Code
- Md Abdul Basit (abasit@isrt.ac.bd)

## Software Versions
- R version: 4.3.2 (2023-10-31)
- Platform: x86_64-pc-linux-gnu (64-bit)
- Running under: Ubuntu 22.04.3 LTS

## Attached Packages
- kableExtra_1.3.4
- lubridate_1.9.3
- forcats_1.0.0
- stringr_1.5.0
- dplyr_1.1.3
- purrr_1.0.2
- readr_2.1.4
- tidyr_1.3.0
- tibble_3.2.1
- ggplot2_3.4.4
- tidyverse_2.0.0

## Loaded via a Namespace
- (Refer to the code for the full list)

---

## Content Description

This repository contains supplementary files for reproducing the analysis and figures presented in the manuscript. Please set the working directory of R to this folder before executing the files.

### `./code/`:
- `sensitivity_analysis_functions.R`: An R script containing functions for conducting sensitivity analysis.
- `simulation_functions.R`: An R script containing functions for conducting the simulation study.

### `./data/`:
- Subfolder containing `nhanes.fish.rda` extracted from the archived version of the R package CrossScreening.

### `./intermediate_results/`:
- Folder containing intermediate results.
- `results_sim.R`: An R script for generating results tables and figures.

### `./simulation_study.Rmd`:
- R Markdown file for performing simulations reported in the manuscript (Section 4).

### `./application.Rmd`:
- R Markdown file for reproducing tables and figures demonstrating the proposed framework (Section 5).

### `./simulation_study.pdf` and `./application.Rmd`:
- Rendered PDF outputs of respective R Markdown files.

### `./simulation_study.R` and `./application.R`:
- R scripts containing code chunks and comments from the corresponding R Markdown files.

