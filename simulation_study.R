#' ---
#' title: "Simulation Study"
#' author: "Md Abdul Basit"
#' output:
#'   pdf_document: default
#' ---
#' 
## ----setup, include=FALSE--------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE,
                      cache = TRUE, cache.path = "caches/simulation_study_cache/")

#' 
#' In this document, we reproduce the results presented in Table 1 (Section 4) of the submitted manuscript. Before executing the following code chunks, please note that--
#' 
#' 1. The simulation results include Monte Carlo errors, but we have verified that the fluctuation in results is negligible and does not affect any substantive conclusions drawn from them.
#' 
#' 2. Reproducing the simulation results is computationally expensive. The time required to render this document was about 9 hours on an Intel(R) Xeon(R) CPU E5-2650 server with 20 cores and 40 threads.
#' 
#' 3. As it takes a long time to run all the codes, some code chunks are not evaluated (`eval=FALSE`), and the results for those code chunks are loaded from `/intermediate_results/` instead. Please set `eval=TRUE` if it is needed to reproduce the results from those chunks.
#' 
#' 4. The function `simulation_final_results` doesn't return any result object. Instead, it pastes the results in `/intermediate_results/simulation_results_scenario-I` (for Scenario-I) or `/intermediate_results/simulation_results_scenario-II` (for Scenario-II). The main advantage of storing results in this approach is that even if the simulation is still running, we may observe the results for the completed simulation settings from the csv files stored in /intermediate_results/ instead of waiting for all the settings to finish.
#' 
#' 5. In order to reduce the computational time, the number of replications within each simulation setting and the number of bootstrap samples used to calculate the confidence intervals could be reduced to 500 by setting `R = 500` and `B = 500` in `simulation_final_results()`, respectively.
#' 
#' 6. If this Rmd file fails to render, delete the `./cache/` folder and then try rendering it again.
#' 
#' # Loading required packages and functions
#' 
#' First we load the necessary packages and functions to conduct the simulation study.
## ----loading-packages------------------------------------------------------------------------------
library(tidyverse)
library(parallel)
library(kableExtra)
source("code/simulation_functions.R")

#' 
#' # Obtaining true (population level) partially identified intervals
#' 
#' Next, we obtain the true partially identified intervals for each simulation setting using large scale numerical approximations. The computation of these true intervals should take up about 40 minutes to complete. The output is saved in  `/intermediate_results/population_level_partial_intervals.Rdata`.
#' 
## ----pop-level-intervals, eval=FALSE---------------------------------------------------------------
## ## The simulation settings
## contrast <- rep(rep(list(c(1,-1,0), c(1,0,-1), c(0,1,-1)), each = 6), 2)
## lambda <- rep(c(0, 0.1, 0.2, 0.5, 1, 2), 6)
## overlap <- rep(c(TRUE, FALSE), each = 18)
## 
## 
## partial_intervals <- round(
##   matrix(
##     unlist(
##       Map(f = function(x, y, z) {
##         pop_level_partial_intervals(seed = 1660, contrast = x, overlap = z,
##                                     lambda = y,size = 1e6)
##       }, x = contrast, y = lambda, z = overlap)
##     ), ncol = 2, byrow = T
##   ), 3
## )
## 
## colnames(partial_intervals) <- c("true_lower","true_upper")
## 
## save(partial_intervals,
##      file = "intermediate_results/population_level_partial_intervals.Rdata")

#' 
#' # Simulation results
#' 
#' Next, we set up a total of 36 simulation settings (18 under each scenario) and then reproduce the simulation results for these settings under each scenario separately. The intermediate results are stored in `/intermediate_results/`.
#' 
#' ## Simulation settings
#' 
## ----sim-settings----------------------------------------------------------------------------------
contrast <- rep(rep(list(c(1,-1,0), c(1,0,-1), c(0,1,-1)), each = 6), 1)
lambda <- rep(c(0, 0.1, 0.2, 0.5, 1, 2), 3)

## Loading the true (population-level) partially identified intervals
load("intermediate_results/population_level_partial_intervals.Rdata")

#' 
#' ## Simulation of Scenario-I (adequate overlap)
#' 
## ----scenario-I, eval=TRUE-------------------------------------------------------------------------
simulation_final_results(seed = 1660, overlap = TRUE, contrast = contrast,
                        lambda = lambda, size = 750, B = 1000, R = 1000,
                        parallel = TRUE, alpha = 0.1,
                        partial_intervals = partial_intervals)


#' 
#' 
#' 
#' ## Simulation of Scenario-II (lack of overlap)
#' 
## ----scenario-II, eval=FALSE-----------------------------------------------------------------------
## simulation_final_results(seed = 1660, overlap = FALSE, contrast = contrast,
##                         lambda = lambda, size = 750, B = 1000, R = 1000,
##                         parallel = TRUE, alpha = 0.1,
##                         partial_intervals = partial_intervals)
## 

#' 
#' 
#' ## Table 1 (Section 4)
#' 
## ---- eval=TRUE------------------------------------------------------------------------------------
scenario_I <- read_csv("intermediate_results/simulation_results_scenario-I.csv")
scenario_II <- read_csv("intermediate_results/simulation_results_scenario-II.csv")

combined_results <- rbind(scenario_I, scenario_II)

  
combined_results %>%
  mutate(
    Overlap = case_when(
      Overlap == 1 ~ "I",
      Overlap == 0 ~ "II"
    ),
    ATE = rep(rep(c("$\\tau_{1, 2}$", "$\\tau_{1, 3}$", "$\\tau_{2, 3}$"), each = 6), 2),
    bias_lower = paste0("$", sprintf("%.2f", bias_lower), "$"),
    bias_upper = paste0("$", sprintf("%.2f", bias_upper), "$"),
    partial_interval = paste0("$", "(", sprintf("%.3f", combined_results$true_lower), ", ",
                              sprintf("%.3f", combined_results$true_upper), ")", "$"),
    point_interval =  paste0("$", "(", sprintf("%.3f", combined_results$point_lower), ", ",
                              sprintf("%.3f", combined_results$point_upper), ")", "$"),
    conf_interval = paste0("$", "(", sprintf("%.3f", combined_results$conf_lower), ", ",
                              sprintf("%.3f", combined_results$conf_upper), ")", "$"),
  ) %>%
  select(Overlap, ATE, lambda, Lambda, bias_lower, bias_upper, non_coverage,
         partial_interval, point_interval, conf_interval) %>%
  kbl(booktabs = T, digits = 3, align = "ccrrrrrrrr", escape = F, format = "latex",
      col.names = linebreak(c("Scenario", "ATE", "$\\lambda$", "$\\Lambda$", "Lower", "Upper",
                    "Non-\n coverage", "Partially identified \n interval",
                    "Point estimate \n interval", "Confidence \n interval"), align = "c"),
      caption = "Simulation results presented in Table 1 (Section 4).") %>%
  add_header_above(c(" " = 4, "% Bias" = 2, " " = 4)) %>%
  kable_styling(latex_options = c("scale_down", "hold_position")) %>%
  collapse_rows(columns = 1:2, valign = "middle", latex_hline = "custom",
  custom_latex_hline = 1:2)


#' 
