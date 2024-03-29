source("code/sensitivity_analysis_functions.R")

library(tidyverse)

## Function for generaitng data for simulation
generate_data <- function(size, overlap) {
  X_1 <- rbinom(n = size, size = 1, prob = 0.5)
  X_2 <- runif(n = size, min = -1, max = 1)
  X_3 <- rnorm(n = size, mean = 0, sd = 0.5)
  X <- cbind(1, X_1, X_2, X_3)
  
  if (overlap) { ## Adequate overlap
    k2 <- 0.1
    k3 <- -0.1
  } else {       ## Lack of overlap
    k2 <- 3
    k3 <- 3
  }
  
  B_1 <- matrix(rep(0, 4))
  B_2 <- matrix(k2 * c(0, 1, 1, 1))
  B_3 <- matrix(k3 * c(0, 1, 1, -1))
  
  B <- cbind(B_1, B_2, B_3)
  p_A <- exp((X %*% B)) / rowSums((exp(X %*% B)))
  
  
  D <- matrix(rep(0, 3 * size), ncol = 3)
  for (i in 1:size) {
    D[i, ] <- rmultinom(n = 1, size = 1, prob = p_A[i, ])
  }
  
  g_1 <- matrix(c(1, 1, 1, 1))
  g_2 <- matrix(c(1, 1, -1, 1))
  g_3 <- matrix(c(1, 1, 1, -1))
  
  g <- cbind(g_1, g_2, g_3)
  
  p_Y <- exp((X %*% g)) / rowSums((exp(X %*% g)))
  
  Y <- matrix(rep(0, 3 * size), ncol = 3)
  for (i in 1:size) {
    Y[i, ] <- rmultinom(n = 1, size = 1, prob = p_Y[i, ])
  }
  
  
  dat_full <- data.frame(Y, D, X[, -1])
  names(dat_full) <- c("Y1", "Y2", "Y3", "D1", "D2", "D3", "X1", "X2", "X3")
  
  dat_obs <- dat_full %>%
    mutate(
      Y_obs = case_when(
        D1 == 1 ~ Y1,
        D2 == 1 ~ Y2,
        D3 == 1 ~ Y3
      ),
      A = case_when(
        D1 == 1 ~ 1,
        D2 == 1 ~ 2,
        D3 == 1 ~ 3
      ),
      A = factor(A)
    ) %>%
    select(Y_obs, A, starts_with("D"), starts_with("X"))
  
  return(dat_obs)
}



## Function for finding true (population level) partially identified intervals
pop_level_partial_intervals <- function(seed, contrast, overlap,
                                        lambda, size) {
  set.seed(seed)
  dat_obs <- generate_data(size, overlap)
  
  ## Estimating GPS
  fitted_probs <- nnet::multinom(
    formula = A ~ X1 + X2 + X3,
    data = dat_obs, trace = "FALSE"
  )$fitted.values
  
  
  partial_interval <-
    extrema.os.unified(
      A_name = "A", Y_name = "Y_obs",
      gps.formula = A ~ X1 + X2 + X3,
      contrast = contrast, data = dat_obs, lambda = lambda
    )
  
  return(partial_interval)
}


## Function for obtaining simulation results (single replicate) with
##  B = 1000 bootstrap samples
simulation_single_replicate <- function(data, contrast, lambda, B, alpha) {
  point_est <-
    extrema.os.unified(
      data = data, A_name = "A", Y_name = "Y_obs",
      gps.formula = A ~ X1 + X2 + X3, contrast = contrast,
      lambda = lambda
    )
  conf_int <-
    bootsens.os.unified(
      data = data, A_name = "A", Y_name = "Y_obs",
      gps.formula = A ~ X1 + X2 + X3, contrast = contrast,
      lambda = lambda, B = B, parallel = FALSE, alpha = alpha
    )
  return(c(point_est, conf_int))
}

# Function for obtaining simulation results (R replicate)
simulation_multiple_replicates <- function(datasets, contrast, lambda, B,
                                           parallel, true.lower, true.upper,
                                           alpha) {
  no.cores <- ifelse(parallel, detectCores(), 1)
  res <- mclapply(
    X = datasets, FUN = function(X) {
      simulation_single_replicate(X, contrast, lambda, B, alpha)
    },
    mc.cores = no.cores
  )
  res <- do.call(rbind, res)
  point_boot <- res[, c(1, 2)]
  conf_boot <- res[, c(3, 4)]
  
  R <- length(datasets)
  
  stdev_l <- sqrt((R - 1) / R * var(point_boot[, 1]))
  stdev_u <- sqrt((R - 1) / R * var(point_boot[, 2]))
  bias_l <- mean((point_boot[, 1] - true.lower) / stdev_l) * 100
  bias_u <- mean((point_boot[, 2] - true.upper) / stdev_u) * 100
  
  return(list(
    point_estimate = apply(point_boot, 2, median),
    stdev_l = round(stdev_l, 4),
    stdev_u = round(stdev_u, 4),
    bias_l = round(bias_l, 2),
    bias_u = round(bias_u, 2),
    confidence_interval = apply(conf_boot, 2, median),
    non_coverage = mean(ifelse(round(conf_boot[, 1], 3) <= true.lower &
                                 round(conf_boot[, 2], 3) >= true.upper, 0, 1))
  ))
}

## Function for obtaining the final simulation results under
## a given scenario (Scenario-I (overlap) or Scenario-II (lack of overlap))

simulation_final_results <- function(seed, overlap, contrast, lambda, size,
                                     B = 1000, R = 1000, parallel = TRUE,
                                     alpha = 0.1, partial_intervals) {
  set.seed(seed)
  datasets <- replicate(
    n = R,
    expr = generate_data(size = size, overlap = overlap),
    simplify = FALSE
  )
  
  row_names <-
    c(
      "Overlap", "C1", "C2", "C3", "lambda", "Lambda", "sd_lower", "sd_upper", "bias_lower",
      "bias_upper", "non_coverage", "true_lower", "true_upper", "point_lower",
      "point_upper", "conf_lower", "conf_upper"
    )
  
  file_name <- ifelse(overlap, "intermediate_results/simulation_results_scenario-I.csv",
                      "intermediate_results/simulation_results_scenario-II.csv"
  )
  
  if (overlap) {
    rows <- 1:18
  } else {
    rows <- 19:36
  }
  
  cat(row_names, "\n", file = file_name, sep = ",")
  
  out <-
    Map(f = function(x, y, rownum) {
      result <- simulation_multiple_replicates(
        datasets = datasets, contrast = x,
        lambda = y, B = B, parallel = parallel,
        alpha = alpha,
        true.lower = partial_intervals[rownum, 1],
        true.upper = partial_intervals[rownum, 2]
      )
      table_output <-
        c(
          overlap, x, y, round(exp(y), 2), result$stdev_l, result$stdev_u,
          result$bias_l, result$bias_u, result$non_coverage,
          partial_intervals[rownum, 1],
          partial_intervals[rownum, 2],
          round(result$point_estimate, 3),
          round(result$confidence_interval, 3)
        )
      cat(table_output, "\n", file = file_name, sep = ",", append = T)
    }, x = contrast, y = lambda, rownum = rows)
}
