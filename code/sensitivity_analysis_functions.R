## Function to obtain the partially identified point estimate intervals
## of counterfactual means under the proposed framework
cf.mean.estimate <- function(A, Y, lambda = 0, fitted.probs) {
  fitted.logit <- qlogis(fitted.probs)

  eg <- exp(-fitted.logit)
  Y <- Y[A == 1]
  eg <- eg[A == 1]
  eg <- eg[order(-Y)]
  Y <- Y[order(-Y)]

  ## maximization
  num.each.low <- Y * (1 + exp(-lambda) * eg)
  num.each.up <- Y * (1 + exp(lambda) * eg)
  num <-
    c(0, cumsum(num.each.up)) + c(rev(cumsum(rev(num.each.low))), 0)

  den.each.low <- (1 + exp(-lambda) * eg)
  den.each.up <- (1 + exp(lambda) * eg)
  den <-
    c(0, cumsum(den.each.up)) + c(rev(cumsum(rev(den.each.low))), 0)

  maximum <- max(num / den)
  ## print(den[which.max(num/den)] / n)

  ## minimization
  num <-
    c(0, cumsum(num.each.low)) + c(rev(cumsum(rev(num.each.up))), 0)
  den <-
    c(0, cumsum(den.each.low)) + c(rev(cumsum(rev(den.each.up))), 0)
  minimum <- min(num / den)
  ## print(den[which.min(num/den)] / n)

  point_estimate <- c(minimum, maximum)
  names(point_estimate) <- c("lower", "upper")
  return(point_estimate)
}


## Function to estimate the generalized propensity scores
gps.estimate <- function(data, A_name, gps.formula, ordinal = FALSE) {
  if (typeof(gps.formula) != "character") {
    gps.formula <- Reduce(paste0, deparse(gps.formula))
  }
  gps.formula <- gsub(" ", "", gps.formula)

  if (ordinal) {
    fitted.probs <- VGAM::vglm(
      formula = gps.formula,
      family = cratio(parallel = F, reverse = F),
      data = data
    )@fitted.values
  } else {
    fitted.probs <- nnet::multinom(
      formula = gps.formula,
      data = data,
      trace = F
    )$fitted.values
  }

  ## Creating the propensity score matrix in the case of binary treatment
  if (dim(fitted.probs)[2] == 1) {
    fitted.probs <- cbind(1 - fitted.probs, fitted.probs)
  }


  return(fitted.probs)
}


## Function to obtain the partially identified interval of the target estimand
extrema.os.unified <-
  function(data,
           A_name,
           Y_name,
           gps.formula,
           contrast,
           lambda = 0,
           ordinal = FALSE) {
    fitted.probs <- gps.estimate(data, A_name, gps.formula, ordinal)

    treatment <- c(0, 0)
    control <- c(0, 0)

    dummies <- fastDummies::dummy_cols(data[, A_name])[, -1]
    colnames(dummies) <- paste0("A", seq_along(contrast))

    data_dummies <- cbind(data, dummies)

    for (i in seq_along(contrast)) {
      if (contrast[i] > 0) {
        cf.mean <-
          cf.mean.estimate(
            data_dummies[, paste0("A", i)],
            data[, Y_name],
            lambda,
            fitted.probs[, i]
          )

        treatment <- treatment + contrast[i] * cf.mean
      } else if (contrast[i] < 0) {
        cf.mean <-
          cf.mean.estimate(
            data_dummies[, paste0("A", i)],
            data[, Y_name],
            lambda,
            fitted.probs[, i]
          )
        control <- control + contrast[i] * cf.mean
      }
    }

    return(treatment + rev(control))
  }

## Function to obtain a (1 - alpha)% confidence interval for the target estimand
bootsens.os.unified <-
  function(data,
           A_name,
           Y_name,
           gps.formula,
           contrast,
           lambda = 0,
           alpha = 0.1,
           parallel = FALSE,
           B = 1000,
           ordinal = FALSE) {
    no.cores <- parallel::detectCores()
    n <- dim(data)[1]

    if (parallel) {
      out <- parallel::mclapply(1:B, function(iter) {
        s <- sample(1:n, n, TRUE)
        res <- tryCatch(
          extrema.os.unified(
            data = data[s, ], A_name, Y_name,
            gps.formula, contrast, lambda, ordinal
          ),
          error = function(e) {
            print(e)
          }
        )
        res
      }, mc.cores = no.cores)
    } else {
      out <- lapply(1:B, function(iter) {
        s <- sample(1:n, n, TRUE)
        res <- tryCatch(
          extrema.os.unified(
            data = data[s, ], A_name, Y_name,
            gps.formula, contrast, lambda, ordinal
          ),
          error = function(e) {
            print(e)
          }
        )
        res
      })
    }

    out <- do.call(rbind, out)

    c(
      quantile(out[, 1], alpha / 2, na.rm = TRUE),
      quantile(out[, 2], 1 - alpha / 2, na.rm = TRUE)
    )
  }

## Function to conduct sensitivity analysis for an array of sensitivity
## parameters lambda
sensitivity_IPW <- function(data, A_name, Y_name, gps.formula, contrast, lambda = 0,
                            alpha = 0.1, parallel = F, B = 1000, ordinal = FALSE) {
  if (length(lambda) == 1) {
    res <- c(
      lambda,
      extrema.os.unified(
        data, A_name, Y_name, gps.formula, contrast,
        lambda, ordinal
      ),
      bootsens.os.unified(
        data, A_name, Y_name, gps.formula, contrast,
        lambda, alpha, parallel, B, ordinal
      )
    )
  } else if (length(lambda) > 1) {
    res <- lapply(
      X = as.list(lambda),
      function(X, ...) {
        c(
          X,
          extrema.os.unified(
            data, A_name, Y_name, gps.formula, contrast,
            X, ordinal
          ),
          bootsens.os.unified(
            data, A_name, Y_name, gps.formula, contrast,
            X, alpha, parallel, B, ordinal
          )
        )
      }
    )
    res <- do.call(rbind, res)
  }
  return(res)
}
