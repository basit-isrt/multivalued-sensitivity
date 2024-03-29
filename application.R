#' ---
#' title: "Application"
#' author: "Md Abdul Basit"
#' output: 
#'   pdf_document: default
#' ---
#' 
## ----setup, include=FALSE--------------------------------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE, warning = FALSE, message = FALSE,
  cache = TRUE, cache.path = "caches/application_cache/"
)

#' 
#' In this document, we reproduce Table 2 and Figure 1 presented in Section 5 of the submitted manuscript. The time required to render this document was less than just two minutes on a Intel(R) Xeon(R) CPU E5-2650 server with 20 cores and 40 threads.
#' 
#' # Loading the dataset
#' 
## ----loading-packages------------------------------------------------------------------------------
library(tidyverse)
library(kableExtra)
library(ggplot2)

source("code/sensitivity_analysis_functions.R")
load(file = "data/nhanes.fish.rda")

#' 
#' # Creating the new treatment variable with three levels
#' 
## ----new-treatment---------------------------------------------------------------------------------
nhanes.fish.new <- nhanes.fish %>%
  mutate(fish.level.new = case_when(
    fish == 0 ~ "no",
    fish > 0 & fish <= 1 ~ "low",
    fish > 1 ~ "high"
  )) %>%
  mutate(A = case_when(
    fish == 0 ~ 1,
    fish > 0 & fish <= 1 ~ 2,
    fish > 1 ~ 3
  )) %>%
  mutate(Y = log2(o.LBXTHG))

#' 
#' # Sensitivity analysis of the pairwise average treatment effects
#' 
## ----sensitivity-analysis--------------------------------------------------------------------------
gps.formula <- A~gender + age + income + income.missing + race + education +
  smoking.ever + smoking.now


no_vs_low <- sensitivity_IPW(
  data = nhanes.fish.new, A_name = "A", Y_name = "Y",
  gps.formula = gps.formula, contrast = c(-1, 1, 0),
  lambda = c(0, 0.5, 0.75, 1, 1.5, 2), alpha = 0.1,
  parallel = T, B = 1000
)

no_vs_high <- sensitivity_IPW(
  data = nhanes.fish.new, A_name = "A", Y_name = "Y",
  gps.formula = gps.formula, contrast = c(-1, 0, 1),
  lambda = c(0, 0.5, 0.75, 1, 1.5, 2), alpha = 0.1,
  parallel = T, B = 1000
)

low_vs_high <- sensitivity_IPW(
  data = nhanes.fish.new, A_name = "A", Y_name = "Y",
  gps.formula = gps.formula, contrast = c(0, -1, 1),
  lambda = c(0, 0.5, 0.75, 1, 1.5, 2), alpha = 0.1,
  parallel = T, B = 1000
)

result <- as_tibble(rbind(no_vs_low, no_vs_high, low_vs_high))
names(result) <- c("lambda", "lower", "upper", "lower_5", "upper_95")

#' 
#' ## Table 2 (Section 5)
#' 
## ----results, eval=TRUE----------------------------------------------------------------------------
result_table <-
  result %>%
  mutate(
    Lambda = as.factor(round(exp(lambda), 2)),
    midpoint = (lower + upper) / 2,
    estimand = factor(
      rep(1:3, each = 6),
      levels = 1:3,
      labels = c(
        "$\\tau_{1, 2}$",
        "$\\tau_{1, 3}$",
        "$\\tau_{2, 3}$"
      )
    ),
    point_interval = paste0(
      "$", "(", sprintf("%.2f", result$lower), ", ",
      sprintf("%.2f", result$upper), ")", "$"
    ),
    conf_interval = paste0(
      "$", "(", sprintf("%.2f", result$lower_5), ", ",
      sprintf("%.2f", result$upper_95), ")", "$"
    ),
  )

result_table %>%
  select(estimand, Lambda, lambda, point_interval, conf_interval) %>%
  kbl(
    booktabs = T,
    align = "crrrr",
    digits = 2,
    escape = F,
    caption = "Sensitivity analysis results presented in Table 2 (Section 5)",
    col.names = linebreak(c(
      "Estimand", "$\\Lambda$", "$\\lambda$",
      "Point estimate \n interval",
      "90\\% confidence \n interval"
    ), align = "c")
  ) %>%
  kable_styling(position = "center", latex_options = "hold_position") %>%
  collapse_rows(
    columns = 1,
    valign = "middle",
    latex_hline = "major"
  )

#' 
#' 
#' ## Figure 1 (Section 5)
#' 
## ----sens-plot, eval=TRUE--------------------------------------------------------------------------
sens_plot <- ggplot(
  data = result_table,
  mapping = aes(x = Lambda, y = midpoint, colour = estimand, shape = estimand)
) +
  geom_errorbar(
    mapping = aes(ymin = lower_5, ymax = upper_95),
    width = 0.5,
    size = 0.5,
    position = position_dodge(width = 0.6),
    linetype = 1,
    show.legend = TRUE
  ) +
  geom_errorbar(
    mapping = aes(ymin = lower_5, ymax = upper_95),
    width = 0.5,
    position = position_dodge(width = 0.6),
    linetype = 2,
    show.legend = TRUE
  ) +
  geom_point(
    position = position_dodge(width = 0.6),
    size = 2.5,
    show.legend = TRUE
  ) +
  geom_hline(yintercept = 0, size = 0.5, linetype = 1, color = "grey") +
  xlab(expression(Sensitivity ~ parameter ~ Lambda)) +
  ylab("Average treatment effect (ATE)") +
  scale_color_brewer(
    palette = "Set1",
    labels = expression(tau[1 * "," * 2], tau[1 * "," * 3], tau[2 * "," * 3])
  ) +
  scale_shape(
    solid = TRUE,
    labels = expression(tau[1 * "," * 2], tau[1 * "," * 3], tau[2 * "," * 3])
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 10),
    legend.position = c(0.01, 0.99),
    legend.title = element_blank(),
    legend.justification = c("left", "top"),
    legend.text = element_text(size = 8),
    legend.background = element_rect(colour = "grey45")
  )

#' 
## ----Figure-2, echo=FALSE,fig.height=3.5, fig.width=6, fig.pos="!H", fig.cap= "\\label{fig:sens_plot}Graphical representation of the sensitivity analysis results  presented in Figure 1 (Section 5)."----
sens_plot

