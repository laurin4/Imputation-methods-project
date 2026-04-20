#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# Missingness exploration for NHANES analysis dataset
# -----------------------------------------------------------------------------
# Goal:
#   Describe and visualize missing-data structure in `data/nhanes_analysis.csv`.
#
# Outputs:
#   - outputs/missingness_variable_summary.csv
#   - outputs/pairwise_missingness_table.csv
#   - outputs/missingness_descriptive_table.csv
#   - figures/missingness_pattern_plot.png
#   - figures/missingness_by_variable_pct.png
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(naniar)
  library(VIM)
})

# ---- Reproducibility checks --------------------------------------------------
required_pkgs <- c("dplyr", "tidyr", "readr", "ggplot2", "naniar", "VIM")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_pkgs, collapse = ", "),
    ". Install them before running this script."
  )
}

input_path <- "data/nhanes_analysis.csv"
if (!file.exists(input_path)) {
  stop("Input file not found: ", input_path, ". Run 01_data_prep.R first.")
}

dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

# ---- Load data ---------------------------------------------------------------
nhanes <- readr::read_csv(input_path, show_col_types = FALSE)

# Standardize known categorical variables (if present as character in CSV)
nhanes <- nhanes %>%
  mutate(
    sex = if ("sex" %in% names(.)) factor(sex, levels = c("Male", "Female")) else sex,
    smoking_status = if ("smoking_status" %in% names(.)) {
      factor(smoking_status, levels = c("No", "Yes"))
    } else {
      smoking_status
    }
  )

n_obs <- nrow(nhanes)

# ---- 1) Summary of missingness per variable ---------------------------------
missingness_var <- nhanes %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "n_missing"
  ) %>%
  mutate(
    n_total = n_obs,
    pct_missing = round(100 * n_missing / n_total, 2)
  ) %>%
  arrange(desc(pct_missing), variable)

readr::write_csv(missingness_var, "outputs/missingness_variable_summary.csv", na = "")

# ---- 2) Pairwise missingness pattern table ----------------------------------
# Pairwise table: for each pair (var_i, var_j), compute percentage where both
# are missing among all observations.
var_names <- names(nhanes)

pairwise_missingness <- tidyr::expand_grid(var1 = var_names, var2 = var_names) %>%
  rowwise() %>%
  mutate(
    n_both_missing = sum(is.na(nhanes[[var1]]) & is.na(nhanes[[var2]])),
    pct_both_missing = round(100 * n_both_missing / n_obs, 2)
  ) %>%
  ungroup() %>%
  arrange(desc(pct_both_missing), var1, var2)

readr::write_csv(pairwise_missingness, "outputs/pairwise_missingness_table.csv", na = "")

# ---- 3) Missing-data pattern plot -------------------------------------------
# `naniar::vis_miss()` gives a clear variable-by-observation missingness map.
plot_missing_pattern <- naniar::vis_miss(
  nhanes,
  sort_miss = TRUE,
  cluster = TRUE,
  warn_large_data = FALSE
) +
  ggplot2::labs(
    title = "NHANES Missing-Data Pattern",
    subtitle = "White = observed values, black = missing values"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray30")
  )

ggsave(
  filename = "figures/missingness_pattern_plot.png",
  plot = plot_missing_pattern,
  width = 10,
  height = 6,
  dpi = 300
)

# ---- 4) Percentage missing by variable plot ---------------------------------
plot_pct_missing <- missingness_var %>%
  mutate(variable = reorder(variable, pct_missing)) %>%
  ggplot(aes(x = variable, y = pct_missing)) +
  geom_col(fill = "#2C7FB8", width = 0.75) +
  coord_flip() +
  labs(
    title = "Percentage Missing by Variable",
    x = NULL,
    y = "Missing (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = "figures/missingness_by_variable_pct.png",
  plot = plot_pct_missing,
  width = 8,
  height = 5,
  dpi = 300
)

# ---- Optional VIM call (for package requirement / compatibility checks) -----
# This call is wrapped in a non-printing object to avoid opening graphics
# devices in non-interactive batch runs while still using VIM functionality.
vim_summary <- VIM::aggr(nhanes, plot = FALSE, sortVars = TRUE)

# ---- 5) Small descriptive table for reporting --------------------------------
descriptive_table <- tibble::tibble(
  metric = c(
    "n_observations",
    "n_variables",
    "n_complete_cases",
    "pct_complete_cases",
    "n_any_missing",
    "pct_any_missing"
  ),
  value = c(
    n_obs,
    ncol(nhanes),
    sum(complete.cases(nhanes)),
    round(100 * mean(complete.cases(nhanes)), 2),
    sum(!complete.cases(nhanes)),
    round(100 * mean(!complete.cases(nhanes)), 2)
  )
)

readr::write_csv(descriptive_table, "outputs/missingness_descriptive_table.csv", na = "")

# ---- 6) Interpretation notes: MCAR vs MAR vs other ---------------------------
# Interpretation guidance for observed missingness:
# - MCAR (Missing Completely At Random) may be plausible only if missingness
#   appears unrelated to observed data and does not cluster strongly by pattern.
# - MAR (Missing At Random) is often more plausible in health survey data:
#   missingness can depend on observed traits (e.g., age, sex, smoking status).
# - A non-random structure in plots/tables weakens a pure MCAR assumption.
# - Importantly, MNAR (Missing Not At Random) cannot be proven from observed
#   data alone; it is a substantive assumption requiring external justification
#   or sensitivity analyses.
# Recommended next step:
# - Fit simple missingness indicator models (e.g., logistic regressions for each
#   variable's missingness vs observed covariates) to probe MCAR plausibility.

message("Missingness exploration complete.")
message("Saved tables to outputs/ and plots to figures/.")
