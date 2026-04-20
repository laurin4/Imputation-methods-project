#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# Generate artificial missingness datasets for benchmarking
# -----------------------------------------------------------------------------
# Goal:
#   Create NHANES-based benchmarking datasets with artificial missingness under:
#   - MCAR (completely at random)
#   - MAR  (depends on observed covariates via logistic probabilities)
#
# Inputs:
#   - data/nhanes_analysis.csv
#
# Outputs (in sim_data/):
#   - datasets/sim_<mechanism>_r<rate>_rep<rep>.csv
#   - masks/mask_<mechanism>_r<rate>_rep<rep>.csv
#   - simulation_manifest.csv
#
# Notes:
#   - Missingness is induced in predictor variables only:
#       age, sex, bmi, systolic_bp, smoking_status, income_ratio
#   - hba1c is preserved as a complete target variable for benchmarking.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
})

# ---- Reproducibility checks --------------------------------------------------
required_pkgs <- c("dplyr", "tidyr", "readr", "purrr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_pkgs, collapse = ", "),
    ". Install them before running this script."
  )
}

# ---- Configuration -----------------------------------------------------------
input_path <- "data/nhanes_analysis.csv"

analysis_vars <- c(
  "hba1c", "age", "sex", "bmi", "systolic_bp", "smoking_status", "income_ratio"
)

# Predictors to be artificially masked
mask_vars <- c("age", "sex", "bmi", "systolic_bp", "smoking_status", "income_ratio")

missing_rates <- c(0.10, 0.20, 0.30)
mechanisms <- c("MCAR", "MAR")
n_reps <- 100

output_root <- "sim_data"
output_data_dir <- file.path(output_root, "datasets")
output_mask_dir <- file.path(output_root, "masks")

dir.create(output_data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_mask_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(input_path)) {
  stop("Input file not found: ", input_path, ". Run 01_data_prep.R first.")
}

# ---- Runtime controls for balanced CPU usage --------------------------------
# Optional environment variables:
#   FAST_MODE        -> "1" enables quicker defaults
#   N_REPS           -> override number of repetitions per mechanism x level
#   REF_MAX_ROWS     -> cap complete-case reference rows used in simulations
fast_mode <- identical(Sys.getenv("FAST_MODE", unset = "0"), "1")
n_reps_env <- as.integer(Sys.getenv("N_REPS", unset = if (fast_mode) "20" else as.character(n_reps)))
ref_max_rows <- as.integer(Sys.getenv("REF_MAX_ROWS", unset = if (fast_mode) "2000" else "0"))

if (is.na(n_reps_env) || n_reps_env <= 0) {
  stop("N_REPS must be a positive integer.")
}
if (is.na(ref_max_rows) || ref_max_rows < 0) {
  stop("REF_MAX_ROWS must be a non-negative integer.")
}
n_reps <- n_reps_env

# ---- Utility helpers ---------------------------------------------------------

# Convert rates to labels suitable for filenames and metadata.
rate_to_label <- function(rate) {
  paste0(sprintf("%02d", as.integer(round(rate * 100))), "pct")
}

# Safe inverse-logit
inv_logit <- function(x) {
  1 / (1 + exp(-x))
}

# Build complete-case reference dataset on the final analysis variables.
build_reference_data <- function(df, vars) {
  missing_cols <- setdiff(vars, names(df))
  if (length(missing_cols) > 0) {
    stop("Required variable(s) missing from input: ", paste(missing_cols, collapse = ", "))
  }

  df %>%
    select(all_of(vars)) %>%
    tidyr::drop_na() %>%
    mutate(row_id = row_number(), .before = 1)
}

# MCAR masking:
# Randomly mask exactly floor(rate * N * P) cells among chosen columns.
apply_mcar_mask <- function(df, vars_to_mask, missing_rate, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n <- nrow(df)
  p <- length(vars_to_mask)
  n_cells <- n * p
  n_mask <- as.integer(round(missing_rate * n_cells))

  if (n_mask <= 0) {
    return(list(data = df, mask_index = tibble::tibble()))
  }

  cell_grid <- tidyr::expand_grid(
    row_id = seq_len(n),
    variable = vars_to_mask
  ) %>%
    mutate(cell_id = row_number())

  sampled_ids <- sample(cell_grid$cell_id, size = n_mask, replace = FALSE)
  mask_index <- cell_grid %>%
    filter(.data$cell_id %in% sampled_ids) %>%
    select(.data$row_id, .data$variable)

  # Capture true values before masking for evaluation.
  mask_index <- mask_index %>%
    mutate(
      true_value = purrr::map2_chr(
        .data$row_id, .data$variable,
        ~as.character(df[[.y]][.x])
      )
    )

  out <- df
  for (j in seq_len(nrow(mask_index))) {
    out[[mask_index$variable[j]]][mask_index$row_id[j]] <- NA
  }

  list(data = out, mask_index = mask_index)
}

# MAR masking:
# Missingness probability per cell depends on observed age, sex, bmi.
# We calibrate an intercept so mean probability approximately equals target rate.
apply_mar_mask <- function(df, vars_to_mask, missing_rate, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n <- nrow(df)
  p <- length(vars_to_mask)
  n_cells <- n * p
  n_target <- as.integer(round(missing_rate * n_cells))

  if (n_target <= 0) {
    return(list(data = df, mask_index = tibble::tibble()))
  }

  # Standardized predictors for stable logistic probabilities.
  z_age <- as.numeric(scale(df$age))
  z_bmi <- as.numeric(scale(df$bmi))
  sex_male <- ifelse(df$sex == "Male", 1, 0)

  # Variable-specific coefficients (moderate effect sizes).
  var_coefs <- tibble::tibble(
    variable = vars_to_mask,
    b_age = c(0.45, 0.35, 0.40, 0.50, 0.30, 0.25),
    b_male = c(0.25, 0.20, 0.20, 0.15, 0.40, 0.10),
    b_bmi = c(0.35, 0.20, 0.45, 0.30, 0.20, 0.40)
  )

  # Build row-variable cell grid with MAR linear predictor.
  mar_grid <- tidyr::expand_grid(
    row_id = seq_len(n),
    variable = vars_to_mask
  ) %>%
    left_join(var_coefs, by = "variable") %>%
    mutate(
      lp_no_intercept = .data$b_age * z_age[.data$row_id] +
        .data$b_male * sex_male[.data$row_id] +
        .data$b_bmi * z_bmi[.data$row_id]
    )

  # Intercept calibration to match the desired average masking probability.
  calibrate_intercept <- function(alpha) {
    mean(inv_logit(alpha + mar_grid$lp_no_intercept)) - missing_rate
  }
  alpha <- uniroot(calibrate_intercept, lower = -10, upper = 10)$root

  mar_grid <- mar_grid %>%
    mutate(prob_miss = inv_logit(alpha + .data$lp_no_intercept))

  # Sample exactly n_target cells weighted by MAR probabilities.
  sampled_rows <- sample(
    x = seq_len(nrow(mar_grid)),
    size = n_target,
    replace = FALSE,
    prob = mar_grid$prob_miss
  )

  mask_index <- mar_grid %>%
    slice(sampled_rows) %>%
    select(.data$row_id, .data$variable)

  mask_index <- mask_index %>%
    mutate(
      true_value = purrr::map2_chr(
        .data$row_id, .data$variable,
        ~as.character(df[[.y]][.x])
      )
    )

  out <- df
  for (j in seq_len(nrow(mask_index))) {
    out[[mask_index$variable[j]]][mask_index$row_id[j]] <- NA
  }

  list(data = out, mask_index = mask_index)
}

# Dispatch wrapper for missingness mechanism.
generate_one_simulation <- function(ref_df, mechanism, missing_rate, repetition, seed_base = 20260420) {
  mechanism <- toupper(mechanism)
  if (!mechanism %in% c("MCAR", "MAR")) {
    stop("Unknown mechanism: ", mechanism, ". Use 'MCAR' or 'MAR'.")
  }

  # Deterministic per-run seed for reproducibility.
  # Keeps each mechanism/rate/rep combination stable across reruns.
  mech_id <- ifelse(mechanism == "MCAR", 1L, 2L)
  run_seed <- seed_base + mech_id * 100000L + as.integer(round(missing_rate * 1000)) * 100L + repetition

  masked <- if (mechanism == "MCAR") {
    apply_mcar_mask(ref_df, mask_vars, missing_rate, seed = run_seed)
  } else {
    apply_mar_mask(ref_df, mask_vars, missing_rate, seed = run_seed)
  }

  sim_data <- masked$data %>%
    mutate(
      repetition = repetition,
      mechanism = mechanism,
      missing_rate = missing_rate,
      .before = 1
    )

  mask_table <- masked$mask_index %>%
    mutate(
      repetition = repetition,
      mechanism = mechanism,
      missing_rate = missing_rate,
      .before = 1
    ) %>%
    arrange(.data$row_id, .data$variable)

  list(data = sim_data, mask = mask_table)
}

# ---- Load and prepare reference subset --------------------------------------
raw_data <- readr::read_csv(input_path, show_col_types = FALSE)

reference_data <- build_reference_data(raw_data, analysis_vars)

if (nrow(reference_data) == 0) {
  stop("Reference dataset has 0 rows after complete-case filtering.")
}

if (ref_max_rows > 0 && nrow(reference_data) > ref_max_rows) {
  set.seed(20260420)
  reference_data <- reference_data %>%
    slice_sample(n = ref_max_rows)
}

# ---- Run simulation grid -----------------------------------------------------
sim_grid <- tidyr::expand_grid(
  mechanism = mechanisms,
  missing_rate = missing_rates,
  repetition = seq_len(n_reps)
) %>%
  mutate(rate_label = purrr::map_chr(missing_rate, rate_to_label))

manifest <- vector("list", length = nrow(sim_grid))

for (i in seq_len(nrow(sim_grid))) {
  mechanism_i <- sim_grid$mechanism[i]
  rate_i <- sim_grid$missing_rate[i]
  rep_i <- sim_grid$repetition[i]
  rate_label_i <- sim_grid$rate_label[i]

  sim_result <- generate_one_simulation(
    ref_df = reference_data,
    mechanism = mechanism_i,
    missing_rate = rate_i,
    repetition = rep_i
  )

  data_file <- file.path(
    output_data_dir,
    sprintf("sim_%s_r%s_rep%03d.csv", tolower(mechanism_i), rate_label_i, rep_i)
  )
  mask_file <- file.path(
    output_mask_dir,
    sprintf("mask_%s_r%s_rep%03d.csv", tolower(mechanism_i), rate_label_i, rep_i)
  )

  readr::write_csv(sim_result$data, data_file, na = "")
  readr::write_csv(sim_result$mask, mask_file, na = "")

  manifest[[i]] <- tibble::tibble(
    mechanism = mechanism_i,
    missing_rate = rate_i,
    repetition = rep_i,
    n_rows = nrow(sim_result$data),
    n_masked_cells = nrow(sim_result$mask),
    dataset_file = data_file,
    mask_file = mask_file
  )
}

simulation_manifest <- bind_rows(manifest) %>%
  arrange(mechanism, missing_rate, repetition)

readr::write_csv(simulation_manifest, file.path(output_root, "simulation_manifest.csv"), na = "")

message("Simulation complete.")
message("Reference complete-case rows: ", nrow(reference_data))
message("Generated datasets: ", nrow(sim_grid))
message("Repetitions per setting: ", n_reps)
message("Outputs written under: ", output_root, "/")
