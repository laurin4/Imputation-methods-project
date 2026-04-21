#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# Evaluate imputation methods:
#   1) Imputation accuracy on artificially masked cells
#   2) Downstream inference performance in linear regression
# -----------------------------------------------------------------------------
# Model:
#   hba1c ~ age + sex + bmi + systolic_bp + smoking_status + income_ratio
#
# Inputs:
#   - data/nhanes_analysis.csv
#   - sim_data/datasets/
#   - sim_data/masks/
#   - imputed_data/<method>/
#
# Outputs:
#   - outputs/imputation_accuracy_by_repetition.csv
#   - outputs/imputation_accuracy_summary.csv
#   - outputs/model_estimates_by_repetition.csv
#   - outputs/model_performance_summary.csv
#   - figures/accuracy_numeric_methods.png
#   - figures/accuracy_categorical_methods.png
#   - figures/inference_bias_methods.png
#   - figures/inference_coverage_methods.png
#   - figures/accuracy_numeric_boxplots.png
#   - figures/accuracy_categorical_boxplots.png
#   - figures/inference_se_vs_empiricalsd.png
#   - figures/inference_ciwidth_methods.png
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(mice)
  library(broom)
})

# ---- Defensive checks --------------------------------------------------------
required_pkgs <- c("dplyr", "tidyr", "readr", "purrr", "stringr", "ggplot2", "mice", "broom")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_pkgs, collapse = ", "),
    ". Install them before running this script."
  )
}

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

if (!file.exists("data/nhanes_analysis.csv")) {
  stop("Missing input file: data/nhanes_analysis.csv. Run 01_data_prep.R first.")
}
if (!dir.exists("sim_data/datasets") || !dir.exists("sim_data/masks")) {
  stop("Missing sim_data inputs. Run 03_generate_missingness.R first.")
}
if (!dir.exists("imputed_data")) {
  stop("Missing imputed_data directory. Run 04_imputation_methods.R first.")
}

# ---- Runtime controls for balanced CPU usage --------------------------------
# Optional environment variables:
#   FAST_MODE        -> "1" enables practical defaults for quicker evaluation
#   MAX_FILES_EVAL   -> evaluate only first N simulation files
#   METHODS_EVAL     -> comma-separated subset of methods to evaluate
#   MICE_M_EVAL      -> m for pooled mice inference in evaluation
#   MICE_MAXIT_EVAL  -> maxit for pooled mice inference in evaluation
#   SKIP_PLOTS       -> "1" skips plot generation (tables only)
fast_mode <- identical(Sys.getenv("FAST_MODE", unset = "0"), "1")
max_files_eval <- as.integer(Sys.getenv("MAX_FILES_EVAL", unset = if (fast_mode) "20" else "0"))
mice_m_eval <- as.integer(Sys.getenv("MICE_M_EVAL", unset = if (fast_mode) "5" else "20"))
mice_maxit_eval <- as.integer(Sys.getenv("MICE_MAXIT_EVAL", unset = if (fast_mode) "5" else "10"))
skip_plots <- identical(Sys.getenv("SKIP_PLOTS", unset = "0"), "1")

if (is.na(max_files_eval) || max_files_eval < 0) stop("MAX_FILES_EVAL must be a non-negative integer.")
if (is.na(mice_m_eval) || mice_m_eval <= 0) stop("MICE_M_EVAL must be a positive integer.")
if (is.na(mice_maxit_eval) || mice_maxit_eval <= 0) stop("MICE_MAXIT_EVAL must be a positive integer.")

# ---- Configuration -----------------------------------------------------------
analysis_vars <- c("hba1c", "age", "sex", "bmi", "systolic_bp", "smoking_status", "income_ratio")
predictor_vars <- c("age", "sex", "bmi", "systolic_bp", "smoking_status", "income_ratio")
numeric_vars <- c("age", "bmi", "systolic_bp", "income_ratio")
categorical_vars <- c("sex", "smoking_status")
methods <- c("complete_case", "mean_mode", "knn", "missforest", "mice")

method_labels <- c(
  complete_case = "Complete-case",
  mean_mode = "Mean/mode",
  knn = "kNN",
  missforest = "missForest",
  mice = "MICE"
)
variable_labels <- c(
  age = "Age (years)",
  bmi = "Body Mass Index (kg/m^2)",
  systolic_bp = "Systolic blood pressure (mmHg)",
  income_ratio = "Income-to-poverty ratio",
  sex = "Sex",
  smoking_status = "Current smoking status"
)

methods_env <- Sys.getenv("METHODS_EVAL", unset = "")
if (nzchar(methods_env)) {
  requested_methods <- strsplit(methods_env, ",", fixed = TRUE)[[1]] %>%
    trimws() %>%
    unique()
  invalid_methods <- setdiff(requested_methods, methods)
  if (length(invalid_methods) > 0) {
    stop(
      "Invalid method(s) in METHODS_EVAL: ",
      paste(invalid_methods, collapse = ", "),
      ". Valid methods are: ",
      paste(methods, collapse = ", ")
    )
  }
  methods <- requested_methods
}

model_formula <- as.formula("hba1c ~ age + sex + bmi + systolic_bp + smoking_status + income_ratio")

# ---- Helpers -----------------------------------------------------------------
stop_if_missing_cols <- function(df, required_cols, context) {
  miss <- setdiff(required_cols, names(df))
  if (length(miss) > 0) {
    stop("Missing required column(s) in ", context, ": ", paste(miss, collapse = ", "))
  }
}

coerce_types <- function(df) {
  out <- df
  for (v in c("hba1c", "age", "bmi", "systolic_bp", "income_ratio")) {
    if (v %in% names(out)) out[[v]] <- suppressWarnings(as.numeric(out[[v]]))
  }
  if ("sex" %in% names(out)) out$sex <- factor(out$sex, levels = c("Male", "Female"))
  if ("smoking_status" %in% names(out)) out$smoking_status <- factor(out$smoking_status, levels = c("No", "Yes"))
  out
}

parse_sim_filename <- function(path) {
  nm <- basename(path)
  # sim_mcar_r10pct_rep001.csv
  m <- stringr::str_match(nm, "^sim_(mcar|mar)_r(\\d+)pct_rep(\\d{3})\\.csv$")
  if (any(is.na(m))) {
    stop("Filename does not match expected simulation pattern: ", nm)
  }
  tibble::tibble(
    mechanism = toupper(m[, 2]),
    missing_rate = as.numeric(m[, 3]) / 100,
    repetition = as.integer(m[, 4]),
    sim_file = nm
  )
}

to_mask_name <- function(sim_file) {
  stringr::str_replace(sim_file, "^sim_", "mask_")
}

to_imputed_name <- function(method, sim_file) {
  paste0(method, "_", sim_file)
}

calc_rmse <- function(truth, estimate) {
  sqrt(mean((estimate - truth)^2, na.rm = TRUE))
}

calc_nrmse <- function(truth, estimate) {
  denom <- stats::sd(truth, na.rm = TRUE)
  if (is.na(denom) || denom == 0) return(NA_real_)
  calc_rmse(truth, estimate) / denom
}

calc_pfc <- function(truth, estimate) {
  mean(truth != estimate, na.rm = TRUE)
}

extract_model_stats_lm <- function(df, method_name) {
  fit <- tryCatch(
    lm(model_formula, data = df),
    error = function(e) stop("Model fit failed for method ", method_name, ": ", e$message)
  )

  broom::tidy(fit, conf.int = TRUE, conf.level = 0.95) %>%
    transmute(
      term = .data$term,
      estimate = .data$estimate,
      std.error = .data$std.error,
      conf.low = .data$conf.low,
      conf.high = .data$conf.high
    )
}

extract_model_stats_mice <- function(df_incomplete, m = 20, maxit = 10, seed = 7777) {
  meth <- mice::make.method(df_incomplete)
  for (v in names(df_incomplete)) {
    x <- df_incomplete[[v]]
    if (is.numeric(x)) {
      meth[v] <- "pmm"
    } else if (is.factor(x) && nlevels(x) == 2) {
      meth[v] <- "logreg"
    } else if (is.factor(x) && is.ordered(x)) {
      meth[v] <- "polr"
    } else if (is.factor(x)) {
      meth[v] <- "polyreg"
    } else {
      stop("Unsupported type in mice inference for variable: ", v)
    }
  }

  pred <- mice::make.predictorMatrix(df_incomplete)
  diag(pred) <- 0

  imp <- tryCatch(
    mice::mice(
      data = df_incomplete,
      m = m,
      method = meth,
      predictorMatrix = pred,
      maxit = maxit,
      seed = seed,
      printFlag = FALSE
    ),
    error = function(e) stop("mice() failed during inference pooling: ", e$message)
  )

  fit <- with(imp, lm(hba1c ~ age + sex + bmi + systolic_bp + smoking_status + income_ratio))
  pooled <- pool(fit)

  summary(pooled, conf.int = TRUE, conf.level = 0.95) %>%
    as_tibble() %>%
    transmute(
      term = .data$term,
      estimate = .data$estimate,
      std.error = .data$std.error,
      conf.low = .data$`2.5 %`,
      conf.high = .data$`97.5 %`
    )
}

# ---- Ground-truth model coefficients -----------------------------------------
nhanes <- readr::read_csv("data/nhanes_analysis.csv", show_col_types = FALSE) %>%
  coerce_types() %>%
  select(all_of(analysis_vars)) %>%
  tidyr::drop_na()

if (nrow(nhanes) == 0) {
  stop("No complete cases found in data/nhanes_analysis.csv for truth model.")
}

truth_model <- lm(model_formula, data = nhanes)
truth_coef <- broom::tidy(truth_model) %>%
  select(term, true_estimate = estimate)

# ---- Discover simulation files ----------------------------------------------
sim_files <- list.files("sim_data/datasets", pattern = "^sim_.*\\.csv$", full.names = TRUE)
if (length(sim_files) == 0) stop("No simulation datasets found in sim_data/datasets.")
if (max_files_eval > 0) {
  sim_files <- sim_files[seq_len(min(max_files_eval, length(sim_files)))]
}

# ---- Main evaluation loop ----------------------------------------------------
accuracy_rows <- list()
model_rows <- list()
row_i <- 0L
model_i <- 0L

for (sim_path in sim_files) {
  sim_info <- parse_sim_filename(sim_path)
  sim_file <- sim_info$sim_file
  mechanism <- sim_info$mechanism
  missing_rate <- sim_info$missing_rate
  repetition <- sim_info$repetition
  message(
    "Evaluating ", basename(sim_path),
    " (", which(sim_files == sim_path), "/", length(sim_files), ")",
    " | methods: ", paste(methods, collapse = ", ")
  )

  sim_df <- readr::read_csv(sim_path, show_col_types = FALSE) %>%
    coerce_types()
  stop_if_missing_cols(sim_df, c("row_id", analysis_vars), sim_file)

  mask_path <- file.path("sim_data/masks", to_mask_name(sim_file))
  if (!file.exists(mask_path)) {
    stop("Missing mask file for ", sim_file, ": ", mask_path)
  }
  mask_df <- readr::read_csv(mask_path, show_col_types = FALSE) %>%
    mutate(true_value = as.character(true_value))
  stop_if_missing_cols(mask_df, c("row_id", "variable", "true_value"), basename(mask_path))

  for (method_name in methods) {
    imputed_path <- file.path("imputed_data", method_name, to_imputed_name(method_name, sim_file))
    if (!file.exists(imputed_path)) {
      stop("Missing imputed file for method ", method_name, ": ", imputed_path)
    }

    imp_raw <- readr::read_csv(imputed_path, show_col_types = FALSE)
    stop_if_missing_cols(imp_raw, c("imputation_method", "imputation_id", analysis_vars), basename(imputed_path))

    # -------------------- Imputation accuracy ---------------------------------
    # For non-complete-case methods, row order is preserved from simulation data.
    # For complete_case, rows are dropped, so masked-cell value matching by row_id
    # is not identifiable from saved outputs; we report NA metrics for accuracy.
    if (method_name == "complete_case") {
      for (v in predictor_vars) {
        metric_set <- if (v %in% numeric_vars) c("rmse", "nrmse") else "pfc"
        for (metric_name in metric_set) {
          row_i <- row_i + 1L
          accuracy_rows[[row_i]] <- tibble::tibble(
            mechanism = mechanism,
            missing_rate = missing_rate,
            repetition = repetition,
            method = method_name,
            variable = v,
            metric = metric_name,
            value = NA_real_,
            note = "not estimable: complete-case output drops rows"
          )
        }
      }
    } else {
      imp_df <- imp_raw %>%
        coerce_types() %>%
        mutate(row_id = row_number()) %>%
        select(row_id, imputation_id, all_of(analysis_vars))

      joined <- mask_df %>%
        left_join(
          imp_df %>% select(row_id, imputation_id, all_of(predictor_vars)),
          by = "row_id"
        )

      # Defensive check: non-CC methods should preserve all rows.
      if (any(is.na(joined$imputation_id))) {
        stop(
          "Row alignment failed for method ", method_name, " in ", sim_file,
          ". Imputed file appears to have dropped rows."
        )
      }

      # Variable-specific comparison table
      comp <- joined %>%
        rowwise() %>%
        mutate(imputed_value = as.character(get(variable))) %>%
        ungroup() %>%
        mutate(
          true_num = suppressWarnings(as.numeric(true_value)),
          imp_num = suppressWarnings(as.numeric(imputed_value))
        )

      # Compute metrics per imputation, then average across imputations.
      for (v in predictor_vars) {
        vv <- comp %>% filter(variable == v)
        if (nrow(vv) == 0) next

        if (v %in% numeric_vars) {
          by_imp <- vv %>%
            group_by(imputation_id) %>%
            summarise(
              rmse = calc_rmse(true_num, imp_num),
              nrmse = calc_nrmse(true_num, imp_num),
              .groups = "drop"
            )

          row_i <- row_i + 1L
          accuracy_rows[[row_i]] <- tibble::tibble(
            mechanism = mechanism,
            missing_rate = missing_rate,
            repetition = repetition,
            method = method_name,
            variable = v,
            metric = "rmse",
            value = mean(by_imp$rmse, na.rm = TRUE),
            note = NA_character_
          )
          row_i <- row_i + 1L
          accuracy_rows[[row_i]] <- tibble::tibble(
            mechanism = mechanism,
            missing_rate = missing_rate,
            repetition = repetition,
            method = method_name,
            variable = v,
            metric = "nrmse",
            value = mean(by_imp$nrmse, na.rm = TRUE),
            note = NA_character_
          )
        } else if (v %in% categorical_vars) {
          by_imp <- vv %>%
            group_by(imputation_id) %>%
            summarise(
              pfc = calc_pfc(true_value, imputed_value),
              .groups = "drop"
            )
          row_i <- row_i + 1L
          accuracy_rows[[row_i]] <- tibble::tibble(
            mechanism = mechanism,
            missing_rate = missing_rate,
            repetition = repetition,
            method = method_name,
            variable = v,
            metric = "pfc",
            value = mean(by_imp$pfc, na.rm = TRUE),
            note = NA_character_
          )
        }
      }
    }

    # -------------------- Downstream inference --------------------------------
    if (method_name == "mice") {
      # Use Rubin's rules with with() and pool() on the original incomplete data.
      model_tbl <- extract_model_stats_mice(
        df_incomplete = sim_df %>% select(all_of(analysis_vars)),
        m = mice_m_eval,
        maxit = mice_maxit_eval,
        seed = 900000 + repetition + as.integer(round(missing_rate * 1000))
      )
    } else {
      # Use imputed dataset for non-mice methods.
      # If multiple imputation_id rows exist unexpectedly, keep id==1.
      model_input <- imp_raw %>%
        filter(imputation_id == 1) %>%
        coerce_types() %>%
        select(all_of(analysis_vars))

      if (nrow(model_input) == 0) {
        stop("No rows available for model fitting in ", basename(imputed_path))
      }
      model_tbl <- extract_model_stats_lm(model_input, method_name)
    }

    model_tbl <- model_tbl %>%
      mutate(
        mechanism = mechanism,
        missing_rate = missing_rate,
        repetition = repetition,
        method = method_name,
        .before = 1
      )

    model_i <- model_i + 1L
    model_rows[[model_i]] <- model_tbl
  }
}

# ---- Bind and summarize ------------------------------------------------------
accuracy_by_rep <- bind_rows(accuracy_rows) %>%
  arrange(method, mechanism, missing_rate, repetition, variable, metric)

model_by_rep <- bind_rows(model_rows) %>%
  left_join(truth_coef, by = "term") %>%
  mutate(
    bias = estimate - true_estimate,
    ci_width = conf.high - conf.low,
    ci_cover = if_else(true_estimate >= conf.low & true_estimate <= conf.high, 1, 0, missing = 0)
  ) %>%
  arrange(method, mechanism, missing_rate, repetition, term)

accuracy_summary <- accuracy_by_rep %>%
  group_by(method, mechanism, missing_rate, variable, metric) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = stats::sd(value, na.rm = TRUE),
    n_repetitions = sum(!is.na(value)),
    .groups = "drop"
  ) %>%
  arrange(metric, variable, mechanism, missing_rate, method)

model_summary <- model_by_rep %>%
  group_by(method, mechanism, missing_rate, term) %>%
  summarise(
    mean_bias = mean(bias, na.rm = TRUE),
    empirical_sd = stats::sd(estimate, na.rm = TRUE),
    mean_model_se = mean(std.error, na.rm = TRUE),
    ci_coverage = mean(ci_cover, na.rm = TRUE),
    avg_ci_width = mean(ci_width, na.rm = TRUE),
    n_repetitions = n(),
    .groups = "drop"
  ) %>%
  arrange(term, mechanism, missing_rate, method)

# ---- Save tidy tables --------------------------------------------------------
readr::write_csv(accuracy_by_rep, "outputs/imputation_accuracy_by_repetition.csv", na = "")
readr::write_csv(accuracy_summary, "outputs/imputation_accuracy_summary.csv", na = "")
readr::write_csv(model_by_rep, "outputs/model_estimates_by_repetition.csv", na = "")
readr::write_csv(model_summary, "outputs/model_performance_summary.csv", na = "")

# ---- Plots -------------------------------------------------------------------
if (!skip_plots) {
  accuracy_summary_plot <- accuracy_summary %>%
    mutate(
      method = dplyr::recode(method, !!!method_labels),
      variable = dplyr::recode(variable, !!!variable_labels),
      mechanism = factor(mechanism, levels = c("MCAR", "MAR")),
      missing_rate_label = paste0(as.integer(100 * missing_rate), "%")
    )

  accuracy_by_rep_plot <- accuracy_by_rep %>%
    mutate(
      method = dplyr::recode(method, !!!method_labels),
      variable = dplyr::recode(variable, !!!variable_labels),
      mechanism = factor(mechanism, levels = c("MCAR", "MAR")),
      missing_rate_label = paste0(as.integer(100 * missing_rate), "%")
    )

  model_summary_plot <- model_summary %>%
    mutate(
      method = dplyr::recode(method, !!!method_labels),
      mechanism = factor(mechanism, levels = c("MCAR", "MAR")),
      missing_rate_label = paste0(as.integer(100 * missing_rate), "%")
    )

  acc_num_plot <- accuracy_summary %>%
    filter(metric %in% c("rmse", "nrmse"), !is.na(mean_value)) %>%
    mutate(
      method = dplyr::recode(method, !!!method_labels),
      variable = dplyr::recode(variable, !!!variable_labels),
      mechanism = factor(mechanism, levels = c("MCAR", "MAR")),
      missing_rate_label = paste0(as.integer(100 * missing_rate), "%")
    ) %>%
    ggplot(aes(x = missing_rate_label, y = mean_value, color = method, group = method)) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.8) +
    facet_grid(metric ~ variable + mechanism, scales = "free_y") +
    labs(
      title = "Numeric Imputation Accuracy Across Methods",
      subtitle = "Lower values indicate better imputation quality",
      x = "Artificial missingness level",
      y = "Mean error across repetitions",
      color = "Method"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "gray30"),
      strip.text = element_text(size = 9),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold")
    )

  ggsave(
    filename = "figures/accuracy_numeric_methods.png",
    plot = acc_num_plot,
    width = 14,
    height = 8,
    dpi = 300
  )

  acc_cat_plot <- accuracy_summary_plot %>%
    filter(metric == "pfc", !is.na(mean_value)) %>%
    ggplot(aes(x = missing_rate_label, y = mean_value, color = method, group = method)) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.8) +
    facet_grid(variable ~ mechanism, scales = "free_y") +
    labs(
      title = "Categorical Imputation Error (PFC) Across Methods",
      subtitle = "PFC = proportion falsely classified (lower is better)",
      x = "Artificial missingness level",
      y = "Mean PFC across repetitions",
      color = "Method"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "gray30"),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold")
    )

  ggsave(
    filename = "figures/accuracy_categorical_methods.png",
    plot = acc_cat_plot,
    width = 10,
    height = 6,
    dpi = 300
  )

  bias_plot <- model_summary_plot %>%
    filter(term != "(Intercept)") %>%
    ggplot(aes(x = missing_rate_label, y = mean_bias, color = method, group = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(size = 1.8) +
    geom_line(linewidth = 0.7) +
    facet_grid(term ~ mechanism, scales = "free_y") +
    labs(
      title = "Mean Coefficient Bias by Method",
      subtitle = "Bias = estimate - reference coefficient",
      x = "Artificial missingness level",
      y = "Mean bias",
      color = "Method"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "gray30"),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold")
    )

  ggsave(
    filename = "figures/inference_bias_methods.png",
    plot = bias_plot,
    width = 11,
    height = 9,
    dpi = 300
  )

  coverage_plot <- model_summary_plot %>%
    filter(term != "(Intercept)") %>%
    ggplot(aes(x = missing_rate_label, y = ci_coverage, color = method, group = method)) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "gray50") +
    geom_point(size = 1.8) +
    geom_line(linewidth = 0.7) +
    facet_grid(term ~ mechanism) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = "95% CI Coverage by Method",
      subtitle = "Dashed line indicates nominal 95% coverage target",
      x = "Artificial missingness level",
      y = "Coverage probability",
      color = "Method"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "gray30"),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold")
    )

  ggsave(
    filename = "figures/inference_coverage_methods.png",
    plot = coverage_plot,
    width = 11,
    height = 9,
    dpi = 300
  )

  # Additional plot: distribution of numeric metrics across repetitions
  acc_num_box <- accuracy_by_rep_plot %>%
    filter(metric %in% c("rmse", "nrmse"), !is.na(value), method != "Complete-case") %>%
    ggplot(aes(x = method, y = value, fill = method)) +
    geom_boxplot(outlier.alpha = 0.25, width = 0.75) +
    facet_grid(metric + mechanism ~ variable + missing_rate_label, scales = "free_y") +
    labs(
      title = "Distribution of Numeric Imputation Error Across Repetitions",
      subtitle = "Boxplots show variability and robustness of each method",
      x = "Method",
      y = "Error metric value",
      fill = "Method"
    ) +
    theme_minimal(base_size = 10.5) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "gray30"),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      axis.text.x = element_text(angle = 20, hjust = 1),
      legend.position = "none"
    )

  ggsave(
    filename = "figures/accuracy_numeric_boxplots.png",
    plot = acc_num_box,
    width = 14,
    height = 9,
    dpi = 300
  )

  # Additional plot: distribution of categorical PFC across repetitions
  acc_cat_box <- accuracy_by_rep_plot %>%
    filter(metric == "pfc", !is.na(value), method != "Complete-case") %>%
    ggplot(aes(x = method, y = value, fill = method)) +
    geom_boxplot(outlier.alpha = 0.25, width = 0.75) +
    facet_grid(mechanism ~ variable + missing_rate_label, scales = "free_y") +
    labs(
      title = "Distribution of Categorical Error (PFC) Across Repetitions",
      subtitle = "Lower PFC indicates better classification recovery",
      x = "Method",
      y = "PFC",
      fill = "Method"
    ) +
    theme_minimal(base_size = 10.5) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "gray30"),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      axis.text.x = element_text(angle = 20, hjust = 1),
      legend.position = "none"
    )

  ggsave(
    filename = "figures/accuracy_categorical_boxplots.png",
    plot = acc_cat_box,
    width = 12,
    height = 7,
    dpi = 300
  )

  # Additional plot: model-based SE vs empirical SD
  se_vs_empirical <- model_summary_plot %>%
    filter(term != "(Intercept)") %>%
    ggplot(aes(x = empirical_sd, y = mean_model_se, color = method)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(alpha = 0.9, size = 2.2) +
    facet_grid(mechanism ~ missing_rate_label) +
    labs(
      title = "Model-Based SE vs Empirical SD",
      subtitle = "Points near diagonal indicate calibrated uncertainty",
      x = "Empirical SD of estimates across repetitions",
      y = "Mean model-based standard error",
      color = "Method"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "gray30"),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold")
    )

  ggsave(
    filename = "figures/inference_se_vs_empiricalsd.png",
    plot = se_vs_empirical,
    width = 10.5,
    height = 6.8,
    dpi = 300
  )

  # Additional plot: average CI width by method
  ci_width_plot <- model_summary_plot %>%
    filter(term != "(Intercept)") %>%
    ggplot(aes(x = method, y = avg_ci_width, fill = method)) +
    geom_col(position = "dodge", width = 0.75) +
    facet_grid(term ~ mechanism + missing_rate_label, scales = "free_y") +
    labs(
      title = "Average 95% CI Width by Method",
      subtitle = "Wider intervals indicate more conservative uncertainty estimates",
      x = "Method",
      y = "Average CI width",
      fill = "Method"
    ) +
    theme_minimal(base_size = 10.5) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "gray30"),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      axis.text.x = element_text(angle = 20, hjust = 1),
      legend.position = "none"
    )

  ggsave(
    filename = "figures/inference_ciwidth_methods.png",
    plot = ci_width_plot,
    width = 14,
    height = 9,
    dpi = 300
  )
}

message("Evaluation complete.")
message("Files evaluated: ", length(sim_files))
message("Methods evaluated: ", paste(methods, collapse = ", "))
message("MICE pooling parameters: m=", mice_m_eval, ", maxit=", mice_maxit_eval)
if (skip_plots) message("Plot generation skipped (SKIP_PLOTS=1).")
message("Saved tables to outputs/ and plots to figures/.")
