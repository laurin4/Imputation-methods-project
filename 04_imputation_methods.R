#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# Apply imputation / missing-data handling methods to simulated datasets
# -----------------------------------------------------------------------------
# Goal:
#   For each dataset in sim_data/datasets/, apply:
#     1) complete-case analysis
#     2) mean/mode imputation
#     3) kNN imputation (VIM::kNN)
#     4) missForest (missForest::missForest)
#     5) multiple imputation (mice::mice, m = 20)
#
# Outputs:
#   - imputed_data/<method>/<method>_<original_file>.csv
#   - imputed_data/imputation_manifest.csv
#
# Standardized output format:
#   - Original simulation metadata columns are retained when present:
#       repetition, mechanism, missing_rate
#   - Added columns:
#       imputation_method, imputation_id
#   - For MICE: imputation_id = 1..m
#   - For single-output methods: imputation_id = 1
# -----------------------------------------------------------------------------

# ---- Reproducibility / defensive checks -------------------------------------
required_pkgs <- c("dplyr", "tidyr", "readr", "purrr", "stringr", "VIM", "missForest", "mice")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_pkgs, collapse = ", "),
    ". Install them before running this script, e.g.:\n",
    "install.packages(c(",
    paste(sprintf('"%s"', missing_pkgs), collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(stringr)
})

input_dir <- "sim_data/datasets"
output_root <- "imputed_data"
method_names <- c("complete_case", "mean_mode", "knn", "missforest", "mice")

expected_vars <- c(
  "hba1c", "age", "sex", "bmi", "systolic_bp", "smoking_status", "income_ratio"
)

binary_factor_vars <- c("sex", "smoking_status")
ordered_factor_vars <- character(0)
unordered_factor_vars <- c("sex", "smoking_status")
numeric_vars <- c("hba1c", "age", "bmi", "systolic_bp", "income_ratio")

for (m in method_names) {
  dir.create(file.path(output_root, m), recursive = TRUE, showWarnings = FALSE)
}

if (!dir.exists(input_dir)) {
  stop("Input directory not found: ", input_dir, ". Run 03_generate_missingness.R first.")
}

dataset_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(dataset_files) == 0) {
  stop("No simulation files found in ", input_dir, ".")
}

# ---- Runtime controls for faster/debug runs ---------------------------------
# Optional environment variables:
#   MAX_FILES           -> process only first N simulation files
#   METHODS             -> comma-separated subset, e.g. "complete_case,mean_mode,mice"
#   KNN_K               -> k for VIM::kNN (default 5)
#   MISSFOREST_NTREE    -> number of trees for missForest (default 100)
#   MICE_M              -> number of imputations for mice (default 20)
#   MICE_MAXIT          -> number of mice iterations (default 10)
max_files <- as.integer(Sys.getenv("MAX_FILES", unset = "0"))
if (!is.na(max_files) && max_files > 0) {
  dataset_files <- dataset_files[seq_len(min(max_files, length(dataset_files)))]
}

methods_env <- Sys.getenv("METHODS", unset = "")
if (nzchar(methods_env)) {
  requested_methods <- strsplit(methods_env, ",", fixed = TRUE)[[1]] %>%
    trimws() %>%
    unique()
  invalid_methods <- setdiff(requested_methods, method_names)
  if (length(invalid_methods) > 0) {
    stop(
      "Invalid method(s) in METHODS: ",
      paste(invalid_methods, collapse = ", "),
      ". Valid methods are: ",
      paste(method_names, collapse = ", ")
    )
  }
  method_names <- requested_methods
}

knn_k <- as.integer(Sys.getenv("KNN_K", unset = "5"))
missforest_ntree <- as.integer(Sys.getenv("MISSFOREST_NTREE", unset = "100"))
mice_m <- as.integer(Sys.getenv("MICE_M", unset = "20"))
mice_maxit <- as.integer(Sys.getenv("MICE_MAXIT", unset = "10"))

if (is.na(knn_k) || knn_k <= 0) stop("KNN_K must be a positive integer.")
if (is.na(missforest_ntree) || missforest_ntree <= 0) stop("MISSFOREST_NTREE must be a positive integer.")
if (is.na(mice_m) || mice_m <= 0) stop("MICE_M must be a positive integer.")
if (is.na(mice_maxit) || mice_maxit <= 0) stop("MICE_MAXIT must be a positive integer.")

# ---- Helper functions --------------------------------------------------------

assert_required_columns <- function(df, required_cols, context = "dataset") {
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required column(s) in ", context, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
}

safe_mode <- function(x) {
  x_non_na <- x[!is.na(x)]
  if (length(x_non_na) == 0) return(NA)
  tab <- table(x_non_na)
  names(tab)[which.max(tab)][1]
}

coerce_analysis_types <- function(df) {
  assert_required_columns(df, expected_vars, context = "input dataframe")

  out <- df

  for (v in numeric_vars) {
    out[[v]] <- suppressWarnings(as.numeric(out[[v]]))
  }

  if ("sex" %in% names(out)) {
    out$sex <- factor(out$sex, levels = c("Male", "Female"))
  }
  if ("smoking_status" %in% names(out)) {
    out$smoking_status <- factor(out$smoking_status, levels = c("No", "Yes"))
  }

  out
}

split_metadata_analysis <- function(df) {
  metadata_cols <- intersect(c("repetition", "mechanism", "missing_rate"), names(df))
  analysis <- df %>% select(all_of(expected_vars))
  metadata <- if (length(metadata_cols) > 0) df %>% select(all_of(metadata_cols)) else tibble::tibble()
  list(metadata = metadata, analysis = analysis)
}

bind_standardized_output <- function(imputed_analysis, original_df, method_name, imputation_id = 1L) {
  parts <- split_metadata_analysis(original_df)
  metadata <- parts$metadata

  if (nrow(imputed_analysis) != nrow(original_df)) {
    stop(
      "Row count mismatch after ", method_name, ": expected ",
      nrow(original_df), ", got ", nrow(imputed_analysis), "."
    )
  }

  out <- if (ncol(metadata) > 0) {
    bind_cols(metadata, imputed_analysis)
  } else {
    imputed_analysis
  }

  out %>%
    mutate(
      imputation_method = method_name,
      imputation_id = imputation_id,
      .after = last_col()
    )
}

# ---- Method 1: Complete-case analysis ---------------------------------------
impute_complete_case <- function(df) {
  typed <- coerce_analysis_types(df)
  parts <- split_metadata_analysis(typed)

  cc_idx <- complete.cases(parts$analysis)
  if (!any(cc_idx)) {
    stop("Complete-case analysis produced 0 rows.")
  }

  metadata_cc <- if (ncol(parts$metadata) > 0) parts$metadata[cc_idx, , drop = FALSE] else tibble::tibble()
  analysis_cc <- parts$analysis[cc_idx, , drop = FALSE]

  out <- if (ncol(metadata_cc) > 0) bind_cols(metadata_cc, analysis_cc) else analysis_cc
  out %>%
    mutate(
      imputation_method = "complete_case",
      imputation_id = 1L,
      .after = last_col()
    )
}

# ---- Method 2: Mean/mode imputation -----------------------------------------
impute_mean_mode <- function(df) {
  typed <- coerce_analysis_types(df)
  analysis <- split_metadata_analysis(typed)$analysis

  # Numeric: mean imputation
  for (v in numeric_vars) {
    if (!v %in% names(analysis)) next
    mean_val <- mean(analysis[[v]], na.rm = TRUE)
    if (is.nan(mean_val)) {
      stop("Cannot mean-impute variable ", v, " because all values are NA.")
    }
    analysis[[v]][is.na(analysis[[v]])] <- mean_val
  }

  # Factors: mode imputation
  for (v in union(binary_factor_vars, setdiff(unordered_factor_vars, binary_factor_vars))) {
    if (!v %in% names(analysis)) next
    mode_val <- safe_mode(analysis[[v]])
    if (is.na(mode_val)) {
      stop("Cannot mode-impute variable ", v, " because all values are NA.")
    }
    analysis[[v]][is.na(analysis[[v]])] <- mode_val
    analysis[[v]] <- droplevels(factor(analysis[[v]], levels = levels(typed[[v]])))
  }

  bind_standardized_output(
    imputed_analysis = analysis,
    original_df = typed,
    method_name = "mean_mode",
    imputation_id = 1L
  )
}

# ---- Method 3: kNN imputation (VIM::kNN) ------------------------------------
impute_knn <- function(df, k = 5, seed = 1234) {
  set.seed(seed)
  typed <- coerce_analysis_types(df)
  analysis <- split_metadata_analysis(typed)$analysis %>% as.data.frame()

  imputed <- tryCatch(
    VIM::kNN(
      data = analysis,
      variable = names(analysis),
      k = k,
      imp_var = FALSE
    ),
    error = function(e) {
      stop("VIM::kNN failed: ", e$message)
    }
  )

  bind_standardized_output(
    imputed_analysis = as_tibble(imputed),
    original_df = typed,
    method_name = "knn",
    imputation_id = 1L
  )
}

# ---- Method 4: missForest ----------------------------------------------------
impute_missforest <- function(df, seed = 5678, ntree = 100) {
  set.seed(seed)
  typed <- coerce_analysis_types(df)
  analysis <- split_metadata_analysis(typed)$analysis %>% as.data.frame()

  mf_out <- tryCatch(
    missForest::missForest(
      xmis = analysis,
      ntree = ntree,
      verbose = FALSE
    ),
    error = function(e) {
      stop("missForest::missForest failed: ", e$message)
    }
  )

  imputed <- as_tibble(mf_out$ximp)

  bind_standardized_output(
    imputed_analysis = imputed,
    original_df = typed,
    method_name = "missforest",
    imputation_id = 1L
  )
}

# ---- Method 5: Multiple imputation (mice) -----------------------------------
build_mice_method_vector <- function(df_analysis) {
  meth <- mice::make.method(df_analysis)

  for (v in names(df_analysis)) {
    x <- df_analysis[[v]]

    if (is.numeric(x)) {
      meth[v] <- "pmm"
    } else if (is.factor(x) && nlevels(x) == 2) {
      meth[v] <- "logreg"
    } else if (is.factor(x) && is.ordered(x)) {
      meth[v] <- "polr"
    } else if (is.factor(x)) {
      meth[v] <- "polyreg"
    } else {
      stop("Unsupported variable type for MICE in variable '", v, "'.")
    }
  }

  meth
}

impute_mice <- function(df, m = 20, maxit = 10, seed = 91011) {
  typed <- coerce_analysis_types(df)
  analysis <- split_metadata_analysis(typed)$analysis

  # Apply explicit factor ordering rules when relevant.
  for (v in ordered_factor_vars) {
    if (v %in% names(analysis)) analysis[[v]] <- as.ordered(analysis[[v]])
  }
  for (v in setdiff(unordered_factor_vars, ordered_factor_vars)) {
    if (v %in% names(analysis) && is.factor(analysis[[v]])) {
      analysis[[v]] <- factor(analysis[[v]], ordered = FALSE)
    }
  }

  mice_methods <- build_mice_method_vector(analysis)
  pred <- mice::make.predictorMatrix(analysis)
  diag(pred) <- 0

  mice_fit <- tryCatch(
    mice::mice(
      data = analysis,
      m = m,
      method = mice_methods,
      predictorMatrix = pred,
      maxit = maxit,
      seed = seed,
      printFlag = FALSE
    ),
    error = function(e) {
      stop("mice::mice failed: ", e$message)
    }
  )

  completed_long <- mice::complete(mice_fit, action = "long", include = FALSE) %>%
    as_tibble() %>%
    rename(imputation_id = ".imp") %>%
    select(all_of("imputation_id"), all_of(names(analysis)))

  # Add simulation metadata for each completed dataset.
  parts <- split_metadata_analysis(typed)
  metadata <- parts$metadata

  if (ncol(metadata) > 0) {
    completed_long <- completed_long %>%
      group_by(.data$imputation_id) %>%
      group_modify(~bind_cols(metadata, .x)) %>%
      ungroup()
  }

  completed_long %>%
    mutate(imputation_method = "mice", .after = last_col()) %>%
    select(
      any_of(c("repetition", "mechanism", "missing_rate")),
      all_of(expected_vars),
      all_of(c("imputation_method", "imputation_id"))
    )
}

# ---- Apply all methods to one dataset ---------------------------------------
apply_all_methods <- function(df, seed_offset = 0L) {
  out <- list()
  if ("complete_case" %in% method_names) {
    out$complete_case <- impute_complete_case(df)
  }
  if ("mean_mode" %in% method_names) {
    out$mean_mode <- impute_mean_mode(df)
  }
  if ("knn" %in% method_names) {
    out$knn <- impute_knn(df, k = knn_k, seed = 2000 + seed_offset)
  }
  if ("missforest" %in% method_names) {
    out$missforest <- impute_missforest(df, ntree = missforest_ntree, seed = 3000 + seed_offset)
  }
  if ("mice" %in% method_names) {
    out$mice <- impute_mice(df, m = mice_m, maxit = mice_maxit, seed = 4000 + seed_offset)
  }
  out
}

# ---- Main loop ---------------------------------------------------------------
manifest_entries <- vector("list", length(dataset_files) * length(method_names))
entry_idx <- 0L

for (file_i in seq_along(dataset_files)) {
  file_path <- dataset_files[file_i]
  file_name <- basename(file_path)

  df_in <- tryCatch(
    readr::read_csv(file_path, show_col_types = FALSE),
    error = function(e) stop("Failed reading simulation file ", file_name, ": ", e$message)
  )

  assert_required_columns(df_in, expected_vars, context = file_name)

  message(
    "Processing file ", file_i, "/", length(dataset_files), ": ", file_name,
    " | methods: ", paste(method_names, collapse = ", ")
  )

  method_results <- apply_all_methods(df_in, seed_offset = file_i)

  for (m in method_names) {
    out_df <- method_results[[m]]
    out_file <- file.path(output_root, m, paste0(m, "_", file_name))
    readr::write_csv(out_df, out_file, na = "")

    entry_idx <- entry_idx + 1L
    manifest_entries[[entry_idx]] <- tibble::tibble(
      source_file = file_name,
      method = m,
      output_file = out_file,
      n_rows = nrow(out_df),
      n_cols = ncol(out_df)
    )
  }
}

imputation_manifest <- bind_rows(manifest_entries)
readr::write_csv(imputation_manifest, file.path(output_root, "imputation_manifest.csv"), na = "")

message("Imputation complete.")
message("Processed simulation files: ", length(dataset_files))
message("Methods run: ", paste(method_names, collapse = ", "))
message(
  "Parameters used: KNN_K=", knn_k,
  ", MISSFOREST_NTREE=", missforest_ntree,
  ", MICE_M=", mice_m,
  ", MICE_MAXIT=", mice_maxit
)
message("Saved outputs in: ", output_root, "/")
