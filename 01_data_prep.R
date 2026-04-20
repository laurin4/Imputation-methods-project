#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# NHANES 2017-2018 data preparation for missing-data imputation project
# -----------------------------------------------------------------------------
# This script:
#   1) Downloads selected NHANES XPT files directly from CDC
#   2) Merges files by SEQN
#   3) Restricts to adults (age >= 20)
#   4) Selects and recodes analysis variables
#   5) Builds a missingness summary
#   6) Saves analysis and summary outputs to CSV
# -----------------------------------------------------------------------------

# Reproducibility: fail early if required packages are missing
required_pkgs <- c("tidyverse", "haven")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_pkgs, collapse = ", "),
    ". Install them before running this script."
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(haven)
})

# ---- Runtime controls for balanced CPU usage --------------------------------
# Optional environment variables:
#   FAST_MODE           -> "1" enables practical defaults for quicker runs
#   MAX_ADULT_ROWS      -> keep only first N adult rows after preprocessing
#   SAMPLE_ADULT_FRAC   -> random sample fraction in (0, 1], applied before MAX
fast_mode <- identical(Sys.getenv("FAST_MODE", unset = "0"), "1")
max_adult_rows <- as.integer(Sys.getenv("MAX_ADULT_ROWS", unset = if (fast_mode) "2500" else "0"))
sample_adult_frac <- as.numeric(Sys.getenv("SAMPLE_ADULT_FRAC", unset = "1"))

if (is.na(max_adult_rows) || max_adult_rows < 0) {
  stop("MAX_ADULT_ROWS must be a non-negative integer.")
}
if (is.na(sample_adult_frac) || sample_adult_frac <= 0 || sample_adult_frac > 1) {
  stop("SAMPLE_ADULT_FRAC must be in the interval (0, 1].")
}

# Create output directories if they do not already exist
dir.create("data", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

# NHANES 2017-2018 source file URLs (CDC)
# Current CDC direct file endpoint for 2017-2018 cycle:
# https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/
base_url <- "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/"
xpt_urls <- list(
  DEMO_J = paste0(base_url, "DEMO_J.xpt"),
  BMX_J  = paste0(base_url, "BMX_J.xpt"),
  BPX_J  = paste0(base_url, "BPX_J.xpt"),
  GHB_J  = paste0(base_url, "GHB_J.xpt"),
  SMQ_J  = paste0(base_url, "SMQ_J.xpt"),
  INQ_J  = paste0(base_url, "INQ_J.xpt")
)

# Alternate endpoint pattern (legacy NHANES path), kept as fallback.
fallback_base_url <- "https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/"
xpt_fallback_urls <- list(
  DEMO_J = paste0(fallback_base_url, "DEMO_J.XPT"),
  BMX_J  = paste0(fallback_base_url, "BMX_J.XPT"),
  BPX_J  = paste0(fallback_base_url, "BPX_J.XPT"),
  GHB_J  = paste0(fallback_base_url, "GHB_J.XPT"),
  SMQ_J  = paste0(fallback_base_url, "SMQ_J.XPT"),
  INQ_J  = paste0(fallback_base_url, "INQ_J.XPT")
)

# Robust XPT reader with fallback URL and explicit diagnostics
read_nhanes_xpt <- function(url, table_name) {
  fallback_url <- xpt_fallback_urls[[table_name]]
  if (is.null(fallback_url)) {
    stop("No fallback URL configured for table: ", table_name)
  }

  # Retry wrapper to handle transient network issues
  read_with_retry <- function(target_url, attempts = 3, sleep_seconds = 1) {
    last_error <- NULL
    for (i in seq_len(attempts)) {
      out <- tryCatch(
        haven::read_xpt(target_url),
        error = function(e) {
          last_error <<- e
          NULL
        }
      )
      if (!is.null(out)) return(out)
      if (i < attempts) Sys.sleep(sleep_seconds)
    }
    stop(last_error)
  }

  primary <- tryCatch(
    read_with_retry(url, attempts = 2, sleep_seconds = 1),
    error = function(e) e
  )

  if (!inherits(primary, "error")) {
    return(primary)
  }

  message("Primary NHANES URL failed for ", table_name, ". Trying alternate CDC endpoint...")

  mirror <- tryCatch(
    read_with_retry(fallback_url, attempts = 2, sleep_seconds = 1),
    error = function(e) e
  )

  if (!inherits(mirror, "error")) {
    return(mirror)
  }

  stop(
    "Failed to read ", table_name, " from both NHANES sources.\n",
    "Primary URL: ", url, "\n",
    "Primary error: ", primary$message, "\n",
    "Fallback URL: ", fallback_url, "\n",
    "Fallback error: ", mirror$message, "\n",
    "If you are behind a proxy/firewall, try running from a network with direct HTTPS access."
  )
}

# Read and minimally subset each table before merging for efficiency
demo <- read_nhanes_xpt(xpt_urls$DEMO_J, "DEMO_J") %>%
  select(SEQN, RIDAGEYR, RIAGENDR)

bmx <- read_nhanes_xpt(xpt_urls$BMX_J, "BMX_J") %>%
  select(SEQN, BMXBMI)

bpx <- read_nhanes_xpt(xpt_urls$BPX_J, "BPX_J") %>%
  select(SEQN, starts_with("BPXSY"))

ghb <- read_nhanes_xpt(xpt_urls$GHB_J, "GHB_J") %>%
  select(SEQN, LBXGH)

smq <- read_nhanes_xpt(xpt_urls$SMQ_J, "SMQ_J") %>%
  select(SEQN, SMQ020)

inq_raw <- read_nhanes_xpt(xpt_urls$INQ_J, "INQ_J")

# Income-to-poverty ratio naming can vary by release.
income_ratio_var <- c("INDFMPIR", "INDFMMPI", "INDFMMPC")
income_ratio_var <- income_ratio_var[income_ratio_var %in% names(inq_raw)][1]
if (is.na(income_ratio_var) || length(income_ratio_var) == 0) {
  stop(
    "Could not find an income ratio variable in INQ_J. Checked: ",
    paste(c("INDFMPIR", "INDFMMPI", "INDFMMPC"), collapse = ", ")
  )
}

inq <- inq_raw %>%
  transmute(
    SEQN,
    INDFMPIR = .data[[income_ratio_var]]
  )

# Merge all selected tables by respondent ID
nhanes_merged <- list(demo, bmx, bpx, ghb, smq, inq) %>%
  reduce(full_join, by = "SEQN")

# Build analysis dataset
nhanes_analysis <- nhanes_merged %>%
  # Restrict to adults
  filter(RIDAGEYR >= 20) %>%
  # Derive average systolic BP across available measures (ignores NAs)
  mutate(
    systolic_bp = rowMeans(
      across(any_of(c("BPXSY1", "BPXSY2", "BPXSY3", "BPXSY4"))),
      na.rm = TRUE
    ),
    # If all systolic readings are missing, rowMeans returns NaN -> convert to NA
    systolic_bp = na_if(systolic_bp, NaN)
  ) %>%
  transmute(
    age = RIDAGEYR,
    sex = case_when(
      RIAGENDR == 1 ~ "Male",
      RIAGENDR == 2 ~ "Female",
      TRUE ~ NA_character_
    ),
    bmi = BMXBMI,
    hba1c = LBXGH,
    income_ratio = INDFMPIR,
    smoking_status = case_when(
      SMQ020 == 1 ~ "Yes",
      SMQ020 == 2 ~ "No",
      # NHANES special/missing codes (e.g., 7 = Refused, 9 = Don't know)
      SMQ020 %in% c(7, 9) ~ NA_character_,
      TRUE ~ NA_character_
    ),
    systolic_bp = systolic_bp
  ) %>%
  # Extra safety: remove known impossible/invalid sentinel values if encountered
  mutate(
    bmi = if_else(bmi <= 0, NA_real_, bmi),
    hba1c = if_else(hba1c <= 0, NA_real_, hba1c),
    income_ratio = if_else(income_ratio < 0, NA_real_, income_ratio),
    systolic_bp = if_else(systolic_bp <= 0, NA_real_, systolic_bp),
    sex = factor(sex, levels = c("Male", "Female")),
    smoking_status = factor(smoking_status, levels = c("No", "Yes"))
  )

# Optional downsampling for faster downstream processing
if (sample_adult_frac < 1) {
  set.seed(20260420)
  n_keep <- max(1L, floor(nrow(nhanes_analysis) * sample_adult_frac))
  nhanes_analysis <- nhanes_analysis %>%
    slice_sample(n = n_keep)
}

if (max_adult_rows > 0 && nrow(nhanes_analysis) > max_adult_rows) {
  nhanes_analysis <- nhanes_analysis %>%
    slice_head(n = max_adult_rows)
}

# Missingness summary for analysis variables
missingness_summary <- nhanes_analysis %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "n_missing"
  ) %>%
  mutate(
    pct_missing = round((n_missing / nrow(nhanes_analysis)) * 100, 2)
  ) %>%
  arrange(desc(pct_missing), variable)

# Save outputs
readr::write_csv(nhanes_analysis, "data/nhanes_analysis.csv", na = "")
readr::write_csv(missingness_summary, "outputs/missingness_summary.csv", na = "")

message("Data preparation complete.")
message("Rows in analysis dataset: ", nrow(nhanes_analysis))
message("Saved: data/nhanes_analysis.csv")
message("Saved: outputs/missingness_summary.csv")
