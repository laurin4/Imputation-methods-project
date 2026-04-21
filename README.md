# Missing-Data Imputation Project (NHANES 2017-2018)

End-to-end R pipeline for a Data Science in Health project:

- Prepare NHANES analysis data
- Explore missingness structure
- Simulate artificial missingness (MCAR + MAR)
- Compare imputation methods
- Evaluate both imputation accuracy and downstream inference
- Produce publication-oriented tables/figures and final report material

## Project Structure

- `01_data_prep.R`  
  Download/merge NHANES XPT files, clean variables, restrict to adults, export analysis dataset.
- `02_missingness_exploration.R`  
  Missingness summaries + figures (pattern, variable %, pairwise heatmap).
- `03_generate_missingness.R`  
  Generate benchmarking datasets and mask indices under MCAR/MAR across missingness levels.
- `04_imputation_methods.R`  
  Run complete-case, mean/mode, kNN, missForest, and MICE imputation methods.
- `05_evaluation.R`  
  Compute RMSE/NRMSE/PFC and inference metrics (bias, coverage, CI width, etc.).
- `report.Rmd`  
  Final scientific report template integrating outputs.

## Current Default Runtime (CPU-friendly and balanced)

Defaults are intentionally limited so the pipeline is feasible on a laptop:

- Data prep: `MAX_ADULT_ROWS=2000`
- Simulation: `N_REPS=8`, `REF_MAX_ROWS=1000`
- Imputation defaults: `MISSFOREST_NTREE=40`, `MICE_M=8`, `MICE_MAXIT=5`
- Evaluation defaults: `MICE_M_EVAL=8`, `MICE_MAXIT_EVAL=5`

Important implementation detail:

- `04_imputation_methods.R` and `05_evaluation.R` use `sim_data/simulation_manifest.csv` as the source of truth, preventing stale-file mixing and ensuring both MCAR and MAR are included.

## Quick Start (recommended full run with current defaults)

```bash
Rscript "01_data_prep.R" && \
Rscript "02_missingness_exploration.R" && \
Rscript "03_generate_missingness.R" && \
Rscript "04_imputation_methods.R" && \
Rscript "05_evaluation.R"
```

## Optional Overrides

Use environment variables only if you want to override defaults:

- `01_data_prep.R`: `MAX_ADULT_ROWS`, `SAMPLE_ADULT_FRAC`
- `02_missingness_exploration.R`: `MAX_ROWS_EXPLORATION`
- `03_generate_missingness.R`: `N_REPS`, `REF_MAX_ROWS`
- `04_imputation_methods.R`: `MAX_FILES_PER_MECH`, `METHODS`, `KNN_K`, `MISSFOREST_NTREE`, `MICE_M`, `MICE_MAXIT`
- `05_evaluation.R`: `MAX_FILES_PER_MECH_EVAL`, `METHODS_EVAL`, `MICE_M_EVAL`, `MICE_MAXIT_EVAL`, `SKIP_PLOTS`

## Key Outputs

### Data and simulation artifacts
- `data/nhanes_analysis.csv`
- `sim_data/simulation_manifest.csv`
- `sim_data/datasets/*.csv`
- `sim_data/masks/*.csv`
- `imputed_data/imputation_manifest.csv`
- `imputed_data/<method>/*.csv`

### Evaluation tables (`outputs/`)
- `missingness_variable_summary.csv`
- `missingness_descriptive_table.csv`
- `pairwise_missingness_table.csv`
- `imputation_accuracy_summary.csv`
- `imputation_accuracy_by_repetition.csv`
- `model_performance_summary.csv`
- `model_estimates_by_repetition.csv`

### Figures (`figures/`)
- Missingness: pattern, variable %, pairwise heatmap
- Accuracy: numeric/categorical method comparison + boxplots
- Inference: bias, coverage, SE vs empirical SD, CI width

## Suggested Main-Text Figures

1. `missingness_by_variable_pct.png`
2. `missingness_pairwise_heatmap.png`
3. `accuracy_numeric_methods.png`
4. `accuracy_categorical_methods.png`
5. `inference_bias_methods.png`
6. `inference_coverage_methods.png`

## Notes

- If `UpSetR` is not installed, the optional UpSet missingness plot is skipped safely.
- `figures/`, `data/`, `sim_data/`, and `outputs/` are ignored via `.gitignore` to avoid committing large artifacts.