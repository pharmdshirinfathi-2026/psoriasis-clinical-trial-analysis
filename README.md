# Psoriasis Clinical Trial Therapeutic-Class Analysis

This repository contains the R code and input metadata used to classify completed phase II–IV psoriasis and psoriatic arthritis clinical trials by therapeutic class and to generate the summary tables and figures reported in the associated publication.

## Repository contents

```text
.
├── Clinical_trial_analysis.R              # Main R analysis script
├── meta_data_psoriasis.csv                # Input clinical trial metadata
├── README.md                              # Project description and reproducibility guide
└── .gitignore                             # Optional: excludes generated outputs/cache files
```

After running the script, the analysis creates a date-stamped output folder:

```text
meta_data_psoriasis_YYYY-MM-DD/
├── tables/
├── figures/
└── QC_manual_classification/
```

## Study aim

The script analyzes ClinicalTrials.gov-derived metadata for completed phase II–IV psoriasis/psoriatic arthritis trials. Interventions are parsed from the input metadata and classified into therapeutic categories using a predefined dictionary and optional manual quality-control overrides.

The main therapeutic classes are:

| Code | Class |
|---:|---|
| 0 | Non-drug / non-relevant |
| 1 | Small molecule |
| 2 | Biologic |
| 3 | Conventional systemic / csDMARD / conventional topical |

Trial-level categories are then assigned based on the combination of intervention classes present in each trial.

## Input data

The input file must be named:

```text
meta_data_psoriasis.csv
```

Required columns:

```text
NCT Number
Study Title
Study URL
Brief Summary
Interventions
Primary Outcome Measures
Phases
Study Design
Completion Date
```

The currently provided input file contains ClinicalTrials.gov-style metadata with one row per trial and pipe-separated intervention records in the `Interventions` column.

## Software requirements

The script uses the native R pipe operator (`|>`), so R version 4.1 or newer is recommended.

Required R packages:

```r
install.packages(c(
  "dplyr",
  "tidyr",
  "stringr",
  "readr",
  "purrr",
  "tibble",
  "ggplot2",
  "lubridate"
))
```

The script also uses the base/grid package `grid`, which is included with R.

## How to run

1. Download or clone the repository.
2. Open R or RStudio in the repository folder.
3. Make sure `Clinical_trial_analysis.R` and `meta_data_psoriasis.csv` are in the same folder.
4. Run:

```r
source("Clinical_trial_analysis.R")
```

or from the terminal:

```bash
Rscript Clinical_trial_analysis.R
```

A new output folder will be created automatically.

## Outputs

### Tables

The `tables/` folder includes:

| File | Description |
|---|---|
| `01_intervention_level_classified.csv` | Intervention-level classification results |
| `02_intervention_matched_summary.csv` | Summary of dictionary-matched intervention names |
| `03_trial_level_classified.csv` | Trial-level classification table |
| `summary_01_trial_class_distribution.csv` | Trial-level therapeutic-class distribution |
| `summary_02_phase_by_trial_class.csv` | Trial class summarized by clinical phase |
| `summary_03_completion_year_by_trial_class.csv` | Trial class summarized by completion year |
| `summary_04_intervention_code_distribution.csv` | Intervention-record distribution by class code |
| `summary_05_trials_with_unmatched_items.csv` | Trials containing unmatched intervention records |

### Figures

The `figures/` folder includes publication-ready PDF and 600-dpi PNG versions of:

| Figure file prefix | Description |
|---|---|
| `Figure_1_trial_class_distribution` | Trial-level therapeutic-class distribution |
| `Figure_2_phase_by_trial_class` | Therapeutic-class distribution across phases |
| `Figure_3_completion_year_by_trial_class` | Temporal distribution by completion year |
| `Figure_4_intervention_code_distribution` | Intervention-level class-code distribution |

## Reproducibility notes

- The script uses a dictionary-based intervention classification approach.
- The dictionary is embedded directly in the R script for transparency.
- The output folder is date-stamped, so repeated runs on different dates will create separate output folders.
- The script should be run from the folder containing the input CSV file.

## Recommended citation
The current work has not been published yet.