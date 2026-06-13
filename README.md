# Recovery from Diaphragm Paralysis After Congenital Heart Surgery

Analysis code for: **Recovery from Diaphragm Paralysis After Congenital Heart Surgery**

This repository contains the R scripts used for data cleaning, descriptive statistics,
primary and secondary survival analyses, sensitivity analyses, and figure generation
for the single-center cohort study.

## Scripts

| Script | Description |
|---|---|
| `00_setup.R` | Environment setup, packages, constants |
| `01_data_cleaning.R` | Data import and cleaning |
| `02_descriptive.R` | Table 1 (descriptive statistics) |
| `04_turnbull_npmle.R` | Figure 1 — Turnbull NPMLE with bootstrap 95% CI (primary analysis) |
| `06_bayesian_aft.R` | Figure 2 — Bayesian log-logistic AFT model (secondary analysis) |
| `07_sensitivity.R` | Supplemental Figures 2–4 — Sensitivity analyses SA1–SA4 |
| `08_supplemental_figure1_flow_diagram.R` | Supplemental Figure 1 — Study flow diagram |

> **Note:** Script numbering reflects the analysis pipeline as finalized.
> Intermediate exploratory scripts (03, 05) were consolidated into later
> scripts during development and are not included separately.

## Environment

- R version 4.6.0
- Required packages: `ggplot2`, `survival`, `flexsurv`, `brms` (cmdstanr backend),
  `icenReg`, `IPDfromKM`, `labelled`

## Data Availability

Patient-level data are not publicly available due to institutional data-sharing
restrictions. Digitized data derived from Smith et al. (2013), used to inform
the Bayesian prior, are described in the manuscript Methods.

## License

This project is licensed under the MIT License.
