# =============================================================
# 00_setup.R
# Project: Recovery Time of Diaphragm Paralysis after
#          Pediatric Cardiac Surgery
# Authors: Yuichi Tabata, Shuhei Fujita, Shinichiro Tsurukawa
# Institution: Kyoto Prefectural University of Medicine
# =============================================================


# -------------------------------------------------------------
# 1. Random seed (for bootstrap and MCMC reproducibility)
# -------------------------------------------------------------
set.seed(20250101)


# -------------------------------------------------------------
# 2. Package loading
# -------------------------------------------------------------

# Data manipulation
library(tidyverse)    # dplyr, ggplot2, tidyr, readr, etc.
library(lubridate)    # Date arithmetic

# Survival analysis — frequentist
library(survival)     # Surv(), survfit(), survreg()
library(flexsurv)     # Flexible parametric survival models (AFT)
library(interval)     # Turnbull NPMLE (icfit)
library(icenReg)      # Alternative interval-censored regression

# Survival analysis — Bayesian
library(brms)         # Bayesian multilevel models via Stan
library(posterior)    # Posterior summaries (rhat, ess, draws)
library(bayesplot)    # MCMC diagnostics and trace plots

# Pseudo-IPD reconstruction (Guyot algorithm)
library(IPDfromKM)    # Optional: alternative to manual implementation

# Tables and reporting
library(gtsummary)    # Descriptive statistics tables (Table 1)
library(gt)           # Table formatting
library(knitr)        # kable()
library(kableExtra)   # Table styling

# Figures
library(ggsurvfit)    # Survival curve plotting
library(patchwork)    # Combining multiple ggplot2 panels
library(scales)       # Axis formatting

# Utilities
library(here)         # Relative file paths (here::here())
library(renv)         # Package version management


# -------------------------------------------------------------
# 3. Global ggplot2 theme
# -------------------------------------------------------------
theme_set(
  theme_classic(base_size = 12) +
    theme(
      axis.text    = element_text(color = "black"),
      legend.position = "bottom",
      plot.title   = element_text(face = "bold", size = 13)
    )
)


# -------------------------------------------------------------
# 4. Factor level definitions
#    (set globally to ensure consistent ordering across scripts)
# -------------------------------------------------------------

# Sex
sex_levels <- c("Male", "Female")

# Procedure type
procedure_levels <- c("Sternotomy", "Thoracotomy")

# Physiology
physiology_levels <- c("Biventricular", "Single ventricle")

# Hypothermia strategy
hypothermia_levels <- c("Normothermia", "Mild", "Moderate", "Deep")

# Paralysis side
paralysis_side_levels <- c("Left", "Right", "Bilateral")

# Paralysis type
paralysis_type_levels <- c("Immobile", "Paradoxical", "Weak")

# Imaging modality
modality_levels <- c("Ultrasonography", "Fluoroscopy")

# Recovery confirmed
recovery_levels <- c("Yes", "No")

# Plication
plication_levels <- c("Yes", "No")

# Follow-up status
followup_levels <- c("Recovered", "Censored", "Lost")


# -------------------------------------------------------------
# 5. File paths
# -------------------------------------------------------------
path_raw       <- here("data", "raw")
path_processed <- here("data", "processed")
path_external  <- here("data", "external")
path_figures   <- here("output", "figures")
path_tables    <- here("output", "tables")


# -------------------------------------------------------------
# 6. Study constants
# -------------------------------------------------------------
study_start    <- as.Date("2009-01-01")
study_end      <- as.Date("2024-12-31")
followup_cutoff <- as.Date("2025-12-31")
bootstrap_n    <- 2000   # Bootstrap iterations for Turnbull CI


# -------------------------------------------------------------
# 7. Session info (saved to docs/session_info.txt on run)
# -------------------------------------------------------------
sink(here("docs", "session_info.txt"))
print(sessionInfo())
sink()

message("00_setup.R complete — packages loaded, seed set.")
