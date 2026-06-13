# ============================================================
#  02_descriptive.R
#  Table 1: Demographic and Clinical Characteristics
#  JAMA Style — CHD Diaphragm Paralysis Recovery Study
#
#  Input : df（01_data_cleaning.R で作成）
#  Output: Table1_CHD_DP_JAMA.docx
#
#  実行方法:
#    source("R/01_data_cleaning.R")
#    source("R/02_descriptive.R")
# ============================================================

library(dplyr)
library(gtsummary)
library(flextable)
library(officer)

# 01_data_cleaning.R が未実行の場合は自動実行
if (!exists("df")) {
  source("R/01_data_cleaning.R")
}


# ============================================================
# 1. 変数変換・ファクター化
# ============================================================
df_tbl <- df %>%
  mutate(
    # age_days → months
    age_months = age_days / 30.4375,

    sex = factor(sex,
                 levels = c("M", "F"),
                 labels = c("Male", "Female")),

    premature = factor(premature,
                       levels = c(1, 0),
                       labels = c("Yes", "No")),

    rachs1_category = factor(rachs1_category,
                             levels = 1:6,
                             labels = paste("Category", 1:6)),

    single_ventricle = factor(single_ventricle,
                              levels = c("Single_ventricle", "Biventricular"),
                              labels = c("Single ventricle", "Biventricular")),

    prior_cardiac_surgery = factor(prior_cardiac_surgery,
                                   levels = c(1, 0),
                                   labels = c("Yes", "No")),

    chromosomal_anomaly = factor(chromosomal_anomaly,
                                 levels = c(1, 0),
                                 labels = c("Yes", "No")),

    syndrome = factor(syndrome,
                      levels = c(1, 0),
                      labels = c("Yes", "No")),

    noncardiac_anomaly = factor(noncardiac_anomaly,
                                levels = c(1, 0),
                                labels = c("Yes", "No")),

    circulatory_arrest = factor(circulatory_arrest,
                                levels = c(1, 0),
                                labels = c("Yes", "No")),

    paralysis_side = factor(paralysis_side,
                            levels = c("Left", "Right", "Bilateral")),

    imaging_modality = factor(imaging_modality,
                              levels = c("Ultrasound", "Fluoroscopy", "Other")),

    paralysis_severity = factor(paralysis_severity,
                                levels = c("Weak", "Immobile", "Paradoxical")),

    plication_performed = factor(plication_performed,
                                 levels = c(1, 0),
                                 labels = c("Yes", "No")),

    recovery_confirmed = factor(recovery_confirmed,
                                levels = c(1, 0),
                                labels = c("Recovered", "Censored"))
  )


# ============================================================
# 2. Table 1 作成
# ============================================================
tbl <- df_tbl %>%
  select(
    # Patient Demographics
    age_months,
    sex,
    weight_kg,
    premature,

    # Congenital Heart Disease
    rachs1_category,
    single_ventricle,
    prior_cardiac_surgery,

    # Genetic and Morphologic Background
    chromosomal_anomaly,
    syndrome,
    noncardiac_anomaly,

    # Surgical Details
    cpb_time_min,
    xclamp_time_min,

    # Diaphragm Paralysis Characteristics
    paralysis_side,
    paralysis_severity,
    days_to_dp_diagnosis,

    # Intervention
    plication_performed,
    days_to_plication,

    # Follow-up and Outcome
    censoring_interval_days,
    followup_days,
    recovery_confirmed
  ) %>%
  tbl_summary(
    label = list(
      age_months 　　　　　　　~ "Age at surgery, mo\u00b2",
      sex                      ~ "Sex",
      weight_kg                ~ "Weight at surgery, kg",
      premature                ~ "Premature birth (<37 weeks)",
      rachs1_category          ~ "RACHS-1 category",
      single_ventricle         ~ "Single ventricle morphology",
      prior_cardiac_surgery    ~ "Prior cardiac surgery",
      chromosomal_anomaly      ~ "Any chromosomal anomaly",
      syndrome                 ~ "Any syndrome",
      noncardiac_anomaly       ~ "Any noncardiac anatomic abnormality",
      cpb_time_min             ~ "Cardiopulmonary bypass time, min",
      xclamp_time_min          ~ "Aortic crossclamp time, min",
      paralysis_side           ~ "Side of diaphragm paralysis",
      paralysis_severity       ~ "Severity of diaphragm paralysis",
      days_to_dp_diagnosis     ~ "Time from surgery to DP diagnosis, d",
      plication_performed      ~ "Diaphragmatic plication performed",
      days_to_plication 　　　 ~ "Time from DP diagnosis to diaphragmatic plication, d\u00b3",
      censoring_interval_days  ~ "Assessment interval (last abnormal to first normal), d\u2074",
      followup_days            ~ "Follow-up duration, d",
      recovery_confirmed       ~ "Recovery status"
    ),
    statistic = list(
      all_continuous()  ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      all_continuous()  ~ 1,
      all_categorical() ~ c(0, 1)
    ),
    missing = "no"
  ) %>%
  modify_header(
    label  ~ "**Characteristic**",
    stat_0 ~ "**Patients (N = {N})**"
  ) %>%
  modify_footnote(
    stat_0 ~ "Data are presented as median (IQR) or No. (%) unless otherwise indicated."
  ) %>%
  bold_labels() %>%
  modify_spanning_header(everything() ~ NA)

print(tbl)


# ============================================================
# 3. Word出力（JAMA投稿用）
# ============================================================
output_path <- "~/Desktop/Table1_CHD_DP_JAMA.docx"

tbl_flex <- tbl %>%
  as_flex_table() %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  bold(part = "header") %>%
  border_remove() %>%
  hline_top(part   = "header",
            border = fp_border(width = 1.5, color = "black")) %>%
  hline_bottom(part  = "header",
               border = fp_border(width = 1.0, color = "black")) %>%
  hline_bottom(part  = "body",
               border = fp_border(width = 1.5, color = "black")) %>%
  set_table_properties(layout = "autofit") %>%
  width(j = 1, width = 4.0) %>%
  width(j = 2, width = 2.0) %>%
  padding(padding.top    = 2, padding.bottom = 2,
          padding.left   = 4, padding.right  = 4,
          part = "all") %>%
  align(align = "left",   part = "all") %>%
  align(align = "center", j = 2, part = "all")

footnote_text <- paste0(
  "Abbreviations: CHD, congenital heart disease; DP, diaphragm paralysis; ",
  "IQR, interquartile range; RACHS-1, Risk Adjustment for Congenital Heart Surgery-1.\n",
  "Chromosomal anomaly includes numeric and structural chromosomal abnormalities ",
  "(e.g., trisomies, 22q11.2 deletion). Syndrome includes recognized clinical syndromes ",
  "without chromosomal abnormality (e.g., CHARGE, Noonan, VACTERL). ",
  "Categories are mutually exclusive.\n",
  "\u00b2Age recorded in days and converted to months (days \u00f7 30.4375) for presentation.\n",
  "\u00b3Among 46 patients who underwent diaphragmatic plication.\n",
  "\u2074Assessment interval is defined as the interval between the last imaging ",
  "at which diaphragm motion was abnormal and the first imaging at which recovery ",
  "was confirmed; calculated only among 84 patients with confirmed recovery. ",
  "Right-censored patients (n = 26) are excluded.\n",
  "Severity of diaphragm paralysis was graded as follows: weak, reduced but present ",
  "diaphragmatic motion; immobile, absent diaphragmatic motion; paradoxical, paradoxical ",
  "diaphragmatic motion (most severe).\n",
  "Assessment interval reflects the actual inter-assessment interval width, ",
  "which defines the censoring intervals used in the Turnbull nonparametric ",
  "maximum likelihood estimator for DP recovery time analysis."
)

doc <- read_docx() %>%
  body_add_par(
    paste0("Table 1. Demographic and Clinical Characteristics of Patients ",
           "With Diaphragm Paralysis After Congenital Heart Surgery"),
    style = "heading 1"
  ) %>%
  body_add_flextable(tbl_flex) %>%
  body_add_par(footnote_text, style = "Normal")

print(doc, target = output_path)
cat("\n\u2713 出力完了:", output_path, "\n")
