# ============================================================
#  07c_sensitivity_SA4.R
#  Supplemental Figure 4: SA4 Cox time-dependent plication
#  （immortal time bias補正の前後比較、Forest plot）
#
#  Input : df（01_data_cleaning.R で作成）
#  Output: SuppFig_SA4_Cox_plication.pdf / .tiff
#
#  注記: これはtime-varying outcome graph（生存曲線）ではなく
#  単一のハザード比を比較するForest plotのため、JTCVSの
#  5要件（年単位X軸・number at risk等）は適用対象外。
#  一般的な図表要件（カラー必須・個別ファイル・凡例埋め込み
#  禁止等）のみ満たせばよい。
#
#  [修正した点]
#  元のスクリプトでは labs() が + チェーンから独立した行として
#  書かれており、fig_sa4に反映されていなかった（X軸ラベルが
#  実際には表示されないバグ）。今回 labs() を ggplot の + チェーン
#  内に統合し、確実に反映されるよう修正した。
# ============================================================

library(survival)
library(ggplot2)
library(dplyr)

set.seed(1234)

if (!exists("df")) source("R/01_data_cleaning.R")

# ============================================================
# 1. データ準備（time-dependent covariate形式に変換）
# ============================================================
df_cox <- df %>%
  mutate(
    time_event = case_when(
      !is.na(right_time) ~ (left_time + right_time) / 2,
      is.na(right_time)  ~ followup_days
    ),
    event = as.integer(!is.na(right_time))
  )

# Naive Cox（immortal time biasあり・参照用）
fit_cox_naive <- coxph(
  Surv(time_event, event) ~ plication_performed,
  data = df_cox, ties = "efron"
)

# time-dependent covariate形式に変換
df_td <- df_cox %>%
  rowwise() %>%
  do({
    row <- .
    if (row$plication_performed == 1 &&
        !is.na(row$days_to_plication) &&
        row$days_to_plication < row$time_event) {
      bind_rows(
        data.frame(patient_id = row$patient_id, tstart = 0,
                   tstop = row$days_to_plication,
                   event_td = 0, plication_td = 0),
        data.frame(patient_id = row$patient_id, tstart = row$days_to_plication,
                   tstop = row$time_event,
                   event_td = row$event, plication_td = 1)
      )
    } else {
      data.frame(patient_id = row$patient_id, tstart = 0,
                 tstop = row$time_event,
                 event_td = row$event, plication_td = 0)
    }
  }) %>%
  ungroup()

# Time-dependent Cox（immortal time bias補正）
fit_cox_td <- coxph(
  Surv(tstart, tstop, event_td) ~ plication_td,
  data = df_td, ties = "efron"
)

hr_naive <- exp(coef(fit_cox_naive)["plication_performed"])
ci_naive <- exp(confint(fit_cox_naive)["plication_performed", ])
hr_td    <- exp(coef(fit_cox_td)["plication_td"])
ci_td    <- exp(confint(fit_cox_td)["plication_td", ])
ph_p     <- cox.zph(fit_cox_td)$table["plication_td", "p"]

cat("=== SA4 Cox sensitivity analysis results ===\n")
cat(sprintf("Naive Cox HR      : %.2f (95%% CI: %.2f\u2013%.2f)\n",
            hr_naive, ci_naive[1], ci_naive[2]))
cat(sprintf("Time-dependent HR : %.2f (95%% CI: %.2f\u2013%.2f)\n",
            hr_td, ci_td[1], ci_td[2]))
cat(sprintf("PH assumption p   : %.3f\n\n", ph_p))

# ============================================================
# 2. Forest plot（labs()のバグを修正）
# ============================================================
forest_df <- data.frame(
  Analysis = c("Naive Cox\n(immortal time bias present)",
               "Time-dependent Cox\n(immortal time bias corrected)"),
  HR    = c(hr_naive, hr_td),
  Lower = c(ci_naive[1], ci_td[1]),
  Upper = c(ci_naive[2], ci_td[2]),
  Bias  = c("Biased", "Corrected")
)
# 表示順を固定（Naiveが上、Time-dependentが下）
forest_df$Analysis <- factor(forest_df$Analysis, levels = rev(forest_df$Analysis))

fig_sa4 <- ggplot(forest_df,
                  aes(x = HR, y = Analysis, color = Bias)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "grey50", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper),
                 height = 0.15, linewidth = 0.9) +
  geom_point(size = 4, shape = 18) +
  geom_text(aes(label = sprintf("HR %.2f\n(95%% CI %.2f\u2013%.2f)", HR, Lower, Upper)),
            hjust = -0.12, size = 3.0, family = "sans") +
  scale_color_manual(values = c("Biased" = "#C0392B", "Corrected" = "#1F3864")) +
  scale_x_continuous(
    name = "Hazard ratio (log scale)\n[HR > 1: faster recovery with plication]",
    trans = "log10"
  ) +
  labs(y = NULL) +
  theme_classic(base_family = "sans", base_size = 11) +
  theme(
    legend.position    = "none",
    axis.title.x       = element_text(size = 10),
    axis.text.x        = element_text(size = 9, color = "grey20"),
    axis.text.y        = element_text(size = 9),
    panel.grid.major.x = element_line(color = "grey92"),
    plot.margin        = margin(10, 60, 10, 10)
  )

print(fig_sa4)

# ============================================================
# 3. 出力
# ============================================================
output_dir <- "~/Desktop"

ggsave(file.path(output_dir, "SuppFig_SA4_Cox_plication.pdf"),
       plot = fig_sa4, width = 7.0, height = 4.0, device = "pdf")
ggsave(file.path(output_dir, "SuppFig_SA4_Cox_plication.tiff"),
       plot = fig_sa4, width = 7.0, height = 4.0,
       dpi = 600, compression = "lzw")

cat("\u2713 SuppFig_SA4 完了（labs()バグ修正済み）\n")
