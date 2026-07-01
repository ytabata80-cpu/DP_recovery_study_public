# ============================================================
#  07_sensitivity.R
#  Pre-specified Sensitivity Analyses（最終版）
#
#  Input : df（01_data_cleaning.R で作成）
#  Output: SuppFig_Sensitivity_SA1_SA2.pdf / .tiff
#          SuppFig_SA3_KM_imputation.pdf / .tiff
#          SuppFig_SA4_Cox_plication.pdf / .tiff
#          ~/Desktop/ に保存
#
#  実行方法:
#    source("R/01_data_cleaning.R")
#    source("R/07_sensitivity.R")
#
#  感度解析一覧:
#    SA1: 両側麻痺例を全例除外
#    SA2: 180日以内右打ち切り除外
#    SA3: KM 3補完（left/mid/right endpoint）
#    SA4: Cox time-dependent plication（immortal time bias）
# ============================================================

library(survival)
library(ggplot2)
library(dplyr)

set.seed(1234)

if (!exists("df")) source("R/01_data_cleaning.R")

output_dir <- "~/Desktop"

cat("=== 主解析（参照）===\n")
cat("n =", nrow(df), "\n")
cat("回復確認:", sum(df$recovery_confirmed == 1, na.rm=TRUE), "\n")
cat("右側打ち切り:", sum(df$recovery_confirmed == 0, na.rm=TRUE), "\n\n")


# ============================================================
# 共通関数
# ============================================================
make_tb_df <- function(fit, label, t_max = 1095) {
  data.frame(time = c(0, fit$time), surv = c(1, fit$surv)) %>%
    mutate(cumrec = 1 - surv, group = label) %>%
    bind_rows(data.frame(
      time = t_max, surv = tail(.$surv, 1),
      cumrec = tail(.$cumrec, 1), group = label
    ))
}


# ============================================================
# SA1: 両側麻痺例を全例除外
# ============================================================
cat("SA1: 両側麻痺例除外\n")
df_sa1 <- df %>% filter(paralysis_side != "Bilateral")
fit_sa1 <- survfit(
  Surv(left_time, right_time, type = "interval2") ~ 1, data = df_sa1)
med_sa1 <- quantile(fit_sa1, probs = 0.5)
cat(sprintf("n=%d, Median=%.0f days\n\n", nrow(df_sa1), med_sa1$quantile))


# ============================================================
# SA2: 180日以内右打ち切り除外
# ============================================================
cat("SA2: 180日以内右打ち切り除外\n")
df_sa2 <- df %>%
  filter(recovery_confirmed == 1 |
         (recovery_confirmed == 0 & followup_days > 180))
fit_sa2 <- survfit(
  Surv(left_time, right_time, type = "interval2") ~ 1, data = df_sa2)
med_sa2 <- quantile(fit_sa2, probs = 0.5)
cat(sprintf("n=%d, Median=%.0f days\n\n", nrow(df_sa2), med_sa2$quantile))


# ============================================================
# SA1・SA2 Figure
# ============================================================
fit_main <- survfit(
  Surv(left_time, right_time, type = "interval2") ~ 1, data = df)

label_main <- "Primary analysis (n=110)"
label_sa1  <- "SA1: Exclude bilateral (n=101)"
label_sa2  <- paste0("SA2: Exclude early censored (n=", nrow(df_sa2), ")")

tb_sa12 <- bind_rows(
  make_tb_df(fit_main, label_main),
  make_tb_df(fit_sa1,  label_sa1),
  make_tb_df(fit_sa2,  label_sa2)
)

COL_MAIN <- "#1F3864"
COL_SA1  <- "#C0392B"
COL_SA2  <- "#2A9D8F"

fig_sa12 <- ggplot(tb_sa12,
                   aes(x = time, y = cumrec * 100,
                       color = group, linetype = group,
                       linewidth = group)) +
  geom_step(direction = "hv") +
  scale_color_manual(values = setNames(
    c(COL_MAIN, COL_SA1, COL_SA2),
    c(label_main, label_sa1, label_sa2)
  )) +
  scale_linetype_manual(values = setNames(
    c("solid", "dashed", "dotted"),
    c(label_main, label_sa1, label_sa2)
  )) +
  scale_linewidth_manual(values = setNames(
    c(1.2, 1.0, 1.0),
    c(label_main, label_sa1, label_sa2)
  ), guide = "none") +
  scale_x_continuous(
    name = "Time from diaphragm paralysis diagnosis (days)",
    limits = c(0, 1095), breaks = seq(0, 1095, by = 182),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    name = "Cumulative recovery probability (%)",
    limits = c(0, 100), breaks = seq(0, 100, by = 25),
    expand = expansion(mult = c(0, 0.02))
  ) +
  theme_classic(base_family = "sans", base_size = 11) +
theme(
  legend.position      = c(0.98, 0.10),
  legend.justification = c(1, 0),
  legend.title         = element_blank(),
  legend.text          = element_text(size = 8.5, family = "sans"),
  legend.key.width     = unit(1.5, "cm"),
  legend.background    = element_rect(fill = "white", color = "grey80",
                                      linewidth = 0.3),
  panel.grid           = element_blank(),
  plot.margin          = margin(10, 15, 10, 10)
)
# labs(title=, caption=) を削除

ggsave(file.path(output_dir, "SuppFig_Sensitivity_SA1_SA2.pdf"),
       plot = fig_sa12, width = 7.0, height = 5.5, device = "pdf")
ggsave(file.path(output_dir, "SuppFig_Sensitivity_SA1_SA2.tiff"),
       plot = fig_sa12, width = 7.0, height = 5.5,
       dpi = 600, compression = "lzw")
cat("SA1/SA2 Figure complete\n\n")


# ============================================================
# SA3: KM 3補完（left / mid / right endpoint）
# ============================================================
cat("SA3: KM 3補完\n")

df_sa3_left <- df %>%
  mutate(time_km = ifelse(!is.na(right_time), left_time, followup_days),
         event = recovery_confirmed)
fit_sa3_left  <- survfit(Surv(time_km, event) ~ 1, data = df_sa3_left)
med_sa3_left  <- quantile(fit_sa3_left, probs = 0.5)

df_sa3_mid <- df %>%
  mutate(time_km = ifelse(!is.na(right_time),
                          (left_time + right_time) / 2, followup_days),
         event = recovery_confirmed)
fit_sa3_mid  <- survfit(Surv(time_km, event) ~ 1, data = df_sa3_mid)
med_sa3_mid  <- quantile(fit_sa3_mid, probs = 0.5)

df_sa3_right <- df %>%
  mutate(time_km = ifelse(!is.na(right_time), right_time, followup_days),
         event = recovery_confirmed)
fit_sa3_right <- survfit(Surv(time_km, event) ~ 1, data = df_sa3_right)
med_sa3_right <- quantile(fit_sa3_right, probs = 0.5)

cat("Left  endpoint:", round(med_sa3_left$quantile),  "days\n")
cat("Mid   endpoint:", round(med_sa3_mid$quantile),   "days\n")
cat("Right endpoint:", round(med_sa3_right$quantile), "days\n\n")

# SA3 Figure
label_turnbull <- "Turnbull NPMLE (median=44d)"
label_left     <- "KM Left endpoint (median=16d)"
label_mid      <- "KM Mid endpoint (median=109d)"
label_right    <- "KM Right endpoint (median=205d)"

tb_sa3 <- bind_rows(
  make_tb_df(fit_main,      label_turnbull),
  make_tb_df(fit_sa3_left,  label_left),
  make_tb_df(fit_sa3_mid,   label_mid),
  make_tb_df(fit_sa3_right, label_right)
)

COL_TURNBULL <- "#1F3864"
COL_LEFT     <- "#C0392B"
COL_MID      <- "#E07B39"
COL_RIGHT    <- "#2A9D8F"

fig_sa3 <- ggplot(tb_sa3,
                  aes(x = time, y = cumrec * 100,
                      color = group, linetype = group,
                      linewidth = group)) +
  geom_step(direction = "hv") +
  geom_hline(yintercept = 50, linetype = "longdash",
             color = "grey40", linewidth = 0.5) +
  scale_color_manual(values = setNames(
    c(COL_TURNBULL, COL_LEFT, COL_MID, COL_RIGHT),
    c(label_turnbull, label_left, label_mid, label_right)
  )) +
  scale_linetype_manual(values = setNames(
    c("solid", "dashed", "dotted", "dotdash"),
    c(label_turnbull, label_left, label_mid, label_right)
  )) +
  scale_linewidth_manual(values = setNames(
    c(1.2, 0.9, 0.9, 0.9),
    c(label_turnbull, label_left, label_mid, label_right)
  ), guide = "none") +
  scale_x_continuous(
    name = "Time from diaphragm paralysis diagnosis (days)",
    limits = c(0, 1095), breaks = seq(0, 1095, by = 182),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    name = "Cumulative recovery probability (%)",
    limits = c(0, 100), breaks = seq(0, 100, by = 25),
    expand = expansion(mult = c(0, 0.02))
  ) +
  theme_classic(base_family = "sans", base_size = 11) +
theme(
  legend.position      = c(0.98, 0.10),
  legend.justification = c(1, 0),
  legend.title         = element_blank(),
  legend.text          = element_text(size = 8.0, family = "sans"),
  legend.key.width     = unit(1.5, "cm"),
  legend.background    = element_rect(fill = "white", color = "grey80",
                                      linewidth = 0.3),
  panel.grid           = element_blank(),
  plot.margin          = margin(10, 15, 10, 10)
)
# labs(title=, caption=) を削除

ggsave(file.path(output_dir, "SuppFig_SA3_KM_imputation.pdf"),
       plot = fig_sa3, width = 7.0, height = 5.5, device = "pdf")
ggsave(file.path(output_dir, "SuppFig_SA3_KM_imputation.tiff"),
       plot = fig_sa3, width = 7.0, height = 5.5,
       dpi = 600, compression = "lzw")
cat("SA3 Figure complete\n\n")


# ============================================================
# SA4: Cox time-dependent plication（immortal time bias補正）
# ============================================================
cat("SA4: Cox time-dependent plication\n")

df_cox <- df %>%
  mutate(
    time_event = case_when(
      !is.na(right_time) ~ (left_time + right_time) / 2,
      is.na(right_time)  ~ followup_days
    ),
    event = as.integer(!is.na(right_time))
  )

# Naive Cox
fit_cox_naive <- coxph(
  Surv(time_event, event) ~ plication_performed,
  data = df_cox, ties = "efron")
hr_naive <- exp(coef(fit_cox_naive)["plication_performed"])
ci_naive <- exp(confint(fit_cox_naive)["plication_performed", ])

# Time-dependent Cox
df_td <- df_cox %>%
  rowwise() %>%
  do({
    row <- .
    if (row$plication_performed == 1 &&
        !is.na(row$days_to_plication) &&
        row$days_to_plication < row$time_event) {
      bind_rows(
        data.frame(patient_id=row$patient_id, tstart=0,
                   tstop=row$days_to_plication, event_td=0, plication_td=0),
        data.frame(patient_id=row$patient_id, tstart=row$days_to_plication,
                   tstop=row$time_event, event_td=row$event, plication_td=1)
      )
    } else {
      data.frame(patient_id=row$patient_id, tstart=0,
                 tstop=row$time_event, event_td=row$event, plication_td=0)
    }
  }) %>% ungroup()

fit_cox_td <- coxph(
  Surv(tstart, tstop, event_td) ~ plication_td,
  data = df_td, ties = "efron")
hr_td <- exp(coef(fit_cox_td)["plication_td"])
ci_td <- exp(confint(fit_cox_td)["plication_td", ])
ph_p  <- cox.zph(fit_cox_td)$table["plication_td", "p"]

cat("Naive Cox HR:", round(hr_naive, 2),
    "(95% CI:", round(ci_naive[1], 2), "-", round(ci_naive[2], 2), ")\n")
cat("TD Cox HR   :", round(hr_td, 2),
    "(95% CI:", round(ci_td[1], 2), "-", round(ci_td[2], 2), ")\n")
cat("PH p        :", round(ph_p, 3), "\n\n")

# SA4 Forest plot
forest_df <- data.frame(
  Analysis = c("Naive Cox\n(immortal time bias present)",
               "Time-dependent Cox\n(immortal time bias corrected)"),
  HR    = c(hr_naive, hr_td),
  Lower = c(ci_naive[1], ci_td[1]),
  Upper = c(ci_naive[2], ci_td[2]),
  Bias  = c("Biased", "Corrected")
)

fig_sa4 <- ggplot(forest_df,
                  aes(x = HR, y = Analysis, color = Bias)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "grey50", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper),
                 height = 0.15, linewidth = 0.9) +
  geom_point(size = 4, shape = 18) +
  geom_text(aes(label = paste0(
  "HR ", round(HR, 2),
  "\n(95% CI ", round(Lower, 2), "-", round(Upper, 2), ")"
)), hjust = -0.12, size = 3.0, family = "sans") +  # serif → sans
  
theme_classic(base_family = "sans", base_size = 11) +
theme(
  legend.position    = "none",
  axis.text.y        = element_text(size = 9),
  panel.grid.major.x = element_line(color = "grey92"),
  plot.margin        = margin(10, 60, 10, 10)
)
# labs(title=, caption=) を削除
# x軸ラベルは残す
labs(
  x = "Hazard ratio (log scale)\n[HR > 1: faster recovery with plication]",
  y = NULL
)

ggsave(file.path(output_dir, "SuppFig_SA4_Cox_plication.pdf"),
       plot = fig_sa4, width = 7.0, height = 4.0, device = "pdf")
ggsave(file.path(output_dir, "SuppFig_SA4_Cox_plication.tiff"),
       plot = fig_sa4, width = 7.0, height = 4.0,
       dpi = 600, compression = "lzw")
cat("SA4 Figure complete\n\n")


# ============================================================
# 結果サマリー
# ============================================================
cat("========================================\n")
cat("SENSITIVITY ANALYSES SUMMARY\n")
cat("========================================\n")
cat("主解析  Median: 44 days (Bootstrap 95% CI: 41-75 days)\n")
cat("SA1     Median:", round(med_sa1$quantile), "days (両側麻痺除外)\n")
cat("SA2     Median:", round(med_sa2$quantile), "days (180日以内打ち切り除外)\n")
cat("SA3 Left  KM :", round(med_sa3_left$quantile),  "days\n")
cat("SA3 Mid   KM :", round(med_sa3_mid$quantile),   "days\n")
cat("SA3 Right KM :", round(med_sa3_right$quantile), "days\n")
cat("SA4 Naive HR :", round(hr_naive, 2),
    "(95% CI:", round(ci_naive[1], 2), "-", round(ci_naive[2], 2), ")\n")
cat("SA4 TD Cox HR:", round(hr_td, 2),
    "(95% CI:", round(ci_td[1], 2), "-", round(ci_td[2], 2), ")\n")
cat("SA4 PH p     :", round(ph_p, 3), "\n")
cat("========================================\n")
cat("\n07_sensitivity.R complete\n")

# ============================================================
# SA1: 両側麻痺例を全例除外
# ============================================================
cat("============================================================\n")
cat("SA1: 両側麻痺例を全例除外\n")
cat("============================================================\n")

df_sa1 <- df %>%
  filter(paralysis_side != "Bilateral")

cat(sprintf("除外症例数: %d（両側麻痺）\n", nrow(df) - nrow(df_sa1)))
cat(sprintf("解析対象: %d例\n", nrow(df_sa1)))
cat(sprintf("回復確認: %d例\n", sum(df_sa1$recovery_confirmed == 1, na.rm=TRUE)))
cat(sprintf("右側打ち切り: %d例\n\n", sum(df_sa1$recovery_confirmed == 0, na.rm=TRUE)))

fit_sa1 <- survfit(
  Surv(left_time, right_time, type = "interval2") ~ 1,
  data = df_sa1
)

med_sa1 <- quantile(fit_sa1, probs = 0.5)
cat(sprintf("SA1 Median recovery: %.0f days\n\n", med_sa1$quantile))


# ============================================================
# SA2: 180日以内右打ち切り除外
# ============================================================
cat("============================================================\n")
cat("SA2: 180日以内右打ち切り除外\n")
cat("============================================================\n")

df_sa2 <- df %>%
  filter(
    recovery_confirmed == 1 |
    (recovery_confirmed == 0 & followup_days > 180)
  )

cat(sprintf("除外症例数: %d（180日以内右打ち切り）\n", nrow(df) - nrow(df_sa2)))
cat(sprintf("解析対象: %d例\n", nrow(df_sa2)))
cat(sprintf("回復確認: %d例\n", sum(df_sa2$recovery_confirmed == 1, na.rm=TRUE)))
cat(sprintf("右側打ち切り: %d例\n\n", sum(df_sa2$recovery_confirmed == 0, na.rm=TRUE)))

fit_sa2 <- survfit(
  Surv(left_time, right_time, type = "interval2") ~ 1,
  data = df_sa2
)

med_sa2 <- quantile(fit_sa2, probs = 0.5)
cat(sprintf("SA2 Median recovery: %.0f days\n\n", med_sa2$quantile))


# ============================================================
# SA3: KM 3補完（left / mid / right endpoint）
# ============================================================
cat("============================================================\n")
cat("SA3: KM 3補完（left / mid / right endpoint）\n")
cat("============================================================\n")

# Left endpoint KM
df_sa3_left <- df %>%
  mutate(
    time_km = ifelse(!is.na(right_time), left_time, followup_days),
    event   = recovery_confirmed
  )
fit_sa3_left  <- survfit(Surv(time_km, event) ~ 1, data = df_sa3_left)
med_sa3_left  <- quantile(fit_sa3_left, probs = 0.5)

# Mid endpoint KM
df_sa3_mid <- df %>%
  mutate(
    time_km = ifelse(!is.na(right_time),
                     (left_time + right_time) / 2,
                     followup_days),
    event   = recovery_confirmed
  )
fit_sa3_mid  <- survfit(Surv(time_km, event) ~ 1, data = df_sa3_mid)
med_sa3_mid  <- quantile(fit_sa3_mid, probs = 0.5)

# Right endpoint KM
df_sa3_right <- df %>%
  mutate(
    time_km = ifelse(!is.na(right_time), right_time, followup_days),
    event   = recovery_confirmed
  )
fit_sa3_right <- survfit(Surv(time_km, event) ~ 1, data = df_sa3_right)
med_sa3_right <- quantile(fit_sa3_right, probs = 0.5)

cat("SA3 KM 3補完 結果:\n")
cat(sprintf("  Left endpoint  : %.0f days\n", med_sa3_left$quantile))
cat(sprintf("  Mid endpoint   : %.0f days\n", med_sa3_mid$quantile))
cat(sprintf("  Right endpoint : %.0f days\n\n", med_sa3_right$quantile))


# ============================================================
# SA4: Cox time-dependent plication（immortal time bias補正）
# ============================================================
cat("============================================================\n")
cat("SA4: Cox time-dependent plication（immortal time bias補正）\n")
cat("============================================================\n")

# midpoint imputation for Cox model
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
  data = df_cox,
  ties = "efron"
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
        data.frame(
          patient_id   = row$patient_id,
          tstart       = 0,
          tstop        = row$days_to_plication,
          event_td     = 0,
          plication_td = 0
        ),
        data.frame(
          patient_id   = row$patient_id,
          tstart       = row$days_to_plication,
          tstop        = row$time_event,
          event_td     = row$event,
          plication_td = 1
        )
      )
    } else {
      data.frame(
        patient_id   = row$patient_id,
        tstart       = 0,
        tstop        = row$time_event,
        event_td     = row$event,
        plication_td = 0
      )
    }
  }) %>%
  ungroup()

# Time-dependent Cox（immortal time bias補正）
fit_cox_td <- coxph(
  Surv(tstart, tstop, event_td) ~ plication_td,
  data = df_td,
  ties = "efron"
)

# PH仮定の検定
ph_test <- cox.zph(fit_cox_td)

hr_naive <- exp(coef(fit_cox_naive)["plication_performed"])
ci_naive <- exp(confint(fit_cox_naive)["plication_performed", ])
hr_td    <- exp(coef(fit_cox_td)["plication_td"])
ci_td    <- exp(confint(fit_cox_td)["plication_td", ])

cat("SA4 Cox sensitivity analysis results:\n")
cat(sprintf("  Naive Cox HR      : %.2f (95%% CI: %.2f\u2013%.2f)\n",
            hr_naive, ci_naive[1], ci_naive[2]))
cat(sprintf("  Time-dependent HR : %.2f (95%% CI: %.2f\u2013%.2f)\n",
            hr_td, ci_td[1], ci_td[2]))
cat(sprintf("  PH assumption p   : %.3f\n\n",
            ph_test$table["plication_td", "p"]))


# ============================================================
# 結果サマリー
# ============================================================
cat("========================================\n")
cat("SENSITIVITY ANALYSES SUMMARY\n")
cat("========================================\n")
cat("主解析  Median: 44 days (Bootstrap 95% CI: 41\u201375 days)\n")
cat(sprintf("SA1     Median: %.0f days（両側麻痺除外）\n",
            med_sa1$quantile))
cat(sprintf("SA2     Median: %.0f days（180日以内打ち切り除外）\n",
            med_sa2$quantile))
cat(sprintf("SA3 Left  KM : %.0f days\n", med_sa3_left$quantile))
cat(sprintf("SA3 Mid   KM : %.0f days\n", med_sa3_mid$quantile))
cat(sprintf("SA3 Right KM : %.0f days\n", med_sa3_right$quantile))
cat(sprintf("SA4 Cox TD HR: %.2f (95%% CI: %.2f\u2013%.2f)\n",
            hr_td, ci_td[1], ci_td[2]))
cat("========================================\n\n")


# ============================================================
# Figure: Supplement感度解析図（SA1・SA2の比較）
# ============================================================

# Turnbull曲線データ作成関数
make_tb_df <- function(fit, label, t_max = 1095) {
  data.frame(
    time = c(0, fit$time),
    surv = c(1, fit$surv)
  ) %>%
    mutate(
      cumrec = 1 - surv,
      group  = label
    ) %>%
    bind_rows(
      data.frame(
        time   = t_max,
        surv   = tail(.$surv, 1),
        cumrec = tail(.$cumrec, 1),
        group  = label
      )
    )
}

# 主解析のTurnbull曲線
fit_main <- survfit(
  Surv(left_time, right_time, type = "interval2") ~ 1,
  data = df
)

# ラベル定義
label_main <- "Primary analysis (n=110)"
label_sa1  <- "SA1: Exclude bilateral (n=101)"
label_sa2  <- sprintf("SA2: Exclude early censored (n=%d)", nrow(df_sa2))

# 曲線データ作成
tb_combined <- bind_rows(
  make_tb_df(fit_main, label_main),
  make_tb_df(fit_sa1,  label_sa1),
  make_tb_df(fit_sa2,  label_sa2)
)

# カラー設定
COL_MAIN <- "#1F3864"
COL_SA1  <- "#C0392B"
COL_SA2  <- "#2A9D8F"

fig_sa <- ggplot(tb_combined,
                 aes(x = time, y = cumrec * 100,
                     color = group, linetype = group,
                     linewidth = group)) +
  geom_step(direction = "hv") +
  scale_color_manual(
    values = setNames(
      c(COL_MAIN, COL_SA1, COL_SA2),
      c(label_main, label_sa1, label_sa2)
    )
  ) +
  scale_linetype_manual(
    values = setNames(
      c("solid", "dashed", "dotted"),
      c(label_main, label_sa1, label_sa2)
    )
  ) +
  scale_linewidth_manual(
    values = setNames(
      c(1.2, 1.0, 1.0),
      c(label_main, label_sa1, label_sa2)
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    name   = "Time from diaphragm paralysis diagnosis (days)",
    limits = c(0, 1095),
    breaks = seq(0, 1095, by = 182),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    name   = "Cumulative recovery probability (%)",
    limits = c(0, 100),
    breaks = seq(0, 100, by = 25),
    expand = expansion(mult = c(0, 0.02))
  ) +
  theme_classic(base_family = "serif", base_size = 11) +
  theme(
    legend.position      = c(0.98, 0.10),
    legend.justification = c(1, 0),
    legend.title         = element_blank(),
    legend.text          = element_text(size = 8.5, family = "serif"),
    legend.key.width     = unit(1.5, "cm"),
    legend.background    = element_rect(fill = "white", color = "grey80",
                                        linewidth = 0.3),
    panel.grid           = element_blank(),
    plot.title           = element_text(size = 11, face = "bold"),
    plot.caption         = element_text(size = 7.5, color = "grey50",
                                        hjust = 0),
    plot.margin          = margin(10, 15, 10, 10)
  ) +
  labs(
    title   = "Supplemental Figure. Sensitivity Analyses: Turnbull NPMLE",
    caption = paste0(
      "SA1: Patients with bilateral diaphragm paralysis excluded (n=9 excluded).\n",
      "SA2: Right-censored patients with follow-up \u2264180 days excluded.\n",
      "Primary analysis shown for reference."
    )
  )

print(fig_sa)


# ============================================================
# 出力（PDF + TIFF 600dpi）
# ============================================================
output_dir <- "~/Desktop"

ggsave(
  filename = file.path(output_dir, "SuppFig_Sensitivity_SA1_SA2.pdf"),
  plot     = fig_sa,
  width    = 7.0,
  height   = 5.5,
  device   = "pdf"
)

ggsave(
  filename    = file.path(output_dir, "SuppFig_Sensitivity_SA1_SA2.tiff"),
  plot        = fig_sa,
  width       = 7.0,
  height      = 5.5,
  dpi         = 600,
  compression = "lzw"
)

cat("\n\u2713 07_sensitivity.R 完了\n")
cat("保存先:", output_dir, "\n")
cat("  - SuppFig_Sensitivity_SA1_SA2.pdf\n")
cat("  - SuppFig_Sensitivity_SA1_SA2.tiff (600 dpi, LZW)\n\n")

# ============================================================
# SA1: 両側麻痺例を全例除外
# ============================================================
cat("============================================================\n")
cat("SA1: 両側麻痺例を全例除外\n")
cat("============================================================\n")

df_sa1 <- df %>%
  filter(paralysis_side != "Bilateral")

cat(sprintf("除外症例数: %d（両側麻痺）\n", nrow(df) - nrow(df_sa1)))
cat(sprintf("解析対象: %d例\n", nrow(df_sa1)))
cat(sprintf("回復確認: %d例\n", sum(df_sa1$recovery_confirmed == 1, na.rm=TRUE)))
cat(sprintf("右側打ち切り: %d例\n\n", sum(df_sa1$recovery_confirmed == 0, na.rm=TRUE)))

# Turnbull推定
fit_sa1 <- survfit(
  Surv(left_time, right_time, type = "interval2") ~ 1,
  data = df_sa1
)

med_sa1 <- quantile(fit_sa1, probs = 0.5)
cat(sprintf("SA1 Median recovery: %.0f days\n\n", med_sa1$quantile))


# ============================================================
# SA2: 180日以内右打ち切り除外
# ============================================================
cat("============================================================\n")
cat("SA2: 180日以内右打ち切り除外\n")
cat("============================================================\n")

df_sa2 <- df %>%
  filter(
    recovery_confirmed == 1 |                    # 回復確認例は全て含む
    (recovery_confirmed == 0 & followup_days > 180)  # 打ち切りは180日超のみ
  )

cat(sprintf("除外症例数: %d（180日以内右打ち切り）\n", nrow(df) - nrow(df_sa2)))
cat(sprintf("解析対象: %d例\n", nrow(df_sa2)))
cat(sprintf("回復確認: %d例\n", sum(df_sa2$recovery_confirmed == 1, na.rm=TRUE)))
cat(sprintf("右側打ち切り: %d例\n\n", sum(df_sa2$recovery_confirmed == 0, na.rm=TRUE)))

# Turnbull推定
fit_sa2 <- survfit(
  Surv(left_time, right_time, type = "interval2") ~ 1,
  data = df_sa2
)

med_sa2 <- quantile(fit_sa2, probs = 0.5)
cat(sprintf("SA2 Median recovery: %.0f days\n\n", med_sa2$quantile))


# ============================================================
# SA3: KM 3補完（left / mid / right endpoint）
# ============================================================
cat("============================================================\n")
cat("SA3: KM 3補完（left / mid / right endpoint）\n")
cat("============================================================\n")

# 回復確認例のみ（区間が存在する症例）
df_sa3 <- df %>% filter(!is.na(right_time))

cat(sprintf("対象症例数（回復確認例）: %d例\n\n", nrow(df_sa3)))

# Left endpoint KM
df_sa3_left <- df %>%
  mutate(
    time_km = ifelse(!is.na(right_time), left_time, followup_days),
    event   = recovery_confirmed
  )

fit_sa3_left <- survfit(Surv(time_km, event) ~ 1, data = df_sa3_left)
med_sa3_left <- quantile(fit_sa3_left, probs = 0.5)

# Mid endpoint KM
df_sa3_mid <- df %>%
  mutate(
    time_km = ifelse(!is.na(right_time),
                     (left_time + right_time) / 2,
                     followup_days),
    event   = recovery_confirmed
  )

fit_sa3_mid <- survfit(Surv(time_km, event) ~ 1, data = df_sa3_mid)
med_sa3_mid <- quantile(fit_sa3_mid, probs = 0.5)

# Right endpoint KM
df_sa3_right <- df %>%
  mutate(
    time_km = ifelse(!is.na(right_time), right_time, followup_days),
    event   = recovery_confirmed
  )

fit_sa3_right <- survfit(Surv(time_km, event) ~ 1, data = df_sa3_right)
med_sa3_right <- quantile(fit_sa3_right, probs = 0.5)

cat("SA3 KM 3補完 結果:\n")
cat(sprintf("  Left endpoint  : %.0f days\n", med_sa3_left$quantile))
cat(sprintf("  Mid endpoint   : %.0f days\n", med_sa3_mid$quantile))
cat(sprintf("  Right endpoint : %.0f days\n\n", med_sa3_right$quantile))


# ============================================================
# SA4: Cox time-dependent plication（immortal time bias補正）
# ============================================================
cat("============================================================\n")
cat("SA4: Cox time-dependent plication（immortal time bias補正）\n")
cat("============================================================\n")

# time-dependent covariate形式に変換
df_td <- df %>%
  mutate(
    time_event = case_when(
      !is.na(right_time) ~ (left_time + right_time) / 2,
      is.na(right_time)  ~ followup_days
    ),
    event = as.integer(!is.na(right_time))
  ) %>%
  rowwise() %>%
  do({
    row <- .
    if (row$plication_performed == 1 &&
        !is.na(row$days_to_plication) &&
        row$days_to_plication < row$time_event) {
      bind_rows(
        data.frame(
          patient_id   = row$patient_id,
          tstart       = 0,
          tstop        = row$days_to_plication,
          event_td     = 0,
          plication_td = 0
        ),
        data.frame(
          patient_id   = row$patient_id,
          tstart       = row$days_to_plication,
          tstop        = row$time_event,
          event_td     = row$event,
          plication_td = 1
        )
      )
    } else {
      data.frame(
        patient_id   = row$patient_id,
        tstart       = 0,
        tstop        = row$time_event,
        event_td     = row$event,
        plication_td = 0
      )
    }
  }) %>%
  ungroup()

# Naive Cox（immortal time biasあり・参照用）
fit_cox_naive <- coxph(
  Surv(time_event, event) ~ plication_performed,
  data = df %>%
    mutate(
      time_event = case_when(
        !is.na(right_time) ~ (left_time + right_time) / 2,
        is.na(right_time)  ~ followup_days
      ),
      event = as.integer(!is.na(right_time))
    ),
  ties = "efron"
)

# Time-dependent Cox（immortal time bias補正）
fit_cox_td <- coxph(
  Surv(tstart, tstop, event_td) ~ plication_td,
  data = df_td,
  ties = "efron"
)

# PH仮定の検定
ph_test <- cox.zph(fit_cox_td)

hr_naive <- exp(coef(fit_cox_naive)["plication_performed"])
ci_naive <- exp(confint(fit_cox_naive)["plication_performed", ])
hr_td    <- exp(coef(fit_cox_td)["plication_td"])
ci_td    <- exp(confint(fit_cox_td)["plication_td", ])

cat("SA4 Cox sensitivity analysis results:\n")
cat(sprintf("  Naive Cox HR      : %.2f (95%% CI: %.2f\u2013%.2f)\n",
            hr_naive, ci_naive[1], ci_naive[2]))
cat(sprintf("  Time-dependent HR : %.2f (95%% CI: %.2f\u2013%.2f)\n",
            hr_td, ci_td[1], ci_td[2]))
cat(sprintf("  PH assumption p   : %.3f\n\n",
            ph_test$table["plication_td", "p"]))


# ============================================================
# 結果サマリー
# ============================================================
cat("========================================\n")
cat("SENSITIVITY ANALYSES SUMMARY\n")
cat("========================================\n")
cat(sprintf("主解析  Median: 44 days (Bootstrap 95%% CI: 41\u201375 days)\n"))
cat(sprintf("SA1     Median: %.0f days（両側麻痺除外）\n",
            med_sa1$quantile))
cat(sprintf("SA2     Median: %.0f days（180日以内打ち切り除外）\n",
            med_sa2$quantile))
cat(sprintf("SA3 Left  KM : %.0f days\n", med_sa3_left$quantile))
cat(sprintf("SA3 Mid   KM : %.0f days\n", med_sa3_mid$quantile))
cat(sprintf("SA3 Right KM : %.0f days\n", med_sa3_right$quantile))
cat(sprintf("SA4 Cox TD HR: %.2f (95%% CI: %.2f\u2013%.2f)\n",
            hr_td, ci_td[1], ci_td[2]))
cat("========================================\n")


# ============================================================
# Figure: Supplement感度解析図（SA1・SA2の比較）
# ============================================================

# 主解析のTurnbull曲線
fit_main <- survfit(
  Surv(left_time, right_time, type = "interval2") ~ 1,
  data = df
)

make_tb_df <- function(fit, label, t_max = 1095) {
  df_out <- data.frame(
    time  = c(0, fit$time),
    surv  = c(1, fit$surv)
  ) %>%
    mutate(
      cumrec = 1 - surv,
      group  = label
    ) %>%
    bind_rows(
      data.frame(
        time   = t_max,
        surv   = tail(.$surv, 1),
        cumrec = tail(.$cumrec, 1),
        group  = label
      )
    )
  df_out
}

tb_main <- make_tb_df(fit_main, "Primary analysis (n=110)")
tb_sa1  <- make_tb_df(fit_sa1,  "SA1: Exclude bilateral (n=101)")
tb_sa2  <- make_tb_df(fit_sa2,  sprintf("SA2: Exclude early censored (n=%d)", nrow(df_sa2)))

tb_combined <- bind_rows(tb_main, tb_sa1, tb_sa2)

COL_MAIN <- "#1F3864"
COL_SA1  <- "#C0392B"
COL_SA2  <- "#2A9D8F"

label_main <- "Primary analysis (n=110)"
label_sa1  <- "SA1: Exclude bilateral (n=101)"
label_sa2  <- sprintf("SA2: Exclude early censored (n=%d)", nrow(df_sa2))

fig_sa <- ggplot(tb_combined,
                 aes(x = time, y = cumrec * 100,
                     color = group, linetype = group)) +
  geom_step(linewidth = 0.9, direction = "hv") +
  scale_color_manual(
    values = setNames(
      c(COL_MAIN, COL_SA1, COL_SA2),
      c(label_main, label_sa1, label_sa2)
    )
  ) +
  scale_linetype_manual(
    values = setNames(
      c("solid", "dashed", "dotted"),
      c(label_main, label_sa1, label_sa2)
    )
  ) +
  scale_x_continuous(
    name   = "Time from diaphragm paralysis diagnosis (days)",
    limits = c(0, 1095),
    breaks = seq(0, 1095, by = 182),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    name   = "Cumulative recovery probability (%)",
    limits = c(0, 100),
    breaks = seq(0, 100, by = 25),
    expand = expansion(mult = c(0, 0.02))
  ) +
  theme_classic(base_family = "sans", base_size = 11) +
  theme(
    legend.position      = c(0.98, 0.10),
    legend.justification = c(1, 0),
    legend.title         = element_blank(),
    legend.text          = element_text(size = 8.5, family = "sans"),
    legend.key.width     = unit(1.5, "cm"),
    legend.background    = element_rect(fill = "white", color = "grey80",
                                        linewidth = 0.3),
    panel.grid           = element_blank(),
    plot.margin          = margin(10, 15, 10, 10)
  )

print(fig_sa)

# 出力
output_dir <- "~/Desktop"

ggsave(
  filename = file.path(output_dir, "SuppFig_Sensitivity_SA1_SA2.pdf"),
  plot     = fig_sa,
  width    = 7.0,
  height   = 5.5,
  device   = "pdf"
)

ggsave(
  filename    = file.path(output_dir, "SuppFig_Sensitivity_SA1_SA2.tiff"),
  plot        = fig_sa,
  width       = 7.0,
  height      = 5.5,
  dpi         = 600,
  compression = "lzw"
)

cat("\n\u2713 07_sensitivity.R 完了\n")
cat("保存先:", output_dir, "\n")
cat("  - SuppFig_Sensitivity_SA1_SA2.pdf\n")
cat("  - SuppFig_Sensitivity_SA1_SA2.tiff\n\n")
# ============================================================
# 共通関数
# ============================================================
make_tb_df <- function(fit, label, t_max = 1095) {
  data.frame(time = c(0, fit$time), surv = c(1, fit$surv)) %>%
    mutate(cumrec = 1 - surv, group = label) %>%
    bind_rows(data.frame(
      time = t_max, surv = tail(.$surv, 1),
      cumrec = tail(.$cumrec, 1), group = label
    ))
}


# ============================================================
# SA1: 両側麻痺例を全例除外
# ============================================================
cat("SA1: 両側麻痺例除外\n")
df_sa1 <- df %>% filter(paralysis_side != "Bilateral")
fit_sa1 <- survfit(
  Surv(left_time, right_time, type = "interval2") ~ 1, data = df_sa1)
med_sa1 <- quantile(fit_sa1, probs = 0.5)
cat(sprintf("n=%d, Median=%.0f days\n\n", nrow(df_sa1), med_sa1$quantile))


# ============================================================
# SA2: 180日以内右打ち切り除外
# ============================================================
cat("SA2: 180日以内右打ち切り除外\n")
df_sa2 <- df %>%
  filter(recovery_confirmed == 1 |
         (recovery_confirmed == 0 & followup_days > 180))
fit_sa2 <- survfit(
  Surv(left_time, right_time, type = "interval2") ~ 1, data = df_sa2)
med_sa2 <- quantile(fit_sa2, probs = 0.5)
cat(sprintf("n=%d, Median=%.0f days\n\n", nrow(df_sa2), med_sa2$quantile))


# ============================================================
# SA1・SA2 Figure
# ============================================================
fit_main <- survfit(
  Surv(left_time, right_time, type = "interval2") ~ 1, data = df)

label_main <- "Primary analysis (n=110)"
label_sa1  <- "SA1: Exclude bilateral (n=101)"
label_sa2  <- paste0("SA2: Exclude early censored (n=", nrow(df_sa2), ")")

tb_sa12 <- bind_rows(
  make_tb_df(fit_main, label_main),
  make_tb_df(fit_sa1,  label_sa1),
  make_tb_df(fit_sa2,  label_sa2)
)

COL_MAIN <- "#1F3864"
COL_SA1  <- "#C0392B"
COL_SA2  <- "#2A9D8F"

fig_sa12 <- ggplot(tb_sa12,
                   aes(x = time, y = cumrec * 100,
                       color = group, linetype = group,
                       linewidth = group)) +
  geom_step(direction = "hv") +
  scale_color_manual(values = setNames(
    c(COL_MAIN, COL_SA1, COL_SA2),
    c(label_main, label_sa1, label_sa2)
  )) +
  scale_linetype_manual(values = setNames(
    c("solid", "dashed", "dotted"),
    c(label_main, label_sa1, label_sa2)
  )) +
  scale_linewidth_manual(values = setNames(
    c(1.2, 1.0, 1.0),
    c(label_main, label_sa1, label_sa2)
  ), guide = "none") +
  scale_x_continuous(
    name = "Time from diaphragm paralysis diagnosis (days)",
    limits = c(0, 1095), breaks = seq(0, 1095, by = 182),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    name = "Cumulative recovery probability (%)",
    limits = c(0, 100), breaks = seq(0, 100, by = 25),
    expand = expansion(mult = c(0, 0.02))
  ) +
  theme_classic(base_family = "serif", base_size = 11) +
  theme(
    legend.position      = c(0.98, 0.10),
    legend.justification = c(1, 0),
    legend.title         = element_blank(),
    legend.text          = element_text(size = 8.5, family = "serif"),
    legend.key.width     = unit(1.5, "cm"),
    legend.background    = element_rect(fill = "white", color = "grey80",
                                        linewidth = 0.3),
    panel.grid  = element_blank(),
    plot.title  = element_text(size = 11, face = "bold"),
    plot.caption = element_text(size = 7.5, color = "grey50", hjust = 0),
    plot.margin = margin(10, 15, 10, 10)
  ) +
  labs(
    title = "Supplemental Figure. Sensitivity Analyses: Turnbull NPMLE",
    caption = paste0(
      "SA1: Patients with bilateral diaphragm paralysis excluded (n=9 excluded).\n",
      "SA2: Right-censored patients with follow-up \u2264180 days excluded.\n",
      "Primary analysis shown for reference."
    )
  )

ggsave(file.path(output_dir, "SuppFig_Sensitivity_SA1_SA2.pdf"),
       plot = fig_sa12, width = 7.0, height = 5.5, device = "pdf")
ggsave(file.path(output_dir, "SuppFig_Sensitivity_SA1_SA2.tiff"),
       plot = fig_sa12, width = 7.0, height = 5.5,
       dpi = 600, compression = "lzw")
cat("SA1/SA2 Figure complete\n\n")


# ============================================================
# SA3: KM 3補完（left / mid / right endpoint）
# ============================================================
cat("SA3: KM 3補完\n")

df_sa3_left <- df %>%
  mutate(time_km = ifelse(!is.na(right_time), left_time, followup_days),
         event = recovery_confirmed)
fit_sa3_left  <- survfit(Surv(time_km, event) ~ 1, data = df_sa3_left)
med_sa3_left  <- quantile(fit_sa3_left, probs = 0.5)

df_sa3_mid <- df %>%
  mutate(time_km = ifelse(!is.na(right_time),
                          (left_time + right_time) / 2, followup_days),
         event = recovery_confirmed)
fit_sa3_mid  <- survfit(Surv(time_km, event) ~ 1, data = df_sa3_mid)
med_sa3_mid  <- quantile(fit_sa3_mid, probs = 0.5)

df_sa3_right <- df %>%
  mutate(time_km = ifelse(!is.na(right_time), right_time, followup_days),
         event = recovery_confirmed)
fit_sa3_right <- survfit(Surv(time_km, event) ~ 1, data = df_sa3_right)
med_sa3_right <- quantile(fit_sa3_right, probs = 0.5)

cat("Left  endpoint:", round(med_sa3_left$quantile),  "days\n")
cat("Mid   endpoint:", round(med_sa3_mid$quantile),   "days\n")
cat("Right endpoint:", round(med_sa3_right$quantile), "days\n\n")

# SA3 Figure
label_turnbull <- "Turnbull NPMLE (median=44d)"
label_left     <- "KM Left endpoint (median=16d)"
label_mid      <- "KM Mid endpoint (median=109d)"
label_right    <- "KM Right endpoint (median=205d)"

tb_sa3 <- bind_rows(
  make_tb_df(fit_main,      label_turnbull),
  make_tb_df(fit_sa3_left,  label_left),
  make_tb_df(fit_sa3_mid,   label_mid),
  make_tb_df(fit_sa3_right, label_right)
)

COL_TURNBULL <- "#1F3864"
COL_LEFT     <- "#C0392B"
COL_MID      <- "#E07B39"
COL_RIGHT    <- "#2A9D8F"

fig_sa3 <- ggplot(tb_sa3,
                  aes(x = time, y = cumrec * 100,
                      color = group, linetype = group,
                      linewidth = group)) +
  geom_step(direction = "hv") +
  geom_hline(yintercept = 50, linetype = "longdash",
             color = "grey40", linewidth = 0.5) +
  scale_color_manual(values = setNames(
    c(COL_TURNBULL, COL_LEFT, COL_MID, COL_RIGHT),
    c(label_turnbull, label_left, label_mid, label_right)
  )) +
  scale_linetype_manual(values = setNames(
    c("solid", "dashed", "dotted", "dotdash"),
    c(label_turnbull, label_left, label_mid, label_right)
  )) +
  scale_linewidth_manual(values = setNames(
    c(1.2, 0.9, 0.9, 0.9),
    c(label_turnbull, label_left, label_mid, label_right)
  ), guide = "none") +
  scale_x_continuous(
    name = "Time from diaphragm paralysis diagnosis (days)",
    limits = c(0, 1095), breaks = seq(0, 1095, by = 182),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    name = "Cumulative recovery probability (%)",
    limits = c(0, 100), breaks = seq(0, 100, by = 25),
    expand = expansion(mult = c(0, 0.02))
  ) +
  theme_classic(base_family = "serif", base_size = 11) +
  theme(
    legend.position      = c(0.98, 0.10),
    legend.justification = c(1, 0),
    legend.title         = element_blank(),
    legend.text          = element_text(size = 8.0, family = "serif"),
    legend.key.width     = unit(1.5, "cm"),
    legend.background    = element_rect(fill = "white", color = "grey80",
                                        linewidth = 0.3),
    panel.grid  = element_blank(),
    plot.title  = element_text(size = 11, face = "bold"),
    plot.caption = element_text(size = 7.5, color = "grey50", hjust = 0),
    plot.margin = margin(10, 15, 10, 10)
  ) +
  labs(
    title = "Supplemental Figure. SA3: KM Three-point Imputation",
    caption = paste0(
      "Dashed horizontal line indicates 50% cumulative recovery probability.\n",
      "Left endpoint: recovery assumed at last abnormal imaging date.\n",
      "Mid endpoint: recovery assumed at midpoint of censoring interval.\n",
      "Right endpoint: recovery assumed at first normal imaging date.\n",
      "Turnbull NPMLE (primary analysis) shown for reference."
    )
  )

ggsave(file.path(output_dir, "SuppFig_SA3_KM_imputation.pdf"),
       plot = fig_sa3, width = 7.0, height = 5.5, device = "pdf")
ggsave(file.path(output_dir, "SuppFig_SA3_KM_imputation.tiff"),
       plot = fig_sa3, width = 7.0, height = 5.5,
       dpi = 600, compression = "lzw")
cat("SA3 Figure complete\n\n")


# ============================================================
# SA4: Cox time-dependent plication（immortal time bias補正）
# ============================================================
cat("SA4: Cox time-dependent plication\n")

df_cox <- df %>%
  mutate(
    time_event = case_when(
      !is.na(right_time) ~ (left_time + right_time) / 2,
      is.na(right_time)  ~ followup_days
    ),
    event = as.integer(!is.na(right_time))
  )

# Naive Cox
fit_cox_naive <- coxph(
  Surv(time_event, event) ~ plication_performed,
  data = df_cox, ties = "efron")
hr_naive <- exp(coef(fit_cox_naive)["plication_performed"])
ci_naive <- exp(confint(fit_cox_naive)["plication_performed", ])

# Time-dependent Cox
df_td <- df_cox %>%
  rowwise() %>%
  do({
    row <- .
    if (row$plication_performed == 1 &&
        !is.na(row$days_to_plication) &&
        row$days_to_plication < row$time_event) {
      bind_rows(
        data.frame(patient_id=row$patient_id, tstart=0,
                   tstop=row$days_to_plication, event_td=0, plication_td=0),
        data.frame(patient_id=row$patient_id, tstart=row$days_to_plication,
                   tstop=row$time_event, event_td=row$event, plication_td=1)
      )
    } else {
      data.frame(patient_id=row$patient_id, tstart=0,
                 tstop=row$time_event, event_td=row$event, plication_td=0)
    }
  }) %>% ungroup()

fit_cox_td <- coxph(
  Surv(tstart, tstop, event_td) ~ plication_td,
  data = df_td, ties = "efron")
hr_td <- exp(coef(fit_cox_td)["plication_td"])
ci_td <- exp(confint(fit_cox_td)["plication_td", ])
ph_p  <- cox.zph(fit_cox_td)$table["plication_td", "p"]

cat("Naive Cox HR:", round(hr_naive, 2),
    "(95% CI:", round(ci_naive[1], 2), "-", round(ci_naive[2], 2), ")\n")
cat("TD Cox HR   :", round(hr_td, 2),
    "(95% CI:", round(ci_td[1], 2), "-", round(ci_td[2], 2), ")\n")
cat("PH p        :", round(ph_p, 3), "\n\n")

# SA4 Forest plot
forest_df <- data.frame(
  Analysis = c("Naive Cox\n(immortal time bias present)",
               "Time-dependent Cox\n(immortal time bias corrected)"),
  HR    = c(hr_naive, hr_td),
  Lower = c(ci_naive[1], ci_td[1]),
  Upper = c(ci_naive[2], ci_td[2]),
  Bias  = c("Biased", "Corrected")
)

fig_sa4 <- ggplot(forest_df,
                  aes(x = HR, y = Analysis, color = Bias)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "grey50", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper),
                 height = 0.15, linewidth = 0.9) +
  geom_point(size = 4, shape = 18) +
  geom_text(aes(label = paste0(
    "HR ", round(HR, 2),
    "\n(95% CI ", round(Lower, 2), "-", round(Upper, 2), ")"
  )), hjust = -0.12, size = 3.0, family = "serif") +
  scale_color_manual(values = c("Biased"    = "#C0392B",
                                "Corrected" = "#1F3864")) +
  scale_x_log10(limits = c(0.2, 3.0),
                breaks  = c(0.25, 0.5, 1, 2),
                labels  = c("0.25", "0.5", "1", "2")) +
  theme_classic(base_family = "serif", base_size = 11) +
  theme(
    legend.position    = "none",
    axis.text.y        = element_text(size = 9),
    panel.grid.major.x = element_line(color = "grey92"),
    plot.title   = element_text(size = 11, face = "bold"),
    plot.caption = element_text(size = 7.5, color = "grey50", hjust = 0),
    plot.margin  = margin(10, 60, 10, 10)
  ) +
  labs(
    title = "Supplemental Figure. SA4: Effect of Diaphragmatic Plication on Recovery",
    x = "Hazard ratio (log scale)\n[HR > 1: faster recovery with plication]",
    y = NULL,
    caption = paste0(
      "Red: naive Cox model with immortal time bias.\n",
      "Blue: time-dependent Cox model correcting for immortal time bias.\n",
      sprintf("PH assumption: p = %.3f (proportional hazards assumption holds).", ph_p)
    )
  )

ggsave(file.path(output_dir, "SuppFig_SA4_Cox_plication.pdf"),
       plot = fig_sa4, width = 7.0, height = 4.0, device = "pdf")
ggsave(file.path(output_dir, "SuppFig_SA4_Cox_plication.tiff"),
       plot = fig_sa4, width = 7.0, height = 4.0,
       dpi = 600, compression = "lzw")
cat("SA4 Figure complete\n\n")


# ============================================================
# 結果サマリー
# ============================================================
cat("========================================\n")
cat("SENSITIVITY ANALYSES SUMMARY\n")
cat("========================================\n")
cat("主解析  Median: 44 days (Bootstrap 95% CI: 41-75 days)\n")
cat("SA1     Median:", round(med_sa1$quantile), "days (両側麻痺除外)\n")
cat("SA2     Median:", round(med_sa2$quantile), "days (180日以内打ち切り除外)\n")
cat("SA3 Left  KM :", round(med_sa3_left$quantile),  "days\n")
cat("SA3 Mid   KM :", round(med_sa3_mid$quantile),   "days\n")
cat("SA3 Right KM :", round(med_sa3_right$quantile), "days\n")
cat("SA4 Naive HR :", round(hr_naive, 2),
    "(95% CI:", round(ci_naive[1], 2), "-", round(ci_naive[2], 2), ")\n")
cat("SA4 TD Cox HR:", round(hr_td, 2),
    "(95% CI:", round(ci_td[1], 2), "-", round(ci_td[2], 2), ")\n")
cat("SA4 PH p     :", round(ph_p, 3), "\n")
cat("========================================\n")
cat("\n07_sensitivity.R complete\n")

# ============================================================
# SA1: 両側麻痺例を全例除外
# ============================================================
cat("============================================================\n")
cat("SA1: 両側麻痺例を全例除外\n")
cat("============================================================\n")

df_sa1 <- df %>%
  filter(paralysis_side != "Bilateral")

cat(sprintf("除外症例数: %d（両側麻痺）\n", nrow(df) - nrow(df_sa1)))
cat(sprintf("解析対象: %d例\n", nrow(df_sa1)))
cat(sprintf("回復確認: %d例\n", sum(df_sa1$recovery_confirmed == 1, na.rm=TRUE)))
cat(sprintf("右側打ち切り: %d例\n\n", sum(df_sa1$recovery_confirmed == 0, na.rm=TRUE)))

fit_sa1 <- survfit(
  Surv(left_time, right_time, type = "interval2") ~ 1,
  data = df_sa1
)

med_sa1 <- quantile(fit_sa1, probs = 0.5)
cat(sprintf("SA1 Median recovery: %.0f days\n\n", med_sa1$quantile))


# ============================================================
# SA2: 180日以内右打ち切り除外
# ============================================================
cat("============================================================\n")
cat("SA2: 180日以内右打ち切り除外\n")
cat("============================================================\n")

df_sa2 <- df %>%
  filter(
    recovery_confirmed == 1 |
    (recovery_confirmed == 0 & followup_days > 180)
  )

cat(sprintf("除外症例数: %d（180日以内右打ち切り）\n", nrow(df) - nrow(df_sa2)))
cat(sprintf("解析対象: %d例\n", nrow(df_sa2)))
cat(sprintf("回復確認: %d例\n", sum(df_sa2$recovery_confirmed == 1, na.rm=TRUE)))
cat(sprintf("右側打ち切り: %d例\n\n", sum(df_sa2$recovery_confirmed == 0, na.rm=TRUE)))

fit_sa2 <- survfit(
  Surv(left_time, right_time, type = "interval2") ~ 1,
  data = df_sa2
)

med_sa2 <- quantile(fit_sa2, probs = 0.5)
cat(sprintf("SA2 Median recovery: %.0f days\n\n", med_sa2$quantile))


# ============================================================
# SA3: KM 3補完（left / mid / right endpoint）
# ============================================================
cat("============================================================\n")
cat("SA3: KM 3補完（left / mid / right endpoint）\n")
cat("============================================================\n")

# Left endpoint KM
df_sa3_left <- df %>%
  mutate(
    time_km = ifelse(!is.na(right_time), left_time, followup_days),
    event   = recovery_confirmed
  )
fit_sa3_left  <- survfit(Surv(time_km, event) ~ 1, data = df_sa3_left)
med_sa3_left  <- quantile(fit_sa3_left, probs = 0.5)

# Mid endpoint KM
df_sa3_mid <- df %>%
  mutate(
    time_km = ifelse(!is.na(right_time),
                     (left_time + right_time) / 2,
                     followup_days),
    event   = recovery_confirmed
  )
fit_sa3_mid  <- survfit(Surv(time_km, event) ~ 1, data = df_sa3_mid)
med_sa3_mid  <- quantile(fit_sa3_mid, probs = 0.5)

# Right endpoint KM
df_sa3_right <- df %>%
  mutate(
    time_km = ifelse(!is.na(right_time), right_time, followup_days),
    event   = recovery_confirmed
  )
fit_sa3_right <- survfit(Surv(time_km, event) ~ 1, data = df_sa3_right)
med_sa3_right <- quantile(fit_sa3_right, probs = 0.5)

cat("SA3 KM 3補完 結果:\n")
cat(sprintf("  Left endpoint  : %.0f days\n", med_sa3_left$quantile))
cat(sprintf("  Mid endpoint   : %.0f days\n", med_sa3_mid$quantile))
cat(sprintf("  Right endpoint : %.0f days\n\n", med_sa3_right$quantile))


# ============================================================
# SA4: Cox time-dependent plication（immortal time bias補正）
# ============================================================
cat("============================================================\n")
cat("SA4: Cox time-dependent plication（immortal time bias補正）\n")
cat("============================================================\n")

# midpoint imputation for Cox model
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
  data = df_cox,
  ties = "efron"
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
        data.frame(
          patient_id   = row$patient_id,
          tstart       = 0,
          tstop        = row$days_to_plication,
          event_td     = 0,
          plication_td = 0
        ),
        data.frame(
          patient_id   = row$patient_id,
          tstart       = row$days_to_plication,
          tstop        = row$time_event,
          event_td     = row$event,
          plication_td = 1
        )
      )
    } else {
      data.frame(
        patient_id   = row$patient_id,
        tstart       = 0,
        tstop        = row$time_event,
        event_td     = row$event,
        plication_td = 0
      )
    }
  }) %>%
  ungroup()

# Time-dependent Cox（immortal time bias補正）
fit_cox_td <- coxph(
  Surv(tstart, tstop, event_td) ~ plication_td,
  data = df_td,
  ties = "efron"
)

# PH仮定の検定
ph_test <- cox.zph(fit_cox_td)

hr_naive <- exp(coef(fit_cox_naive)["plication_performed"])
ci_naive <- exp(confint(fit_cox_naive)["plication_performed", ])
hr_td    <- exp(coef(fit_cox_td)["plication_td"])
ci_td    <- exp(confint(fit_cox_td)["plication_td", ])

cat("SA4 Cox sensitivity analysis results:\n")
cat(sprintf("  Naive Cox HR      : %.2f (95%% CI: %.2f\u2013%.2f)\n",
            hr_naive, ci_naive[1], ci_naive[2]))
cat(sprintf("  Time-dependent HR : %.2f (95%% CI: %.2f\u2013%.2f)\n",
            hr_td, ci_td[1], ci_td[2]))
cat(sprintf("  PH assumption p   : %.3f\n\n",
            ph_test$table["plication_td", "p"]))


# ============================================================
# 結果サマリー
# ============================================================
cat("========================================\n")
cat("SENSITIVITY ANALYSES SUMMARY\n")
cat("========================================\n")
cat("主解析  Median: 44 days (Bootstrap 95% CI: 41\u201375 days)\n")
cat(sprintf("SA1     Median: %.0f days（両側麻痺除外）\n",
            med_sa1$quantile))
cat(sprintf("SA2     Median: %.0f days（180日以内打ち切り除外）\n",
            med_sa2$quantile))
cat(sprintf("SA3 Left  KM : %.0f days\n", med_sa3_left$quantile))
cat(sprintf("SA3 Mid   KM : %.0f days\n", med_sa3_mid$quantile))
cat(sprintf("SA3 Right KM : %.0f days\n", med_sa3_right$quantile))
cat(sprintf("SA4 Cox TD HR: %.2f (95%% CI: %.2f\u2013%.2f)\n",
            hr_td, ci_td[1], ci_td[2]))
cat("========================================\n\n")


# ============================================================
# Figure: Supplement感度解析図（SA1・SA2の比較）
# ============================================================

# Turnbull曲線データ作成関数
make_tb_df <- function(fit, label, t_max = 1095) {
  data.frame(
    time = c(0, fit$time),
    surv = c(1, fit$surv)
  ) %>%
    mutate(
      cumrec = 1 - surv,
      group  = label
    ) %>%
    bind_rows(
      data.frame(
        time   = t_max,
        surv   = tail(.$surv, 1),
        cumrec = tail(.$cumrec, 1),
        group  = label
      )
    )
}

# 主解析のTurnbull曲線
fit_main <- survfit(
  Surv(left_time, right_time, type = "interval2") ~ 1,
  data = df
)

# ラベル定義
label_main <- "Primary analysis (n=110)"
label_sa1  <- "SA1: Exclude bilateral (n=101)"
label_sa2  <- sprintf("SA2: Exclude early censored (n=%d)", nrow(df_sa2))

# 曲線データ作成
tb_combined <- bind_rows(
  make_tb_df(fit_main, label_main),
  make_tb_df(fit_sa1,  label_sa1),
  make_tb_df(fit_sa2,  label_sa2)
)

# カラー設定
COL_MAIN <- "#1F3864"
COL_SA1  <- "#C0392B"
COL_SA2  <- "#2A9D8F"

fig_sa <- ggplot(tb_combined,
                 aes(x = time, y = cumrec * 100,
                     color = group, linetype = group,
                     linewidth = group)) +
  geom_step(direction = "hv") +
  scale_color_manual(
    values = setNames(
      c(COL_MAIN, COL_SA1, COL_SA2),
      c(label_main, label_sa1, label_sa2)
    )
  ) +
  scale_linetype_manual(
    values = setNames(
      c("solid", "dashed", "dotted"),
      c(label_main, label_sa1, label_sa2)
    )
  ) +
  scale_linewidth_manual(
    values = setNames(
      c(1.2, 1.0, 1.0),
      c(label_main, label_sa1, label_sa2)
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    name   = "Time from diaphragm paralysis diagnosis (days)",
    limits = c(0, 1095),
    breaks = seq(0, 1095, by = 182),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    name   = "Cumulative recovery probability (%)",
    limits = c(0, 100),
    breaks = seq(0, 100, by = 25),
    expand = expansion(mult = c(0, 0.02))
  ) +
  theme_classic(base_family = "serif", base_size = 11) +
  theme(
    legend.position      = c(0.98, 0.10),
    legend.justification = c(1, 0),
    legend.title         = element_blank(),
    legend.text          = element_text(size = 8.5, family = "serif"),
    legend.key.width     = unit(1.5, "cm"),
    legend.background    = element_rect(fill = "white", color = "grey80",
                                        linewidth = 0.3),
    panel.grid           = element_blank(),
    plot.title           = element_text(size = 11, face = "bold"),
    plot.caption         = element_text(size = 7.5, color = "grey50",
                                        hjust = 0),
    plot.margin          = margin(10, 15, 10, 10)
  ) +
  labs(
    title   = "Supplemental Figure. Sensitivity Analyses: Turnbull NPMLE",
    caption = paste0(
      "SA1: Patients with bilateral diaphragm paralysis excluded (n=9 excluded).\n",
      "SA2: Right-censored patients with follow-up \u2264180 days excluded.\n",
      "Primary analysis shown for reference."
    )
  )

print(fig_sa)


# ============================================================
# 出力（PDF + TIFF 600dpi）
# ============================================================
output_dir <- "~/Desktop"

ggsave(
  filename = file.path(output_dir, "SuppFig_Sensitivity_SA1_SA2.pdf"),
  plot     = fig_sa,
  width    = 7.0,
  height   = 5.5,
  device   = "pdf"
)

ggsave(
  filename    = file.path(output_dir, "SuppFig_Sensitivity_SA1_SA2.tiff"),
  plot        = fig_sa,
  width       = 7.0,
  height      = 5.5,
  dpi         = 600,
  compression = "lzw"
)

cat("\n\u2713 07_sensitivity.R 完了\n")
cat("保存先:", output_dir, "\n")
cat("  - SuppFig_Sensitivity_SA1_SA2.pdf\n")
cat("  - SuppFig_Sensitivity_SA1_SA2.tiff (600 dpi, LZW)\n\n")

# ============================================================
# SA1: 両側麻痺例を全例除外
# ============================================================
cat("============================================================\n")
cat("SA1: 両側麻痺例を全例除外\n")
cat("============================================================\n")

df_sa1 <- df %>%
  filter(paralysis_side != "Bilateral")

cat(sprintf("除外症例数: %d（両側麻痺）\n", nrow(df) - nrow(df_sa1)))
cat(sprintf("解析対象: %d例\n", nrow(df_sa1)))
cat(sprintf("回復確認: %d例\n", sum(df_sa1$recovery_confirmed == 1, na.rm=TRUE)))
cat(sprintf("右側打ち切り: %d例\n\n", sum(df_sa1$recovery_confirmed == 0, na.rm=TRUE)))

# Turnbull推定
fit_sa1 <- survfit(
  Surv(left_time, right_time, type = "interval2") ~ 1,
  data = df_sa1
)

med_sa1 <- quantile(fit_sa1, probs = 0.5)
cat(sprintf("SA1 Median recovery: %.0f days\n\n", med_sa1$quantile))


# ============================================================
# SA2: 180日以内右打ち切り除外
# ============================================================
cat("============================================================\n")
cat("SA2: 180日以内右打ち切り除外\n")
cat("============================================================\n")

df_sa2 <- df %>%
  filter(
    recovery_confirmed == 1 |                    # 回復確認例は全て含む
    (recovery_confirmed == 0 & followup_days > 180)  # 打ち切りは180日超のみ
  )

cat(sprintf("除外症例数: %d（180日以内右打ち切り）\n", nrow(df) - nrow(df_sa2)))
cat(sprintf("解析対象: %d例\n", nrow(df_sa2)))
cat(sprintf("回復確認: %d例\n", sum(df_sa2$recovery_confirmed == 1, na.rm=TRUE)))
cat(sprintf("右側打ち切り: %d例\n\n", sum(df_sa2$recovery_confirmed == 0, na.rm=TRUE)))

# Turnbull推定
fit_sa2 <- survfit(
  Surv(left_time, right_time, type = "interval2") ~ 1,
  data = df_sa2
)

med_sa2 <- quantile(fit_sa2, probs = 0.5)
cat(sprintf("SA2 Median recovery: %.0f days\n\n", med_sa2$quantile))


# ============================================================
# SA3: KM 3補完（left / mid / right endpoint）
# ============================================================
cat("============================================================\n")
cat("SA3: KM 3補完（left / mid / right endpoint）\n")
cat("============================================================\n")

# 回復確認例のみ（区間が存在する症例）
df_sa3 <- df %>% filter(!is.na(right_time))

cat(sprintf("対象症例数（回復確認例）: %d例\n\n", nrow(df_sa3)))

# Left endpoint KM
df_sa3_left <- df %>%
  mutate(
    time_km = ifelse(!is.na(right_time), left_time, followup_days),
    event   = recovery_confirmed
  )

fit_sa3_left <- survfit(Surv(time_km, event) ~ 1, data = df_sa3_left)
med_sa3_left <- quantile(fit_sa3_left, probs = 0.5)

# Mid endpoint KM
df_sa3_mid <- df %>%
  mutate(
    time_km = ifelse(!is.na(right_time),
                     (left_time + right_time) / 2,
                     followup_days),
    event   = recovery_confirmed
  )

fit_sa3_mid <- survfit(Surv(time_km, event) ~ 1, data = df_sa3_mid)
med_sa3_mid <- quantile(fit_sa3_mid, probs = 0.5)

# Right endpoint KM
df_sa3_right <- df %>%
  mutate(
    time_km = ifelse(!is.na(right_time), right_time, followup_days),
    event   = recovery_confirmed
  )

fit_sa3_right <- survfit(Surv(time_km, event) ~ 1, data = df_sa3_right)
med_sa3_right <- quantile(fit_sa3_right, probs = 0.5)

cat("SA3 KM 3補完 結果:\n")
cat(sprintf("  Left endpoint  : %.0f days\n", med_sa3_left$quantile))
cat(sprintf("  Mid endpoint   : %.0f days\n", med_sa3_mid$quantile))
cat(sprintf("  Right endpoint : %.0f days\n\n", med_sa3_right$quantile))


# ============================================================
# SA4: Cox time-dependent plication（immortal time bias補正）
# ============================================================
cat("============================================================\n")
cat("SA4: Cox time-dependent plication（immortal time bias補正）\n")
cat("============================================================\n")

# time-dependent covariate形式に変換
df_td <- df %>%
  mutate(
    time_event = case_when(
      !is.na(right_time) ~ (left_time + right_time) / 2,
      is.na(right_time)  ~ followup_days
    ),
    event = as.integer(!is.na(right_time))
  ) %>%
  rowwise() %>%
  do({
    row <- .
    if (row$plication_performed == 1 &&
        !is.na(row$days_to_plication) &&
        row$days_to_plication < row$time_event) {
      bind_rows(
        data.frame(
          patient_id   = row$patient_id,
          tstart       = 0,
          tstop        = row$days_to_plication,
          event_td     = 0,
          plication_td = 0
        ),
        data.frame(
          patient_id   = row$patient_id,
          tstart       = row$days_to_plication,
          tstop        = row$time_event,
          event_td     = row$event,
          plication_td = 1
        )
      )
    } else {
      data.frame(
        patient_id   = row$patient_id,
        tstart       = 0,
        tstop        = row$time_event,
        event_td     = row$event,
        plication_td = 0
      )
    }
  }) %>%
  ungroup()

# Naive Cox（immortal time biasあり・参照用）
fit_cox_naive <- coxph(
  Surv(time_event, event) ~ plication_performed,
  data = df %>%
    mutate(
      time_event = case_when(
        !is.na(right_time) ~ (left_time + right_time) / 2,
        is.na(right_time)  ~ followup_days
      ),
      event = as.integer(!is.na(right_time))
    ),
  ties = "efron"
)

# Time-dependent Cox（immortal time bias補正）
fit_cox_td <- coxph(
  Surv(tstart, tstop, event_td) ~ plication_td,
  data = df_td,
  ties = "efron"
)

# PH仮定の検定
ph_test <- cox.zph(fit_cox_td)

hr_naive <- exp(coef(fit_cox_naive)["plication_performed"])
ci_naive <- exp(confint(fit_cox_naive)["plication_performed", ])
hr_td    <- exp(coef(fit_cox_td)["plication_td"])
ci_td    <- exp(confint(fit_cox_td)["plication_td", ])

cat("SA4 Cox sensitivity analysis results:\n")
cat(sprintf("  Naive Cox HR      : %.2f (95%% CI: %.2f\u2013%.2f)\n",
            hr_naive, ci_naive[1], ci_naive[2]))
cat(sprintf("  Time-dependent HR : %.2f (95%% CI: %.2f\u2013%.2f)\n",
            hr_td, ci_td[1], ci_td[2]))
cat(sprintf("  PH assumption p   : %.3f\n\n",
            ph_test$table["plication_td", "p"]))


# ============================================================
# 結果サマリー
# ============================================================
cat("========================================\n")
cat("SENSITIVITY ANALYSES SUMMARY\n")
cat("========================================\n")
cat(sprintf("主解析  Median: 44 days (Bootstrap 95%% CI: 41\u201375 days)\n"))
cat(sprintf("SA1     Median: %.0f days（両側麻痺除外）\n",
            med_sa1$quantile))
cat(sprintf("SA2     Median: %.0f days（180日以内打ち切り除外）\n",
            med_sa2$quantile))
cat(sprintf("SA3 Left  KM : %.0f days\n", med_sa3_left$quantile))
cat(sprintf("SA3 Mid   KM : %.0f days\n", med_sa3_mid$quantile))
cat(sprintf("SA3 Right KM : %.0f days\n", med_sa3_right$quantile))
cat(sprintf("SA4 Cox TD HR: %.2f (95%% CI: %.2f\u2013%.2f)\n",
            hr_td, ci_td[1], ci_td[2]))
cat("========================================\n")


# ============================================================
# Figure: Supplement感度解析図（SA1・SA2の比較）
# ============================================================

# 主解析のTurnbull曲線
fit_main <- survfit(
  Surv(left_time, right_time, type = "interval2") ~ 1,
  data = df
)

make_tb_df <- function(fit, label, t_max = 1095) {
  df_out <- data.frame(
    time  = c(0, fit$time),
    surv  = c(1, fit$surv)
  ) %>%
    mutate(
      cumrec = 1 - surv,
      group  = label
    ) %>%
    bind_rows(
      data.frame(
        time   = t_max,
        surv   = tail(.$surv, 1),
        cumrec = tail(.$cumrec, 1),
        group  = label
      )
    )
  df_out
}

tb_main <- make_tb_df(fit_main, "Primary analysis (n=110)")
tb_sa1  <- make_tb_df(fit_sa1,  "SA1: Exclude bilateral (n=101)")
tb_sa2  <- make_tb_df(fit_sa2,  sprintf("SA2: Exclude early censored (n=%d)", nrow(df_sa2)))

tb_combined <- bind_rows(tb_main, tb_sa1, tb_sa2)

COL_MAIN <- "#1F3864"
COL_SA1  <- "#C0392B"
COL_SA2  <- "#2A9D8F"

fig_sa <- ggplot(tb_combined,
                 aes(x = time, y = cumrec * 100,
                     color = group, linetype = group)) +
  geom_step(linewidth = 0.9, direction = "hv") +
  scale_color_manual(
    values = c(
      "Primary analysis (n=110)"              = COL_MAIN,
      "SA1: Exclude bilateral (n=101)"        = COL_SA1,
      sprintf("SA2: Exclude early censored (n=%d)", nrow(df_sa2)) = COL_SA2
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Primary analysis (n=110)"              = "solid",
      "SA1: Exclude bilateral (n=101)"        = "dashed",
      sprintf("SA2: Exclude early censored (n=%d)", nrow(df_sa2)) = "dotted"
    )
  ) +
  scale_x_continuous(
    name   = "Time from diaphragm paralysis diagnosis (days)",
    limits = c(0, 1095),
    breaks = seq(0, 1095, by = 182),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    name   = "Cumulative recovery probability (%)",
    limits = c(0, 100),
    breaks = seq(0, 100, by = 25),
    expand = expansion(mult = c(0, 0.02))
  ) +
  theme_classic(base_family = "serif", base_size = 11) +
  theme(
    legend.position      = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.title         = element_blank(),
    legend.text          = element_text(size = 8.5),
    legend.background    = element_rect(fill = "white", color = "grey80",
                                        linewidth = 0.3),
    panel.grid           = element_blank(),
    plot.title           = element_text(size = 11, face = "bold"),
    plot.caption         = element_text(size = 7.5, color = "grey50",
                                        hjust = 0)
  ) +
  labs(
    title   = "Supplemental Figure. Sensitivity Analyses: Turnbull NPMLE",
    caption = paste0(
      "SA1: Patients with bilateral diaphragm paralysis excluded (n=9 excluded).\n",
      "SA2: Right-censored patients with follow-up \u2264180 days excluded.\n",
      "Primary analysis shown for reference."
    )
  )

print(fig_sa)

# 出力
output_dir <- "~/Desktop"

ggsave(
  filename = file.path(output_dir, "SuppFig_Sensitivity_SA1_SA2.pdf"),
  plot     = fig_sa,
  width    = 7.0,
  height   = 5.5,
  device   = "pdf"
)

ggsave(
  filename    = file.path(output_dir, "SuppFig_Sensitivity_SA1_SA2.tiff"),
  plot        = fig_sa,
  width       = 7.0,
  height      = 5.5,
  dpi         = 600,
  compression = "lzw"
)

cat("\n\u2713 07_sensitivity.R 完了\n")
cat("保存先:", output_dir, "\n")
cat("  - SuppFig_Sensitivity_SA1_SA2.pdf\n")
cat("  - SuppFig_Sensitivity_SA1_SA2.tiff\n\n")
