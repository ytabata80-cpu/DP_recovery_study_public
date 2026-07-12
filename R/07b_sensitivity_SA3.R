# ============================================================
#  07b_sensitivity_SA3.R  (JTCVS submission-compliant)
#  Supplemental Figure 3: SA3 KM Three-Point Imputation
#  （Turnbull参照 + 3種のKM補完法、計4群比較）
#
#  Input : df（01_data_cleaning.R で作成）
#  Output: SuppFig_SA3_KM_imputation.pdf / .tiff
#
#  [JTCVS対応で追加した点]
#  - X軸年単位表示
#  - Turnbull（参照曲線）: Figure 1・2と同じNo. followed代替指標
#  - KM 3補完（left/mid/right endpoint）: 全例が単一の
#    imputed event/censoring timeを持つため、標準的な
#    survfit()のnumber-at-risk（summary()$n.risk）がそのまま
#    使える（Turnbullのような方法論的制約はない）
#  - 4群のうち最も早くat risk<10になる群に合わせて表示範囲を統一
# ============================================================

library(survival)
library(ggplot2)
library(dplyr)
library(cowplot)

set.seed(1234)

if (!exists("df")) source("R/01_data_cleaning.R")

# ============================================================
# 1. Turnbull（参照）
# ============================================================
fit_main <- survfit(Surv(left_time, right_time, type = "interval2") ~ 1, data = df)

# ============================================================
# 2. KM 3補完
# ============================================================
df_sa3_left <- df %>%
  mutate(time_km = ifelse(!is.na(right_time), left_time, followup_days),
         event = recovery_confirmed)
fit_sa3_left <- survfit(Surv(time_km, event) ~ 1, data = df_sa3_left)

df_sa3_mid <- df %>%
  mutate(time_km = ifelse(!is.na(right_time),
                          (left_time + right_time) / 2, followup_days),
         event = recovery_confirmed)
fit_sa3_mid <- survfit(Surv(time_km, event) ~ 1, data = df_sa3_mid)

df_sa3_right <- df %>%
  mutate(time_km = ifelse(!is.na(right_time), right_time, followup_days),
         event = recovery_confirmed)
fit_sa3_right <- survfit(Surv(time_km, event) ~ 1, data = df_sa3_right)

med_left  <- quantile(fit_sa3_left,  probs = 0.5)
med_mid   <- quantile(fit_sa3_mid,   probs = 0.5)
med_right <- quantile(fit_sa3_right, probs = 0.5)
cat(sprintf("Left endpoint  median: %.0f days\n", med_left$quantile))
cat(sprintf("Mid endpoint   median: %.0f days\n", med_mid$quantile))
cat(sprintf("Right endpoint median: %.0f days\n\n", med_right$quantile))

# ============================================================
# 3. No. followed（4群：Turnbull=代替指標、KM3種=標準n.risk）
# ============================================================
# 表示範囲（X軸）は3年（1095日）で固定（Figure 1・2との比較容易性のため）。
# 各群の曲線・数値は、投稿規定の文言に従い群ごとに独立して
# at risk<10になる直前で打ち切る。
t_max <- 1095
step_days <- 30
scan_grid <- seq(0, t_max, by = step_days)

last_obs_main <- with(df, ifelse(is.na(right_time), left_time, right_time))

get_km_nrisk <- function(fit, times) {
  s <- summary(fit, times = times, extend = TRUE)
  s$n.risk
}

find_group_cutoff_alt <- function(last_obs, grid) {
  n <- sapply(grid, function(t) sum(last_obs >= t))
  below10 <- which(n < 10)
  if (length(below10) > 0) grid[below10[1] - 1] else max(grid)
}
find_group_cutoff_km <- function(fit, grid) {
  n <- get_km_nrisk(fit, grid)
  below10 <- which(n < 10)
  if (length(below10) > 0) grid[below10[1] - 1] else max(grid)
}

cutoff_turnbull <- find_group_cutoff_alt(last_obs_main, scan_grid)
cutoff_left     <- find_group_cutoff_km(fit_sa3_left,  scan_grid)
cutoff_mid      <- find_group_cutoff_km(fit_sa3_mid,   scan_grid)
cutoff_right    <- find_group_cutoff_km(fit_sa3_right, scan_grid)

cat("群ごとの打ち切り時点（at risk>=10を満たす最大値）:\n")
cat(sprintf("  Turnbull : %.0f days (%.2f years)\n", cutoff_turnbull, cutoff_turnbull/365.25))
cat(sprintf("  KM-left  : %.0f days (%.2f years)\n", cutoff_left,     cutoff_left/365.25))
cat(sprintf("  KM-mid   : %.0f days (%.2f years)\n", cutoff_mid,      cutoff_mid/365.25))
cat(sprintf("  KM-right : %.0f days (%.2f years)\n\n", cutoff_right,    cutoff_right/365.25))

risk_breaks_years <- seq(0, floor(t_max / 365.25 * 2) / 2, by = 0.5)
t_max_years <- t_max / 365.25
if (max(risk_breaks_years) < t_max_years - 0.01) {
  risk_breaks_years <- c(risk_breaks_years, round(t_max_years, 2))
}
risk_breaks_days <- pmin(risk_breaks_years * 365.25, t_max)
x_labels <- as.character(risk_breaks_years)

label_turnbull <- sprintf("Turnbull NPMLE (median=%.0fd)", quantile(fit_main, 0.5)$quantile)
label_left     <- sprintf("KM Left endpoint (median=%.0fd)", med_left$quantile)
label_mid      <- sprintf("KM Mid endpoint (median=%.0fd)",  med_mid$quantile)
label_right    <- sprintf("KM Right endpoint (median=%.0fd)", med_right$quantile)

risk_table_long <- bind_rows(
  data.frame(group = label_turnbull, time_days = risk_breaks_days,
             n = ifelse(risk_breaks_days <= cutoff_turnbull,
                        sapply(risk_breaks_days, function(t) sum(last_obs_main >= t)), NA)),
  data.frame(group = label_left, time_days = risk_breaks_days,
             n = ifelse(risk_breaks_days <= cutoff_left,
                        get_km_nrisk(fit_sa3_left, risk_breaks_days), NA)),
  data.frame(group = label_mid, time_days = risk_breaks_days,
             n = ifelse(risk_breaks_days <= cutoff_mid,
                        get_km_nrisk(fit_sa3_mid, risk_breaks_days), NA)),
  data.frame(group = label_right, time_days = risk_breaks_days,
             n = ifelse(risk_breaks_days <= cutoff_right,
                        get_km_nrisk(fit_sa3_right, risk_breaks_days), NA))
)
cat("=== No. followed / at risk（4群） ===\n")
print(risk_table_long)
cat("\n")

# ============================================================
# 4. 曲線データ
# ============================================================
make_tb_df <- function(fit, label, t_max) {
  data.frame(time = c(0, fit$time), surv = c(1, fit$surv)) %>%
    mutate(cumrec = 1 - surv, group = label) %>%
    bind_rows(data.frame(
      time = t_max, surv = tail(.$surv, 1),
      cumrec = tail(.$cumrec, 1), group = label
    )) %>%
    filter(time <= t_max)
}

tb_sa3 <- bind_rows(
  make_tb_df(fit_main,      label_turnbull, cutoff_turnbull),
  make_tb_df(fit_sa3_left,  label_left,     cutoff_left),
  make_tb_df(fit_sa3_mid,   label_mid,      cutoff_mid),
  make_tb_df(fit_sa3_right, label_right,    cutoff_right)
)

# ============================================================
# 5. Figure（cowplot 2段構成）
# ============================================================
COL_TURNBULL <- "#1F3864"
COL_LEFT     <- "#C0392B"
COL_MID      <- "#E07B39"
COL_RIGHT    <- "#2A9D8F"

shared_x_scale <- scale_x_continuous(
  limits = c(0, t_max),
  breaks = risk_breaks_days,
  labels = x_labels,
  expand = expansion(mult = c(0, 0.02))
)

fig_main <- ggplot(tb_sa3,
                   aes(x = time, y = cumrec * 100,
                       color = group, linetype = group,
                       linewidth = group)) +
  geom_step(direction = "hv") +
  geom_hline(yintercept = 50, linetype = "longdash",
             color = "grey40", linewidth = 0.5) +
  scale_color_manual(values = setNames(
    c(COL_TURNBULL, COL_LEFT, COL_MID, COL_RIGHT),
    c(label_turnbull, label_left, label_mid, label_right))) +
  scale_linetype_manual(values = setNames(
    c("solid", "dashed", "dotted", "dotdash"),
    c(label_turnbull, label_left, label_mid, label_right))) +
  scale_linewidth_manual(values = setNames(
    c(1.2, 0.9, 0.9, 0.9),
    c(label_turnbull, label_left, label_mid, label_right)), guide = "none") +
  shared_x_scale +
  scale_y_continuous(
    name = "Cumulative recovery probability (%)",
    limits = c(0, 100), breaks = seq(0, 100, by = 25),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = "Time from diaphragm paralysis diagnosis (years)") +
  theme_classic(base_family = "sans", base_size = 11) +
  theme(
    axis.title.x      = element_text(size = 10),
    axis.text.x       = element_text(size = 9, color = "grey20"),
    axis.ticks.x      = element_line(color = "grey40"),
    axis.title.y      = element_text(size = 10),
    axis.text.y       = element_text(size = 9, color = "grey20"),
    axis.line.y       = element_line(color = "black", linewidth = 0.5),
    axis.line.x       = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y      = element_line(color = "grey40"),
    legend.position      = c(0.98, 0.06),
    legend.justification = c(1, 0),
    legend.title          = element_blank(),
    legend.text            = element_text(size = 8.0, family = "sans"),
    legend.key.width        = unit(1.5, "cm"),
    legend.background         = element_rect(fill = "white", color = "grey80",
                                             linewidth = 0.3),
    panel.grid                  = element_blank(),
    plot.margin                   = margin(10, 15, 2, 10)
  )

# --- No. followed / at risk 4行パネル ---
group_order <- c(label_turnbull, label_left, label_mid, label_right)
risk_table_long$group <- factor(risk_table_long$group, levels = rev(group_order))
risk_table_long$y_pos  <- as.numeric(risk_table_long$group)

risk_panel <- ggplot(risk_table_long, aes(x = time_days, y = y_pos)) +
  geom_text(aes(label = n), size = 2.6, family = "sans", color = "grey20") +
  scale_y_continuous(
    limits = c(0.4, 4.6), breaks = 1:4,
    labels = rev(c("Turnbull", "KM-left", "KM-mid", "KM-right")),
    name = "No. followed /\nat risk*"
  ) +
  shared_x_scale +
  coord_cartesian(clip = "off") +
  theme_classic(base_family = "sans", base_size = 11) +
  theme(
    axis.title.y  = element_text(size = 8, color = "grey20",
                                 angle = 0, vjust = 0.5, hjust = 1,
                                 margin = margin(r = 10)),
    axis.text.y   = element_text(size = 7, color = "grey20",
                                 margin = margin(r = 10)),
    axis.ticks.y  = element_blank(),
    axis.line.y   = element_blank(),
    axis.title.x  = element_blank(),
    axis.text.x   = element_blank(),
    axis.line.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    panel.grid      = element_blank(),
    plot.margin       = margin(2, 15, 10, 35)
  )

fig_sa3 <- plot_grid(
  fig_main, risk_panel,
  ncol = 1, align = "v", axis = "lr",
  rel_heights = c(4, 1.6)
)

print(fig_sa3)

# ============================================================
# 6. 出力
# ============================================================
output_dir <- "~/Desktop"

ggsave(file.path(output_dir, "SuppFig_SA3_KM_imputation.pdf"),
       plot = fig_sa3, width = 7.0, height = 6.3, device = "pdf")
ggsave(file.path(output_dir, "SuppFig_SA3_KM_imputation.tiff"),
       plot = fig_sa3, width = 7.0, height = 6.3,
       dpi = 600, compression = "lzw")

cat("\u2713 SuppFig_SA3 完了（JTCVS対応: 年単位X軸 + No. followed/at risk 4行）\n")
