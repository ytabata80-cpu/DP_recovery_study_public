# ============================================================
#  07a_sensitivity_SA1_SA2.R  (JTCVS submission-compliant)
#  Supplemental Figure 2: Sensitivity Analyses SA1 & SA2
#  （Turnbull NPMLE、コホート定義違いによる3群比較）
#
#  Input : df（01_data_cleaning.R で作成）
#  Output: SuppFig_SA1_SA2.pdf / .tiff
#
#  [JTCVS対応で追加した点]
#  - X軸年単位表示（Figure 1・2と統一）
#  - No. followed代替指標を3群それぞれについて計算し、
#    risk panelに3行で表示（Turnbull由来のため標準at riskは
#    方法論的に定義不可、Figure 1・2と同じ理由・同じ代替指標）
#  - 3群のうち「最も早くat risk<10になる群」に合わせて
#    表示範囲を統一的に打ち切る
# ============================================================

library(survival)
library(ggplot2)
library(dplyr)
library(cowplot)

set.seed(1234)

if (!exists("df")) source("R/01_data_cleaning.R")

cat("=== 主解析（参照）===\n")
cat("n =", nrow(df), "\n\n")

# ============================================================
# 1. 3群のデータ定義とTurnbull推定
# ============================================================
df_main <- df
df_sa1  <- df %>% filter(paralysis_side != "Bilateral")
df_sa2  <- df %>%
  filter(recovery_confirmed == 1 |
         (recovery_confirmed == 0 & followup_days > 180))

fit_main <- survfit(Surv(left_time, right_time, type = "interval2") ~ 1, data = df_main)
fit_sa1  <- survfit(Surv(left_time, right_time, type = "interval2") ~ 1, data = df_sa1)
fit_sa2  <- survfit(Surv(left_time, right_time, type = "interval2") ~ 1, data = df_sa2)

med_main <- quantile(fit_main, probs = 0.5)
med_sa1  <- quantile(fit_sa1,  probs = 0.5)
med_sa2  <- quantile(fit_sa2,  probs = 0.5)

cat(sprintf("Primary : n=%d, Median=%.0f days\n", nrow(df_main), med_main$quantile))
cat(sprintf("SA1     : n=%d, Median=%.0f days (%d excluded: bilateral)\n",
            nrow(df_sa1), med_sa1$quantile, nrow(df) - nrow(df_sa1)))
cat(sprintf("SA2     : n=%d, Median=%.0f days (%d excluded: early right-censored)\n\n",
            nrow(df_sa2), med_sa2$quantile, nrow(df) - nrow(df_sa2)))

# ============================================================
# 2. No. followed代替指標（3群それぞれで計算）
#    Turnbullでは標準at riskが方法論的に定義不可なため、
#    「最終観察時点 >= t」を満たす患者数を代替指標とする
#    （Figure 1・2と同一ロジック）
# ============================================================
compute_last_obs <- function(d) {
  with(d, ifelse(is.na(right_time), left_time, right_time))
}
last_obs_main <- compute_last_obs(df_main)
last_obs_sa1  <- compute_last_obs(df_sa1)
last_obs_sa2  <- compute_last_obs(df_sa2)

# 表示範囲（X軸）は3年（1095日）で固定（Figure 1・2との比較容易性のため）。
# ただし各群の曲線・No. followed数値は、投稿規定の文言
# 「stop before at risk falls below 10 for that group」に従い、
# 群ごとに独立してat risk<10になる直前で打ち切る。
t_max <- 1095
step_days <- 30
scan_grid <- seq(0, t_max, by = step_days)

find_group_cutoff <- function(last_obs, grid) {
  n <- sapply(grid, function(t) sum(last_obs >= t))
  below10 <- which(n < 10)
  if (length(below10) > 0) grid[below10[1] - 1] else max(grid)
}
cutoff_main <- find_group_cutoff(last_obs_main, scan_grid)
cutoff_sa1  <- find_group_cutoff(last_obs_sa1,  scan_grid)
cutoff_sa2  <- find_group_cutoff(last_obs_sa2,  scan_grid)

cat(sprintf("群ごとの打ち切り時点（at risk>=10を満たす最大値）:\n"))
cat(sprintf("  Primary: %.0f days (%.2f years)\n", cutoff_main, cutoff_main/365.25))
cat(sprintf("  SA1    : %.0f days (%.2f years)\n", cutoff_sa1,  cutoff_sa1/365.25))
cat(sprintf("  SA2    : %.0f days (%.2f years)\n\n", cutoff_sa2,  cutoff_sa2/365.25))

risk_breaks_years <- seq(0, floor(t_max / 365.25 * 2) / 2, by = 0.5)
t_max_years <- t_max / 365.25
if (max(risk_breaks_years) < t_max_years - 0.01) {
  risk_breaks_years <- c(risk_breaks_years, round(t_max_years, 2))
}
risk_breaks_days <- pmin(risk_breaks_years * 365.25, t_max)
x_labels <- as.character(risk_breaks_years)

risk_table_long <- bind_rows(
  data.frame(group = "Primary analysis (n=110)",
             time_days = risk_breaks_days,
             n = ifelse(risk_breaks_days <= cutoff_main,
                        sapply(risk_breaks_days, function(t) sum(last_obs_main >= t)), NA)),
  data.frame(group = "SA1: Exclude bilateral (n=101)",
             time_days = risk_breaks_days,
             n = ifelse(risk_breaks_days <= cutoff_sa1,
                        sapply(risk_breaks_days, function(t) sum(last_obs_sa1 >= t)), NA)),
  data.frame(group = sprintf("SA2: Exclude early censored (n=%d)", nrow(df_sa2)),
             time_days = risk_breaks_days,
             n = ifelse(risk_breaks_days <= cutoff_sa2,
                        sapply(risk_breaks_days, function(t) sum(last_obs_sa2 >= t)), NA))
)
cat("=== No. followed（3群） ===\n")
print(risk_table_long)
cat("\n")

# ============================================================
# 3. Turnbull曲線データ
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

label_main <- "Primary analysis (n=110)"
label_sa1  <- "SA1: Exclude bilateral (n=101)"
label_sa2  <- sprintf("SA2: Exclude early censored (n=%d)", nrow(df_sa2))

tb_sa12 <- bind_rows(
  make_tb_df(fit_main, label_main, cutoff_main),
  make_tb_df(fit_sa1,  label_sa1,  cutoff_sa1),
  make_tb_df(fit_sa2,  label_sa2,  cutoff_sa2)
)

# ============================================================
# 4. Figure（cowplot 2段構成：曲線+X軸 / No. followed 3行）
# ============================================================
COL_MAIN <- "#1F3864"
COL_SA1  <- "#C0392B"
COL_SA2  <- "#2A9D8F"

shared_x_scale <- scale_x_continuous(
  limits = c(0, t_max),
  breaks = risk_breaks_days,
  labels = x_labels,
  expand = expansion(mult = c(0, 0.02))
)

fig_main <- ggplot(tb_sa12,
                   aes(x = time, y = cumrec * 100,
                       color = group, linetype = group,
                       linewidth = group)) +
  geom_step(direction = "hv") +
  scale_color_manual(values = setNames(
    c(COL_MAIN, COL_SA1, COL_SA2), c(label_main, label_sa1, label_sa2))) +
  scale_linetype_manual(values = setNames(
    c("solid", "dashed", "dotted"), c(label_main, label_sa1, label_sa2))) +
  scale_linewidth_manual(values = setNames(
    c(1.2, 1.0, 1.0), c(label_main, label_sa1, label_sa2)), guide = "none") +
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
    legend.position      = c(0.98, 0.08),
    legend.justification = c(1, 0),
    legend.title          = element_blank(),
    legend.text            = element_text(size = 8.5, family = "sans"),
    legend.key.width        = unit(1.5, "cm"),
    legend.background         = element_rect(fill = "white", color = "grey80",
                                             linewidth = 0.3),
    panel.grid                  = element_blank(),
    plot.margin                   = margin(10, 15, 2, 10)
  )

# --- No. followed 3行パネル ---
group_order <- c(label_main, label_sa1, label_sa2)
risk_table_long$group <- factor(risk_table_long$group, levels = rev(group_order))
risk_table_long$y_pos  <- as.numeric(risk_table_long$group)

risk_panel <- ggplot(risk_table_long, aes(x = time_days, y = y_pos)) +
  geom_text(aes(label = n), size = 2.8, family = "sans", color = "grey20") +
  scale_y_continuous(
    limits = c(0.4, 3.6), breaks = 1:3,
    labels = rev(c("Primary", "SA1", "SA2")),
    name = "No. followed*"
  ) +
  shared_x_scale +
  coord_cartesian(clip = "off") +
  theme_classic(base_family = "sans", base_size = 11) +
  theme(
    axis.title.y  = element_text(size = 9, color = "grey20",
                                 angle = 0, vjust = 0.5, hjust = 1,
                                 margin = margin(r = 8)),
    axis.text.y   = element_text(size = 7.5, color = "grey20",
                                 margin = margin(r = 6)),
    axis.ticks.y  = element_blank(),
    axis.line.y   = element_blank(),
    axis.title.x  = element_blank(),
    axis.text.x   = element_blank(),
    axis.line.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    panel.grid      = element_blank(),
    plot.margin       = margin(2, 15, 10, 35)
  )

fig_sa12 <- plot_grid(
  fig_main, risk_panel,
  ncol = 1, align = "v", axis = "lr",
  rel_heights = c(4, 1.3)
)

print(fig_sa12)

# ============================================================
# 5. 出力
# ============================================================
output_dir <- "~/Desktop"

ggsave(file.path(output_dir, "SuppFig_SA1_SA2.pdf"),
       plot = fig_sa12, width = 7.0, height = 6.0, device = "pdf")
ggsave(file.path(output_dir, "SuppFig_SA1_SA2.tiff"),
       plot = fig_sa12, width = 7.0, height = 6.0,
       dpi = 600, compression = "lzw")

cat("\u2713 SuppFig_SA1_SA2 完了（JTCVS対応: 年単位X軸 + No. followed 3行）\n")
