# ============================================================
#  04_turnbull_npmle.R
#  Figure 1: Turnbull NPMLE + Bootstrap 95% CI
#            overlaid with Parametric AFT Models
#            (Weibull / Log-normal / Log-logistic)
#
#  Input : df（01_data_cleaning.R で作成）
#  Output: Figure1_DP_Recovery_Turnbull_AFT.pdf / .tiff
#          ~/Desktop/ に保存
#
#  実行方法:
#    source("R/01_data_cleaning.R")
#    source("R/04_turnbull_npmle.R")
#
#  設計上の注意点:
#  - left_time = 0（術当日診断例）は
#    Turnbull・Bootstrapはそのまま使用
#  - AFTモデル比較時のみ left_time = 0 → 0.5 に変換（df_aft）
#  - Log-logisticが選択されたため最終AFTもdf_aftで実行
#  - Bootstrap用データ（df_boot）は right_time=NA かつ
#    left_time > 1095 の症例を1095日に制限
#  - 表示上限: t_max = 1095日（3年）
#  - P0082（2527日で回復確認）は解析に含むが表示範囲外
# ============================================================

library(survival)
library(flexsurv)
library(ggplot2)
library(dplyr)

set.seed(1234)

# 01_data_cleaning.R が未実行の場合は自動実行
if (!exists("df")) {
  source("R/01_data_cleaning.R")
}

# ============================================================
# 0. 解析用データ確認
# ============================================================
cat("=== 解析対象 ===\n")
cat("n =", nrow(df), "\n")
cat("回復確認例:", sum(df$recovery_confirmed == 1, na.rm=TRUE), "\n")
cat("右側打ち切り例:", sum(df$recovery_confirmed == 0, na.rm=TRUE), "\n")
cat("left_time  範囲:", min(df$left_time), "–",
    max(df$left_time), "days\n")
cat("right_time 範囲:", min(df$right_time, na.rm=TRUE), "–",
    max(df$right_time, na.rm=TRUE), "days\n\n")

# 表示上限: 3年（1095日）
t_max     <- 1095
time_grid <- seq(0, t_max, by = 1)
cat(sprintf("time_grid: 0 – %d days (display limit: 3 years)\n\n",
            t_max))


# ============================================================
# 1. Turnbull NPMLE（主解析）
# ============================================================
surv_obj <- Surv(
  time  = df$left_time,
  time2 = df$right_time,
  type  = "interval2"
)
fit_km <- survfit(surv_obj ~ 1)

# Turnbull曲線のデータフレーム化
tb_df <- data.frame(
  time   = c(0, fit_km$time),
  surv   = c(1, fit_km$surv),
  lower  = c(1, fit_km$lower),
  upper  = c(1, fit_km$upper)
) %>%
  mutate(
    cumrec       = 1 - surv,
    cumrec_lower = 1 - upper,
    cumrec_upper = 1 - lower
  )

# 1095日まで延長（最後の値を水平に延長）
tb_df_extended <- tb_df %>%
  bind_rows(
    data.frame(
      time         = t_max,
      surv         = tail(tb_df$surv, 1),
      lower        = tail(tb_df$lower, 1),
      upper        = tail(tb_df$upper, 1),
      cumrec       = tail(tb_df$cumrec, 1),
      cumrec_lower = tail(tb_df$cumrec_lower, 1),
      cumrec_upper = tail(tb_df$cumrec_upper, 1)
    )
  )

# 凡例統合用にmodel列追加
tb_df_line <- tb_df_extended %>%
  mutate(model = "Turnbull NPMLE")


# ============================================================
# 2. Bootstrap 95% CI（2000 iterations）
# ============================================================
# right_time=NA かつ left_time > 1095 の症例は1095日に制限
df_boot <- df %>%
  mutate(
    left_time = ifelse(
      is.na(right_time) & left_time > t_max,
      t_max, left_time
    )
  )

cat("Bootstrap CI計算中（2000回）...\n")
B        <- 2000
boot_mat <- matrix(NA, nrow = B, ncol = length(time_grid))

for (b in seq_len(B)) {
  idx  <- sample(nrow(df_boot), replace = TRUE)
  df_b <- df_boot[idx, ]
  tryCatch({
    fit_b <- survfit(
      Surv(left_time, right_time, type = "interval2") ~ 1,
      data = df_b
    )
    S_b <- stepfun(fit_b$time, c(1, fit_b$surv))
    boot_mat[b, ] <- 1 - S_b(time_grid)
  }, error = function(e) NULL)
}

boot_lo  <- apply(boot_mat, 2, quantile, probs = 0.025, na.rm = TRUE)
boot_hi  <- apply(boot_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
boot_med <- apply(boot_mat, 2, median,   na.rm = TRUE)

boot_df <- data.frame(
  time       = time_grid,
  cumrec     = boot_med,
  boot_lower = boot_lo,
  boot_upper = boot_hi
)
cat("Bootstrap完了\n\n")


# ============================================================
# 3. 中央値回復時間
# ============================================================
# Turnbull点推定
med_est <- quantile(fit_km, probs = 0.5)
med_val <- round(med_est$quantile)

# Bootstrap由来の95% CI（50%到達時点）
med_lo_boot <- time_grid[which.min(abs(boot_hi - 0.5))]
med_hi_boot <- time_grid[which.min(abs(boot_lo - 0.5))]

cat(sprintf(
  "Median recovery: %d days (95%% bootstrap CI: %d\u2013%d days)\n\n",
  med_val, med_lo_boot, med_hi_boot))


# ============================================================
# 4. Parametric AFT Models
#    left_time = 0 → 0.5 に変換（モデル比較・最終推定用）
# ============================================================
df_aft <- df %>%
  mutate(
    left_time  = ifelse(left_time == 0, 0.5, left_time),
    right_time = ifelse(!is.na(right_time) & right_time == 0,
                        0.5, right_time)
  )

cat("AFTモデル推定中...\n")
fit_wei <- flexsurvreg(
  Surv(left_time, right_time, type = "interval2") ~ 1,
  data = df_aft, dist = "weibull"
)
fit_ln <- flexsurvreg(
  Surv(left_time, right_time, type = "interval2") ~ 1,
  data = df_aft, dist = "lnorm"
)
fit_ll_final <- flexsurvreg(
  Surv(left_time, right_time, type = "interval2") ~ 1,
  data = df_aft, dist = "llogis"
)

# AIC / BIC 比較
aic_tbl <- data.frame(
  Model = c("Weibull", "Log-normal", "Log-logistic"),
  AIC   = c(AIC(fit_wei), AIC(fit_ln),  AIC(fit_ll_final)),
  BIC   = c(BIC(fit_wei), BIC(fit_ln),  BIC(fit_ll_final))
) %>%
  mutate(
    delta_AIC = AIC - min(AIC),
    delta_BIC = BIC - min(BIC)
  ) %>%
  arrange(AIC)

cat("\n=== Model Selection (AIC / BIC) ===\n")
print(aic_tbl)
best_model <- aic_tbl$Model[1]
cat(sprintf("\nSelected model (lowest AIC): %s\n\n", best_model))

# AFT fitted curves
make_aft_df <- function(fit, label, time_grid) {
  S <- summary(fit, t = time_grid, type = "survival")[[1]]
  data.frame(
    time   = time_grid,
    cumrec = 1 - S$est,
    model  = label
  )
}

aft_df <- bind_rows(
  make_aft_df(fit_wei,      "Weibull",      time_grid),
  make_aft_df(fit_ln,       "Log-normal",   time_grid),
  make_aft_df(fit_ll_final, "Log-logistic", time_grid)
) %>%
  mutate(selected = model == best_model)


# ============================================================
# 5. Figure 1（JTCVS対応版）
# ============================================================
COL_TURNBULL <- "#1F3864"
COL_SHADE    <- "#1F3864"
COL_WEI      <- "#E07B39"
COL_LN       <- "#2A9D8F"
COL_LL       <- "#C0392B"

legend_order <- c("Turnbull NPMLE", "Log-logistic",
                  "Log-normal", "Weibull")

all_colors <- c(
  "Turnbull NPMLE" = COL_TURNBULL,
  "Log-logistic"   = COL_LL,
  "Log-normal"     = COL_LN,
  "Weibull"        = COL_WEI
)
all_ltypes <- c(
  "Turnbull NPMLE" = "solid",
  "Log-logistic"   = "dotdash",
  "Log-normal"     = "dotted",
  "Weibull"        = "dashed"
)
all_labels <- c(
  "Turnbull NPMLE" = "Turnbull NPMLE (Bootstrap 95% CI)",
  "Log-logistic"   = sprintf("Log-logistic (AIC = %.1f)",
                             AIC(fit_ll_final)),
  "Log-normal"     = sprintf("Log-normal (AIC = %.1f)",
                             AIC(fit_ln)),
  "Weibull"        = sprintf("Weibull (AIC = %.1f)",
                             AIC(fit_wei))
)

annot_text <- sprintf(
  "Median recovery: %d days\n(95%% bootstrap CI: %d\u2013%d days)",
  med_val, med_lo_boot, med_hi_boot
)

fig1 <- ggplot() +

  geom_ribbon(
    data = boot_df,
    aes(x = time, ymin = boot_lower * 100, ymax = boot_upper * 100),
    fill = COL_SHADE, alpha = 0.15, show.legend = FALSE
  ) +

  geom_line(
    data = aft_df,
    aes(x = time, y = cumrec * 100,
        color = model, linetype = model,
        linewidth = selected),
    alpha = 0.85
  ) +
  scale_linewidth_manual(
    values = c("TRUE" = 0.7, "FALSE" = 0.55),
    guide  = "none"
  ) +

  geom_step(
    data = tb_df_line,
    aes(x = time, y = cumrec * 100,
        color = model, linetype = model),
    linewidth = 1.0, direction = "hv"
  ) +

  geom_vline(xintercept = med_val,
             linetype = "longdash", color = "grey40",
             linewidth = 0.5) +
  geom_hline(yintercept = 50,
             linetype = "longdash", color = "grey40",
             linewidth = 0.5) +

  annotate("label",
           x = 200, y = 35,
           label = annot_text,
           size = 3.0, hjust = 0,
           fill = "white", color = "grey30",
           label.size = 0.3, family = "sans") +   # ← serif から sans に変更

  scale_color_manual(
    name   = NULL,
    values = all_colors,
    labels = all_labels,
    breaks = legend_order
  ) +
  scale_linetype_manual(
    name   = NULL,
    values = all_ltypes,
    labels = all_labels,
    breaks = legend_order
  ) +
  scale_x_continuous(
    name   = "Time from diaphragm paralysis diagnosis (days)",
    limits = c(0, t_max),
    breaks = seq(0, t_max, by = 182),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    name   = "Cumulative recovery probability (%)",
    limits = c(0, 100),
    breaks = seq(0, 100, by = 25),
    expand = expansion(mult = c(0, 0.02))
  ) +

  # JTCVS-style theme（sans-serif、タイトル・サブタイトル・キャプションなし）
  theme_classic(base_family = "sans", base_size = 11) +
  theme(
    axis.title        = element_text(size = 10),
    axis.text         = element_text(size = 9, color = "grey20"),
    axis.line         = element_line(color = "black", linewidth = 0.5),
    axis.ticks        = element_line(color = "grey40"),
    legend.position   = c(0.98, 0.25),
    legend.justification = c(1, 0),
    legend.title      = element_blank(),
    legend.text       = element_text(size = 8.5, family = "sans"),
    legend.key.width  = unit(1.8, "cm"),
    legend.background = element_rect(fill = "white", color = "grey80",
                                     linewidth = 0.3),
    panel.grid        = element_blank(),
    plot.margin       = margin(10, 15, 10, 10)
  )
  # labs(title=..., subtitle=..., caption=...) は削除
  # → Figure Legendは本文のFigure Legends欄に記載

print(fig1)

# ============================================================
# 6. 出力（PDF + TIFF 600dpi）
# ============================================================
output_dir <- "~/Desktop"

# PDF（確認・共著者レビュー用）
ggsave(
  filename = file.path(output_dir,
                       "Figure1_DP_Recovery_Turnbull_AFT.pdf"),
  plot     = fig1,
  width    = 7.0,
  height   = 5.5,
  device   = "pdf"
)

# TIFF 600dpi LZW圧縮（投稿用）
ggsave(
  filename    = file.path(output_dir,
                          "Figure1_DP_Recovery_Turnbull_AFT.tiff"),
  plot        = fig1,
  width       = 7.0,
  height      = 5.5,
  dpi         = 600,
  compression = "lzw"
)

cat("\n\u2713 Figure 1 完成\n")
cat("保存先:", output_dir, "\n")
cat("  - Figure1_DP_Recovery_Turnbull_AFT.pdf\n")
cat("  - Figure1_DP_Recovery_Turnbull_AFT.tiff (600 dpi, LZW)\n\n")


# ============================================================
# 7. Results記述用サマリー
# ============================================================
cat("========================================\n")
cat("RESULTS SUMMARY\n")
cat("========================================\n")
cat(sprintf("N                        : %d\n", nrow(df)))
cat(sprintf("Recovered                : %d (%.1f%%)\n",
            sum(df$recovery_confirmed == 1, na.rm=TRUE),
            sum(df$recovery_confirmed == 1, na.rm=TRUE) /
              nrow(df) * 100))
cat(sprintf("Right-censored           : %d (%.1f%%)\n",
            sum(df$recovery_confirmed == 0, na.rm=TRUE),
            sum(df$recovery_confirmed == 0, na.rm=TRUE) /
              nrow(df) * 100))
cat(sprintf("Median recovery          : %d days\n", med_val))
cat(sprintf("95%% bootstrap CI         : %d\u2013%d days\n",
            med_lo_boot, med_hi_boot))
cat(sprintf("Best-fit AFT model (AIC) : %s (AIC = %.1f)\n",
            best_model, min(aic_tbl$AIC)))
cat("\n=== Cumulative Recovery Rates ===\n")
for (tp in c(30, 60, 90, 180, 365, 730, 1095)) {
  idx <- which.min(abs(time_grid - tp))
  cat(sprintf("%4d days: %5.1f%% (95%% CI: %5.1f%%\u2013%5.1f%%)\n",
              tp,
              round(boot_med[idx] * 100, 1),
              round(boot_lo[idx]  * 100, 1),
              round(boot_hi[idx]  * 100, 1)))
}
cat("========================================\n")  mutate(
    cumrec       = 1 - surv,
    cumrec_lower = 1 - upper,
    cumrec_upper = 1 - lower
  )


# ============================================================
# 2. Bootstrap 95% CI（2000 iterations）
# ============================================================
cat("Bootstrap CI計算中（2000回）...\n")
B        <- 2000
boot_mat <- matrix(NA, nrow = B, ncol = length(time_grid))

for (b in seq_len(B)) {
  idx  <- sample(nrow(df), replace = TRUE)
  df_b <- df[idx, ]
  tryCatch({
    fit_b <- survfit(
      Surv(left_time, right_time, type = "interval2") ~ 1,
      data = df_b
    )
    S_b <- stepfun(fit_b$time, c(1, fit_b$surv))
    boot_mat[b, ] <- 1 - S_b(time_grid)
  }, error = function(e) NULL)
}

boot_lo  <- apply(boot_mat, 2, quantile, probs = 0.025, na.rm = TRUE)
boot_hi  <- apply(boot_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
boot_med <- apply(boot_mat, 2, median,   na.rm = TRUE)

boot_df <- data.frame(
  time       = time_grid,
  cumrec     = boot_med,
  boot_lower = boot_lo,
  boot_upper = boot_hi
)
cat("Bootstrap完了\n\n")


# ============================================================
# 3. 中央値回復時間（Turnbull）
# ============================================================
med_est <- quantile(fit_km, probs = 0.5)
med_val <- round(med_est$quantile)
med_lo  <- round(med_est$lower)
med_hi  <- round(med_est$upper)
cat(sprintf("Median recovery: %d days (95%% CI: %d\u2013%d)\n\n",
            med_val, med_lo, med_hi))


# ============================================================
# 4. Parametric AFT Models
#    Weibull / Log-normal / Log-logistic
# ============================================================
cat("AFTモデル推定中...\n")
fit_wei <- flexsurvreg(
  Surv(left_time, right_time, type = "interval2") ~ 1,
  data = df, dist = "weibull"
)
fit_ln <- flexsurvreg(
  Surv(left_time, right_time, type = "interval2") ~ 1,
  data = df, dist = "lnorm"
)
fit_ll <- flexsurvreg(
  Surv(left_time, right_time, type = "interval2") ~ 1,
  data = df, dist = "llogis"
)

# AIC / BIC 比較
aic_tbl <- data.frame(
  Model = c("Weibull", "Log-normal", "Log-logistic"),
  AIC   = c(AIC(fit_wei), AIC(fit_ln),  AIC(fit_ll)),
  BIC   = c(BIC(fit_wei), BIC(fit_ln),  BIC(fit_ll))
) %>%
  mutate(
    delta_AIC = AIC - min(AIC),
    delta_BIC = BIC - min(BIC)
  ) %>%
  arrange(AIC)

cat("\n=== Model Selection (AIC / BIC) ===\n")
print(aic_tbl)
best_model <- aic_tbl$Model[1]
cat(sprintf("Selected model (lowest AIC): %s\n\n", best_model))

# AFT fitted curves
make_aft_df <- function(fit, label, time_grid) {
  S <- summary(fit, t = time_grid, type = "survival")[[1]]
  data.frame(
    time   = time_grid,
    cumrec = 1 - S$est,
    model  = label
  )
}

aft_df <- bind_rows(
  make_aft_df(fit_wei, "Weibull",      time_grid),
  make_aft_df(fit_ln,  "Log-normal",   time_grid),
  make_aft_df(fit_ll,  "Log-logistic", time_grid)
) %>%
  mutate(selected = model == best_model)


# ============================================================
# 5. Figure 1
# ============================================================
COL_TURNBULL <- "#1F3864"
COL_SHADE    <- "#1F3864"
COL_WEI      <- "#E07B39"
COL_LN       <- "#2A9D8F"
COL_LL       <- "#C0392B"

model_colors <- setNames(
  c(COL_WEI, COL_LN, COL_LL),
  c("Weibull", "Log-normal", "Log-logistic")
)
model_labels <- setNames(
  sprintf("%s (AIC = %.1f)", aic_tbl$Model, aic_tbl$AIC),
  aic_tbl$Model
)
model_ltypes <- c(
  "Weibull"      = "dashed",
  "Log-normal"   = "dotted",
  "Log-logistic" = "dotdash"
)

annot_text <- sprintf(
  "Median recovery: %d days\n(95%% bootstrap CI: %d\u2013%d days)",
  med_val, med_lo, med_hi
)

# x軸のbreaksを実データ範囲に合わせて設定
x_breaks <- seq(0, t_max, by = ifelse(t_max > 365, 56, 28))

fig1 <- ggplot() +

  # Bootstrap 95% CI ribbon
  geom_ribbon(
    data = boot_df,
    aes(x = time, ymin = boot_lower * 100, ymax = boot_upper * 100),
    fill = COL_SHADE, alpha = 0.15
  ) +

  # Parametric AFT curves
  geom_line(
    data = aft_df,
    aes(x = time, y = cumrec * 100,
        color = model, linetype = model,
        linewidth = selected),
    alpha = 0.85
  ) +
  scale_linewidth_manual(
    values = c("TRUE" = 1.0, "FALSE" = 0.75),
    guide  = "none"
  ) +

  # Turnbull NPMLE step curve
  geom_step(
    data = tb_df,
    aes(x = time, y = cumrec * 100),
    color = COL_TURNBULL, linewidth = 1.4, direction = "hv"
  ) +

  # Median reference lines
  geom_vline(xintercept = med_val,
             linetype = "longdash", color = "grey40", linewidth = 0.6) +
  geom_hline(yintercept = 50,
             linetype = "longdash", color = "grey40", linewidth = 0.6) +

  # Median annotation
  annotate("label",
           x = med_val + t_max * 0.03, y = 12,
           label = annot_text,
           size = 3.0, hjust = 0,
           fill = "white", color = "grey30",
           label.size = 0.3, family = "serif") +

  # Manual legend: Turnbull
  annotate("segment",
           x = 0, xend = t_max * 0.06,
           y = 97, yend = 97,
           color = COL_TURNBULL, linewidth = 1.4) +
  annotate("rect",
           xmin = 0, xmax = t_max * 0.06,
           ymin = 93, ymax = 101,
           fill = COL_SHADE, alpha = 0.15, color = NA) +
  annotate("text",
           x = t_max * 0.065, y = 97,
           label = "Turnbull NPMLE (Bootstrap 95% CI)",
           hjust = 0, size = 3.0,
           color = "grey20", family = "serif") +

  # Scales
  scale_color_manual(
    name   = "Parametric AFT models",
    values = model_colors,
    labels = model_labels
  ) +
  scale_linetype_manual(
    name   = "Parametric AFT models",
    values = model_ltypes,
    labels = model_labels
  ) +
  scale_x_continuous(
    name   = "Time from diaphragm paralysis diagnosis (days)",
    limits = c(0, t_max),
    breaks = x_breaks,
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    name   = "Cumulative recovery probability (%)",
    limits = c(0, 100),
    breaks = seq(0, 100, by = 25),
    expand = expansion(mult = c(0, 0.02))
  ) +

  # JAMA-style theme
  theme_classic(base_family = "serif", base_size = 11) +
  theme(
    axis.title        = element_text(size = 10),
    axis.text         = element_text(size = 9, color = "grey20"),
    axis.line         = element_line(color = "black", linewidth = 0.5),
    axis.ticks        = element_line(color = "grey40"),
    legend.position   = c(0.98, 0.20),
    legend.justification = c(1, 0),
    legend.title      = element_text(size = 9, face = "bold"),
    legend.text       = element_text(size = 8.5),
    legend.key.width  = unit(1.8, "cm"),
    legend.background = element_rect(fill = "white", color = "grey80",
                                     linewidth = 0.3),
    panel.grid        = element_blank(),
    plot.title        = element_text(size = 11, face = "bold"),
    plot.subtitle     = element_text(size = 9, color = "grey40"),
    plot.caption      = element_text(size = 7.5, color = "grey50",
                                     hjust = 0),
    plot.margin       = margin(10, 15, 10, 10)
  ) +

  labs(
    title    = "Figure 1. Recovery from Diaphragm Paralysis After Congenital Heart Surgery",
    subtitle = paste0(
      "Turnbull nonparametric maximum likelihood estimate with 2000-iteration ",
      "bootstrap 95% CI;\n",
      "overlaid with Weibull, log-normal, and log-logistic accelerated failure ",
      "time model fits"
    ),
    caption  = paste0(
      "Shaded band: bootstrap 95% pointwise confidence interval (2000 resamples). ",
      "Dashed lines indicate selected model (lowest AIC).\n",
      "All models fitted using interval-censored likelihood. ",
      "AIC = Akaike information criterion."
    )
  )

print(fig1)


# ============================================================
# 6. 出力
# ============================================================
output_dir <- "~/Desktop"

ggsave(
  filename = file.path(output_dir, "Figure1_DP_Recovery_Turnbull_AFT.pdf"),
  plot     = fig1,
  width    = 7.0,
  height   = 5.5,
  device   = cairo_pdf
)
ggsave(
  filename    = file.path(output_dir, "Figure1_DP_Recovery_Turnbull_AFT.tiff"),
  plot        = fig1,
  width       = 7.0,
  height      = 5.5,
  dpi         = 600,
  compression = "lzw"
)

cat("\n\u2713 Figure1_DP_Recovery_Turnbull_AFT.pdf / .tiff を出力しました\n")
cat("出力先:", output_dir, "\n\n")

# ============================================================
# 7. 結果サマリー（Results記述用）
# ============================================================
cat("========================================\n")
cat("RESULTS SUMMARY\n")
cat("========================================\n")
cat(sprintf("N                        : %d\n", nrow(df)))
cat(sprintf("Recovered                : %d (%.1f%%)\n",
            sum(!is.na(df$right_time)),
            sum(!is.na(df$right_time)) / nrow(df) * 100))
cat(sprintf("Right-censored           : %d (%.1f%%)\n",
            sum(is.na(df$right_time)),
            sum(is.na(df$right_time)) / nrow(df) * 100))
cat(sprintf("Median recovery          : %d days (95%% CI: %d\u2013%d)\n",
            med_val, med_lo, med_hi))
cat(sprintf("Best-fit AFT model (AIC) : %s\n", best_model))
cat("Model AIC/BIC comparison:\n")
print(aic_tbl)
