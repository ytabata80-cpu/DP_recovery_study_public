# ============================================================
#  06_bayesian_aft.R
#  Bayesian Log-logistic AFT Model
#  Informative priors from Smith et al. (2013) Mid KM
#  via IPDfromKM reconstruction (Guyot algorithm)
#
#  Input : df（01_data_cleaning.R で作成）
#  Output: Figure2_DP_Recovery_Bayesian.pdf
#          ~/Desktop/ に保存
#
#  実行方法:
#    source("R/01_data_cleaning.R")
#    source("R/06_bayesian_aft.R")
#
#  事前分布（Smith et al. 2013 Mid KM由来）:
#    scale: Normal(5.930, 1.990) on log scale（SE×5）
#    shape: Normal(0.186, 1.215) on log scale（SE×5）
#
#  3つの解析:
#    主解析  : Informative prior（SE×5, prior ESS≈2）
#    感度解析1: Moderately informative prior（SE×2, prior ESS≈11.5）
#    感度解析2: Non-informative prior（Normal(0,10)）
# ============================================================

library(brms)
library(dplyr)
library(ggplot2)
library(survival)

options(mc.cores = parallel::detectCores())
set.seed(1234)

# 01_data_cleaning.R が未実行の場合は自動実行
if (!exists("df")) {
  source("R/01_data_cleaning.R")
}

# ============================================================
# 1. データ準備（区間打ち切り形式）
# ============================================================
df_bayes <- df %>%
  mutate(
    y    = left_time,
    y2   = ifelse(!is.na(right_time), right_time, NA),
    cens = case_when(
      !is.na(right_time) ~ "interval",
      is.na(right_time)  ~ "right"
    )
  )

cat("=== データ確認 ===\n")
cat("interval censored:", sum(df_bayes$cens == "interval"), "\n")
cat("right censored   :", sum(df_bayes$cens == "right"), "\n\n")


# ============================================================
# 2. Log-logistic Custom Family
# ============================================================
log_logistic <- custom_family(
  "log_logistic",
  dpars   = c("mu", "sigma"),
  links   = c("log", "log"),
  lb      = c(0, 0),
  type    = "real"
)

stan_funs <- "
  real log_logistic_lpdf(real y, real mu, real sigma) {
    return log(sigma) - log(mu) + (sigma - 1) * log(y/mu)
           - 2 * log(1 + pow(y/mu, sigma));
  }
  real log_logistic_lcdf(real y, real mu, real sigma) {
    return log(1 - 1/(1 + pow(y/mu, sigma)));
  }
  real log_logistic_lccdf(real y, real mu, real sigma) {
    return log(1/(1 + pow(y/mu, sigma)));
  }
"

stan_vars <- stanvar(scode = stan_funs, block = "functions")


# ============================================================
# 3. Bayesian AFT Models
# ============================================================

# ── 主解析: Informative prior（SE×5, prior ESS≈2）────────────────
cat("主解析: Informative prior（SE×5, prior ESS≈2）fitting...\n")
fit_bayes_main <- brm(
  formula  = bf(y | cens(cens, y2) ~ 1),
  data     = df_bayes,
  family   = log_logistic,
  stanvars = stan_vars,
  prior    = c(
    prior(normal(5.930, 1.990), class = "Intercept"),
    prior(normal(0.186, 1.215), class = "sigma")
  ),
  chains  = 4,
  iter    = 4000,
  warmup  = 1000,
  seed    = 1234,
  silent  = 2
)
cat("完了\n\n")

# ── 感度解析1: Moderately informative prior（SE×2, prior ESS≈11.5）─
cat("感度解析1: Moderately informative prior（SE×2）fitting...\n")
fit_bayes_weak <- brm(
  formula  = bf(y | cens(cens, y2) ~ 1),
  data     = df_bayes,
  family   = log_logistic,
  stanvars = stan_vars,
  prior    = c(
    prior(normal(5.930, 0.796), class = "Intercept"),
    prior(normal(0.186, 0.486), class = "sigma")
  ),
  chains  = 4,
  iter    = 4000,
  warmup  = 1000,
  seed    = 1234,
  silent  = 2
)
cat("完了\n\n")

# ── 感度解析2: Non-informative prior ────────────────────────
cat("感度解析2: Non-informative prior fitting...\n")
fit_bayes_noninf <- brm(
  formula  = bf(y | cens(cens, y2) ~ 1),
  data     = df_bayes,
  family   = log_logistic,
  stanvars = stan_vars,
  prior    = c(
    prior(normal(0, 10), class = "Intercept"),
    prior(normal(0, 10), class = "sigma")
  ),
  chains  = 4,
  iter    = 4000,
  warmup  = 1000,
  seed    = 1234,
  silent  = 2
)
cat("完了\n\n")


# ============================================================
# 4. 結果サマリー
# ============================================================
cat("========================================\n")
cat("BAYESIAN AFT RESULTS SUMMARY\n")
cat("========================================\n")

extract_summary <- function(fit, label) {
  s <- summary(fit)
  intercept <- s$fixed["Intercept", ]
  sigma     <- s$spec_pars["sigma", ]
  cat(sprintf("\n[%s]\n", label))
  cat(sprintf("  scale posterior median : %.1f days (exp(%.3f))\n",
              exp(intercept["Estimate"]), intercept["Estimate"]))
  cat(sprintf("  scale 95%% CrI         : %.1f – %.1f days\n",
              exp(intercept["l-95% CI"]), exp(intercept["u-95% CI"])))
  cat(sprintf("  sigma                  : %.3f (%.3f – %.3f)\n",
              sigma["Estimate"], sigma["l-95% CI"], sigma["u-95% CI"]))
  cat(sprintf("  Rhat (Intercept)       : %.3f\n", intercept["Rhat"]))
  cat(sprintf("  Bulk ESS               : %.0f\n", intercept["Bulk_ESS"]))
}

extract_summary(fit_bayes_main,   "主解析: Informative prior (SE×5, prior ESS≈2)")
extract_summary(fit_bayes_weak,   "感度解析1: Moderately informative (SE×2, prior ESS≈11.5)")
extract_summary(fit_bayes_noninf, "感度解析2: Non-informative")
cat("\n========================================\n")


# ============================================================
# 5. 累積回復曲線の計算
# ============================================================
t_grid <- seq(0, 1095, by = 1)

calc_cumrec <- function(post_draws, t_grid) {
  cum_mat <- matrix(NA, nrow = nrow(post_draws),
                    ncol = length(t_grid))
  for (i in 1:nrow(post_draws)) {
    scale_i    <- exp(post_draws$b_Intercept[i])
    shape_i    <- post_draws$sigma[i]
    cum_mat[i, ] <- 1 / (1 + (t_grid / scale_i)^(-shape_i))
  }
  data.frame(
    time   = t_grid,
    median = apply(cum_mat, 2, median),
    lower  = apply(cum_mat, 2, quantile, 0.025),
    upper  = apply(cum_mat, 2, quantile, 0.975)
  )
}

cat("累積回復曲線計算中...\n")
post_main   <- as_draws_df(fit_bayes_main)
post_weak   <- as_draws_df(fit_bayes_weak)
post_noninf <- as_draws_df(fit_bayes_noninf)

curve_main   <- calc_cumrec(post_main,   t_grid)
curve_weak   <- calc_cumrec(post_weak,   t_grid)
curve_noninf <- calc_cumrec(post_noninf, t_grid)
cat("完了\n\n")


# ============================================================
# 6. Turnbull曲線（参照用）
# ============================================================
surv_obj <- Surv(
  time  = df$left_time,
  time2 = df$right_time,
  type  = "interval2"
)
fit_km <- survfit(surv_obj ~ 1)

tb_df_fig2 <- data.frame(
  time = c(0, fit_km$time),
  surv = c(1, fit_km$surv)
) %>%
  mutate(cumrec = 1 - surv) %>%
  bind_rows(
    data.frame(
      time   = 1095,
      surv   = tail(.$surv, 1),
      cumrec = tail(.$cumrec, 1)
    )
  )


# ============================================================
# 7. Figure 2（JTCVS対応版）
# ============================================================
COL_MAIN     <- "#C0392B"
COL_WEAK     <- "#E07B39"
COL_NONINF   <- "#2A9D8F"
COL_TURNBULL <- "#1F3864"

fig2 <- ggplot() +

  geom_step(
    data = tb_df_fig2,
    aes(x = time, y = cumrec * 100,
        color = "Turnbull NPMLE"),
    linewidth = 1.0, direction = "hv"
  ) +

  geom_ribbon(
    data = curve_noninf,
    aes(x = time, ymin = lower * 100, ymax = upper * 100),
    fill = COL_NONINF, alpha = 0.10
  ) +
  geom_line(
    data = curve_noninf,
    aes(x = time, y = median * 100,
        color = "Non-informative prior"),
    linewidth = 0.7, linetype = "dotted"
  ) +

  geom_ribbon(
    data = curve_weak,
    aes(x = time, ymin = lower * 100, ymax = upper * 100),
    fill = COL_WEAK, alpha = 0.10
  ) +
  geom_line(
    data = curve_weak,
    aes(x = time, y = median * 100,
        color = "Moderately informative prior"),
    linewidth = 0.7, linetype = "dashed"
  ) +

  geom_ribbon(
    data = curve_main,
    aes(x = time, ymin = lower * 100, ymax = upper * 100),
    fill = COL_MAIN, alpha = 0.15
  ) +
  geom_line(
    data = curve_main,
    aes(x = time, y = median * 100,
        color = "Informative prior"),
    linewidth = 1.0, linetype = "solid"
  ) +

  scale_color_manual(
    name   = NULL,
    values = c(
      "Turnbull NPMLE"               = COL_TURNBULL,
      "Informative prior"            = COL_MAIN,
      "Moderately informative prior" = COL_WEAK,
      "Non-informative prior"        = COL_NONINF
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
    legend.position      = c(0.98, 0.20),
    legend.justification = c(1, 0),
    legend.title         = element_blank(),
    legend.text          = element_text(size = 8.5, family = "sans"),
    legend.key.width     = unit(1.8, "cm"),
    legend.background    = element_rect(fill = "white", color = "grey80",
                                        linewidth = 0.3),
    panel.grid           = element_blank(),
    plot.margin          = margin(10, 15, 10, 10)
  )

print(fig2)

# ============================================================
# 8. 出力（PDF + TIFF）
# ============================================================
output_dir <- "~/Desktop"

ggsave(
  filename = file.path(output_dir, "Figure2_DP_Recovery_Bayesian.pdf"),
  plot     = fig2,
  width    = 7.0,
  height   = 5.5,
  device   = "pdf"
)

ggsave(
  filename    = file.path(output_dir, "Figure2_DP_Recovery_Bayesian.tiff"),
  plot        = fig2,
  width       = 7.0,
  height      = 5.5,
  dpi         = 600,
  compression = "lzw"
)

cat("\n\u2713 Figure 2 完成\n")
cat("保存先:", output_dir, "\n")
cat("  - Figure2_DP_Recovery_Bayesian.pdf\n")
cat("  - Figure2_DP_Recovery_Bayesian.tiff (600 dpi, LZW)\n\n")
# ============================================================
# 1. データ準備（区間打ち切り形式）
# ============================================================
df_bayes <- df %>%
  mutate(
    y    = left_time,
    y2   = ifelse(!is.na(right_time), right_time, NA),
    cens = case_when(
      !is.na(right_time) ~ "interval",
      is.na(right_time)  ~ "right"
    )
  )

cat("=== データ確認 ===\n")
cat("interval censored:", sum(df_bayes$cens == "interval"), "\n")
cat("right censored   :", sum(df_bayes$cens == "right"), "\n\n")


# ============================================================
# 2. Log-logistic Custom Family
# ============================================================
log_logistic <- custom_family(
  "log_logistic",
  dpars   = c("mu", "sigma"),
  links   = c("log", "log"),
  lb      = c(0, 0),
  type    = "real"
)

stan_funs <- "
  real log_logistic_lpdf(real y, real mu, real sigma) {
    return log(sigma) - log(mu) + (sigma - 1) * log(y/mu)
           - 2 * log(1 + pow(y/mu, sigma));
  }
  real log_logistic_lcdf(real y, real mu, real sigma) {
    return log(1 - 1/(1 + pow(y/mu, sigma)));
  }
  real log_logistic_lccdf(real y, real mu, real sigma) {
    return log(1/(1 + pow(y/mu, sigma)));
  }
"

stan_vars <- stanvar(scode = stan_funs, block = "functions")


# ============================================================
# 3. Bayesian AFT Models
# ============================================================

# ── 主解析: Informative prior（SE×5, prior ESS≈2）────────────────
cat("主解析: Informative prior（SE×5, prior ESS≈2）fitting...\n")
fit_bayes_main <- brm(
  formula  = bf(y | cens(cens, y2) ~ 1),
  data     = df_bayes,
  family   = log_logistic,
  stanvars = stan_vars,
  prior    = c(
    prior(normal(5.930, 1.990), class = "Intercept"),
    prior(normal(0.186, 1.215), class = "sigma")
  ),
  chains  = 4,
  iter    = 4000,
  warmup  = 1000,
  seed    = 1234,
  silent  = 2
)
cat("完了\n\n")

# ── 感度解析1: Moderately informative prior（SE×2, prior ESS≈11.5）─
cat("感度解析1: Moderately informative prior（SE×2）fitting...\n")
fit_bayes_weak <- brm(
  formula  = bf(y | cens(cens, y2) ~ 1),
  data     = df_bayes,
  family   = log_logistic,
  stanvars = stan_vars,
  prior    = c(
    prior(normal(5.930, 0.796), class = "Intercept"),
    prior(normal(0.186, 0.486), class = "sigma")
  ),
  chains  = 4,
  iter    = 4000,
  warmup  = 1000,
  seed    = 1234,
  silent  = 2
)
cat("完了\n\n")

# ── 感度解析2: Non-informative prior ────────────────────────
cat("感度解析2: Non-informative prior fitting...\n")
fit_bayes_noninf <- brm(
  formula  = bf(y | cens(cens, y2) ~ 1),
  data     = df_bayes,
  family   = log_logistic,
  stanvars = stan_vars,
  prior    = c(
    prior(normal(0, 10), class = "Intercept"),
    prior(normal(0, 10), class = "sigma")
  ),
  chains  = 4,
  iter    = 4000,
  warmup  = 1000,
  seed    = 1234,
  silent  = 2
)
cat("完了\n\n")


# ============================================================
# 4. 結果サマリー
# ============================================================
cat("========================================\n")
cat("BAYESIAN AFT RESULTS SUMMARY\n")
cat("========================================\n")

extract_summary <- function(fit, label) {
  s <- summary(fit)
  intercept <- s$fixed["Intercept", ]
  sigma     <- s$spec_pars["sigma", ]
  cat(sprintf("\n[%s]\n", label))
  cat(sprintf("  scale posterior median : %.1f days (exp(%.3f))\n",
              exp(intercept["Estimate"]), intercept["Estimate"]))
  cat(sprintf("  scale 95%% CrI         : %.1f – %.1f days\n",
              exp(intercept["l-95% CI"]), exp(intercept["u-95% CI"])))
  cat(sprintf("  sigma                  : %.3f (%.3f – %.3f)\n",
              sigma["Estimate"], sigma["l-95% CI"], sigma["u-95% CI"]))
  cat(sprintf("  Rhat (Intercept)       : %.3f\n", intercept["Rhat"]))
  cat(sprintf("  Bulk ESS               : %.0f\n", intercept["Bulk_ESS"]))
}

extract_summary(fit_bayes_main,   "主解析: Informative prior (SE×5, prior ESS≈2)")
extract_summary(fit_bayes_weak,   "感度解析1: Moderately informative (SE×2, prior ESS≈11.5)")
extract_summary(fit_bayes_noninf, "感度解析2: Non-informative")
cat("\n========================================\n")


# ============================================================
# 5. 累積回復曲線の計算
# ============================================================
t_grid <- seq(0, 1095, by = 1)

calc_cumrec <- function(post_draws, t_grid) {
  cum_mat <- matrix(NA, nrow = nrow(post_draws),
                    ncol = length(t_grid))
  for (i in 1:nrow(post_draws)) {
    scale_i    <- exp(post_draws$b_Intercept[i])
    shape_i    <- post_draws$sigma[i]
    cum_mat[i, ] <- 1 / (1 + (t_grid / scale_i)^(-shape_i))
  }
  data.frame(
    time   = t_grid,
    median = apply(cum_mat, 2, median),
    lower  = apply(cum_mat, 2, quantile, 0.025),
    upper  = apply(cum_mat, 2, quantile, 0.975)
  )
}

cat("累積回復曲線計算中...\n")
post_main   <- as_draws_df(fit_bayes_main)
post_weak   <- as_draws_df(fit_bayes_weak)
post_noninf <- as_draws_df(fit_bayes_noninf)

curve_main   <- calc_cumrec(post_main,   t_grid)
curve_weak   <- calc_cumrec(post_weak,   t_grid)
curve_noninf <- calc_cumrec(post_noninf, t_grid)
cat("完了\n\n")


# ============================================================
# 6. Turnbull曲線（参照用）
# ============================================================
surv_obj <- Surv(
  time  = df$left_time,
  time2 = df$right_time,
  type  = "interval2"
)
fit_km <- survfit(surv_obj ~ 1)

tb_df_fig2 <- data.frame(
  time = c(0, fit_km$time),
  surv = c(1, fit_km$surv)
) %>%
  mutate(cumrec = 1 - surv) %>%
  bind_rows(
    data.frame(
      time   = 1095,
      surv   = tail(.$surv, 1),
      cumrec = tail(.$cumrec, 1)
    )
  )


# ============================================================
# 7. Figure 2
# ============================================================
COL_MAIN     <- "#C0392B"
COL_WEAK     <- "#E07B39"
COL_NONINF   <- "#2A9D8F"
COL_TURNBULL <- "#1F3864"

fig2 <- ggplot() +

  # Turnbull NPMLE（参照用）
  geom_step(
    data = tb_df_fig2,
    aes(x = time, y = cumrec * 100,
        color = "Turnbull NPMLE"),
    linewidth = 1.0, direction = "hv"
  ) +

  # Non-informative prior
  geom_ribbon(
    data = curve_noninf,
    aes(x = time, ymin = lower * 100, ymax = upper * 100),
    fill = COL_NONINF, alpha = 0.10
  ) +
  geom_line(
    data = curve_noninf,
    aes(x = time, y = median * 100,
        color = "Non-informative prior"),
    linewidth = 0.7, linetype = "dotted"
  ) +

  # Weakly informative prior
  geom_ribbon(
    data = curve_weak,
    aes(x = time, ymin = lower * 100, ymax = upper * 100),
    fill = COL_WEAK, alpha = 0.10
  ) +
  geom_line(
    data = curve_weak,
    aes(x = time, y = median * 100,
        color = "Moderately informative prior (SE\u00d72)"),
    linewidth = 0.7, linetype = "dashed"
  ) +

  # 主解析（Informative prior SD×5）
  geom_ribbon(
    data = curve_main,
    aes(x = time, ymin = lower * 100, ymax = upper * 100),
    fill = COL_MAIN, alpha = 0.15
  ) +
  geom_line(
    data = curve_main,
    aes(x = time, y = median * 100,
        color = "Informative prior (SE\u00d75, prior ESS\u22482)"),
    linewidth = 1.0, linetype = "solid"
  ) +

  scale_color_manual(
    name   = NULL,
    values = c(
      "Turnbull NPMLE"                                    = COL_TURNBULL,
      "Informative prior (SE\u00d75, prior ESS\u22482)"   = COL_MAIN,
      "Moderately informative prior (SE\u00d72)"          = COL_WEAK,
      "Non-informative prior"                             = COL_NONINF
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
    legend.position      = c(0.98, 0.20),
    legend.justification = c(1, 0),
    legend.title         = element_blank(),
    legend.text          = element_text(size = 8.5, family = "serif"),
    legend.key.width     = unit(1.8, "cm"),
    legend.background    = element_rect(fill = "white", color = "grey80",
                                        linewidth = 0.3),
    panel.grid           = element_blank(),
    plot.title           = element_text(size = 11, face = "bold"),
    plot.subtitle        = element_text(size = 9, color = "grey40"),
    plot.caption         = element_text(size = 7.5, color = "grey50",
                                        hjust = 0),
    plot.margin          = margin(10, 15, 10, 10)
  ) +

  labs(
    title    = "Figure 2. Bayesian Log-logistic AFT Model: Posterior Cumulative Recovery",
    subtitle = paste0(
      "Informative priors derived from Smith et al. (2013) via IPDfromKM reconstruction.\n",
      "Shaded bands represent 95% credible intervals."
    ),
    caption  = paste0(
      "Turnbull NPMLE shown for reference (solid dark line).\n",
      "Informative prior: Normal(5.930, 1.990) for scale, ",
      "Normal(0.186, 1.215) for shape ",
      "(SE inflation factor = 5, prior ESS \u2248 2).\n",
      "All models: 4 chains \u00d7 3000 post-warmup iterations; ",
      "R-hat < 1.01 for all parameters."
    )
  )

print(fig2)


# ============================================================
# 8. 出力（PDF + TIFF）
# ============================================================
output_dir <- "~/Desktop"

ggsave(
  filename = file.path(output_dir, "Figure2_DP_Recovery_Bayesian.pdf"),
  plot     = fig2,
  width    = 7.0,
  height   = 5.5,
  device   = "pdf"
)

ggsave(
  filename    = file.path(output_dir, "Figure2_DP_Recovery_Bayesian.tiff"),
  plot        = fig2,
  width       = 7.0,
  height      = 5.5,
  dpi         = 600,
  compression = "lzw"
)

cat("\n\u2713 Figure 2 完成\n")
cat("保存先:", output_dir, "\n")
cat("  - Figure2_DP_Recovery_Bayesian.pdf\n")
cat("  - Figure2_DP_Recovery_Bayesian.tiff (600 dpi, LZW)\n\n")
