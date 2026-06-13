# ============================================================
#  01_data_cleaning.R
#  データ読み込み・除外・派生変数計算
#
#  Input : chd_dp_data_analysis.xlsx / sheet: "Data"
#  Output: df（解析対象データフレーム）
#          excluded（除外症例データフレーム）
#
#  他のスクリプトからは以下で呼び出す:
#    source("R/01_data_cleaning.R")
#
#  除外基準:
#   - left_time が NA（last_abnormal_date 未記録 = 画像なし・文章記述のみ）
#
#  Day 0起算の定義:
#   - Time zero = 手術日（surgery_date）
#   - left_time  = last_abnormal_date - surgery_date
#   - right_time = first_normal_date  - surgery_date（回復未確認はNA）
# ============================================================

library(readxl)
library(dplyr)

# ============================================================
# 1. データ読み込み
# ============================================================
# スプレッドシート構造:
#   Row 1: セクション名  ← スキップ
#   Row 2: 変数名        ← 列名として使用
#   Row 3: ラベル        ← スキップ
#   Row 4: コーディングノート ← スキップ
#   Row 5〜: データ

DATA_PATH <- file.path(
  Sys.getenv("HOME"),
  "Desktop",
  "chd_dp_data_analysis.xlsx"
)
# Step 1: Row 2（変数名行）を取得
col_names_row <- read_excel(
  DATA_PATH,
  sheet     = "Data",
  skip      = 1,
  n_max     = 1,
  col_names = FALSE
) %>% unlist() %>% as.character()

col_names_row <- col_names_row[!is.na(col_names_row)]

# Step 2: Row 5以降のデータを読み込み
df_temp <- read_excel(
  DATA_PATH,
  sheet     = "Data",
  skip      = 4,
  col_names = FALSE
)

col_names_use <- col_names_row[1:ncol(df_temp)]
df_raw        <- df_temp
names(df_raw) <- col_names_use

cat(sprintf("読み込み行数: %d\n", nrow(df_raw)))

# 空行除去
df_raw <- df_raw %>% filter(!is.na(patient_id))
cat(sprintf("有効行数: %d\n\n", nrow(df_raw)))


# ============================================================
# 2. 変数型変換
# ============================================================

# DATE列（文字列 "YYYY-MM-DD" → Date）
date_cols <- c("surgery_date", "last_abnormal_date", "first_normal_date",
               "dp_diagnosis_date", "plication_date")

for (col in date_cols) {
  if (col %in% names(df_raw)) {
    df_raw[[col]] <- as.Date(df_raw[[col]], format = "%Y-%m-%d")
  } else {
    cat(sprintf("注意: 列 '%s' が見つかりません\n", col))
  }
}

# 数値列
num_cols <- c("age_days", "weight_kg", "premature",
              "rachs1_category", "prior_cardiac_surgery",
              "chromosomal_anomaly", "syndrome", "noncardiac_anomaly",
              "cpb_time_min", "xclamp_time_min",
              "circulatory_arrest", "circulatory_arrest_min",
              "days_to_dp_diagnosis",
              "plication_performed", "days_to_plication",
              "followup_days", "recovery_confirmed")

for (col in num_cols) {
  if (col %in% names(df_raw)) {
    df_raw[[col]] <- as.numeric(df_raw[[col]])
  }
}


# ============================================================
# 3. 派生変数の計算
# ============================================================
df_raw <- df_raw %>%
  mutate(
    # Turnbull推定の入力（day 0起算）
    left_time  = as.numeric(last_abnormal_date - surgery_date),
    right_time = as.numeric(first_normal_date  - surgery_date),

    # days_to_dp_diagnosis: 
    days_to_dp_diagnosis = as.numeric(dp_diagnosis_date - surgery_date),

    # Assessment interval（回復確認例のみ・右側打ち切りはNA）
    censoring_interval_days = right_time - left_time,

    # Cox感度解析用（midpoint imputation）
    time_event_midpoint = ifelse(
      !is.na(right_time),
      (left_time + right_time) / 2,
      left_time
    ),
    event_cox = as.integer(!is.na(right_time)),

    # 循環停止なし → 0
    circulatory_arrest_min = ifelse(
      circulatory_arrest == 0 & is.na(circulatory_arrest_min),
      0, circulatory_arrest_min
    )
  )


# ============================================================
# 4. 除外
# ============================================================
excluded <- df_raw %>% filter(is.na(left_time))

cat(sprintf("除外症例数: %d (left_time = NA: 画像なし・文章記述のみ)\n",
            nrow(excluded)))
if (nrow(excluded) > 0) {
  cat("除外症例 patient_id:",
      paste(excluded$patient_id, collapse = ", "), "\n")
}

df <- df_raw %>% filter(!is.na(left_time))
cat(sprintf("解析対象症例数: %d\n\n", nrow(df)))


# ============================================================
# 5. データ確認
# ============================================================
cat("=== 解析対象データ確認 ===\n")
cat("right_time NA（右側打ち切り）:", sum(is.na(df$right_time)), "\n")
cat("recovery_confirmed == 1     :",
    sum(df$recovery_confirmed == 1, na.rm = TRUE), "\n")
cat("plication_performed == 1    :",
    sum(df$plication_performed == 1, na.rm = TRUE), "\n")
cat("plication NA                :",
    sum(is.na(df$days_to_plication)), "\n")
cat("left_time range             :",
    min(df$left_time, na.rm = TRUE), "–",
    max(df$left_time, na.rm = TRUE), "days\n")
cat("right_time range (recovered):",
    min(df$right_time, na.rm = TRUE), "–",
    max(df$right_time, na.rm = TRUE), "days\n\n")

cat("✓ 01_data_cleaning.R 完了: df（n =", nrow(df), "）を作成しました\n")
