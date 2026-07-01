# =============================================================
# Supplemental_Figure1_FlowDiagram.R
# Study flow diagram — DP recovery after congenital heart surgery
# Output: ~/Desktop/figures/SupplFigure1_FlowDiagram.pdf / .tiff
# =============================================================

library(ggplot2)
library(grid)

OUT_DIR <- file.path(path.expand("~"), "Desktop", "figures")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

W_IN <- 7.0
H_IN <- 5.5
FONT <- "sans"

draw_flow <- function() {

  p <- ggplot() +
    coord_fixed(xlim = c(0, 10), ylim = c(0, 10), expand = FALSE) +
    theme_void() +
    theme(
      text             = element_text(family = FONT),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )

  cx      <- 5.0
  sz_main <- 3.2
  sz_sub  <- 2.9

  draw_box <- function(xmin, xmax, ymin, ymax, lwd = 0.6)
    annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
             fill = "white", color = "black", linewidth = lwd)

  p <- p +

    # ── Box 1: 総手術数 ──────────────────────────────
    draw_box(2.0, 8.0, 8.2, 9.6) +
    annotate("text", x = cx, y = 9.15,
      label = "Consecutive patients undergoing\ncongenital heart surgery (2009\u20132024)",
      family = FONT, size = sz_sub, hjust = 0.5, lineheight = 1.2) +
    annotate("text", x = cx, y = 8.55, label = "N = 2,147",
      family = FONT, size = sz_main, fontface = "bold", hjust = 0.5) +

    # Arrow 1 -> 2
    annotate("segment", x = cx, xend = cx, y = 8.2, yend = 7.2,
      arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
      linewidth = 0.5) +

    # ── Box 2: DP診断例 ───────────────────────────────
    draw_box(2.0, 8.0, 5.8, 7.2) +
    annotate("text", x = cx, y = 6.75,
      label = "Diagnosed with diaphragm paralysis",
      family = FONT, size = sz_sub, hjust = 0.5) +
    annotate("text", x = cx, y = 6.15, label = "N = 111",
      family = FONT, size = sz_main, fontface = "bold", hjust = 0.5) +

    # 縦矢印 2 -> 3
    annotate("segment", x = cx, xend = cx, y = 5.8, yend = 4.8,
      arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
      linewidth = 0.5) +

    # 横矢印 -> 除外ボックス
    annotate("segment", x = cx, xend = 8.15, y = 5.3, yend = 5.3,
      arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
      linewidth = 0.4) +

    # ── 除外ボックス（主流ボックスと重ならない位置）──
    draw_box(8.2, 9.9, 4.65, 5.95, lwd = 0.5) +
    annotate("text", x = 9.05, y = 5.60,
      label = "Excluded: n = 1",
      family = FONT, size = sz_sub - 0.3, fontface = "bold", hjust = 0.5) +
    annotate("text", x = 9.05, y = 5.10,
      label = "Insufficient\nimaging follow-up",
      family = FONT, size = sz_sub - 0.4, hjust = 0.5, lineheight = 1.15) +

    # ── Box 3: 解析対象 ───────────────────────────────
    draw_box(2.0, 8.0, 3.4, 4.8) +
    annotate("text", x = cx, y = 4.35,
      label = "Included in primary analysis",
      family = FONT, size = sz_sub, hjust = 0.5) +
    annotate("text", x = cx, y = 3.75, label = "N = 110",
      family = FONT, size = sz_main, fontface = "bold", hjust = 0.5) +

    # Arrow 3 -> 4
    annotate("segment", x = cx, xend = cx, y = 3.4, yend = 2.4,
      arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
      linewidth = 0.5) +

    # ── Box 4: アウトカム ─────────────────────────────
    draw_box(2.0, 8.0, 1.0, 2.4) +
    annotate("text", x = cx, y = 1.95,
      label = "Recovery confirmed: 84 (76.4%)",
      family = FONT, size = sz_sub, hjust = 0.5) +
    annotate("text", x = cx, y = 1.35,
      label = "Right-censored: 26 (23.6%)",
      family = FONT, size = sz_sub, hjust = 0.5)

  p
}

fig <- draw_flow()

# ── PDF出力 ──────────────────────────────────────────────────
pdf(
  file   = file.path(OUT_DIR, "SupplFigure1_FlowDiagram.pdf"),
  width  = W_IN,
  height = H_IN
)
print(fig)
dev.off()

# ── TIFF出力（600 dpi, LZW）─────────────────────────────────
tiff(
  filename    = file.path(OUT_DIR, "SupplFigure1_FlowDiagram.tiff"),
  width       = W_IN,
  height      = H_IN,
  units       = "in",
  res         = 600,
  compression = "lzw"
)
print(fig)
dev.off()

cat("\u2713 SupplFigure1_FlowDiagram.pdf / .tiff \u3092\u51fa\u529b\u3057\u307e\u3057\u305f\n")
cat("  \u51fa\u529b\u5148:", OUT_DIR, "\n")    annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
             fill = "white", color = "black", linewidth = lwd)

  p <- p +

    # ── Box 1: 総手術数 ──────────────────────────────
    draw_box(2.0, 8.0, 8.2, 9.6) +
    annotate("text", x = cx, y = 9.15,
      label = "Consecutive patients undergoing\ncongenital heart surgery (2009\u20132024)",
      family = FONT, size = sz_sub, hjust = 0.5, lineheight = 1.2) +
    annotate("text", x = cx, y = 8.55, label = "N = 2,147",
      family = FONT, size = sz_main, fontface = "bold", hjust = 0.5) +

    # Arrow 1 -> 2
    annotate("segment", x = cx, xend = cx, y = 8.2, yend = 7.2,
      arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
      linewidth = 0.5) +

    # ── Box 2: DP診断例 ───────────────────────────────
    draw_box(2.0, 8.0, 5.8, 7.2) +
    annotate("text", x = cx, y = 6.75,
      label = "Diagnosed with diaphragm paralysis",
      family = FONT, size = sz_sub, hjust = 0.5) +
    annotate("text", x = cx, y = 6.15, label = "N = 111",
      family = FONT, size = sz_main, fontface = "bold", hjust = 0.5) +

    # 縦矢印 2 -> 3
    annotate("segment", x = cx, xend = cx, y = 5.8, yend = 4.8,
      arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
      linewidth = 0.5) +

    # 横矢印 -> 除外ボックス
    annotate("segment", x = cx, xend = 8.15, y = 5.3, yend = 5.3,
      arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
      linewidth = 0.4) +

    # ── 除外ボックス（主流ボックスと重ならない位置）──
    draw_box(8.2, 9.9, 4.65, 5.95, lwd = 0.5) +
    annotate("text", x = 9.05, y = 5.60,
      label = "Excluded: n = 1",
      family = FONT, size = sz_sub - 0.3, fontface = "bold", hjust = 0.5) +
    annotate("text", x = 9.05, y = 5.10,
      label = "Insufficient\nimaging follow-up",
      family = FONT, size = sz_sub - 0.4, hjust = 0.5, lineheight = 1.15) +

    # ── Box 3: 解析対象 ───────────────────────────────
    draw_box(2.0, 8.0, 3.4, 4.8) +
    annotate("text", x = cx, y = 4.35,
      label = "Included in primary analysis",
      family = FONT, size = sz_sub, hjust = 0.5) +
    annotate("text", x = cx, y = 3.75, label = "N = 110",
      family = FONT, size = sz_main, fontface = "bold", hjust = 0.5) +

    # Arrow 3 -> 4
    annotate("segment", x = cx, xend = cx, y = 3.4, yend = 2.4,
      arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
      linewidth = 0.5) +

    # ── Box 4: アウトカム ─────────────────────────────
    draw_box(2.0, 8.0, 1.0, 2.4) +
    annotate("text", x = cx, y = 1.95,
      label = "Recovery confirmed: 84 (76.4%)",
      family = FONT, size = sz_sub, hjust = 0.5) +
    annotate("text", x = cx, y = 1.35,
      label = "Right-censored: 26 (23.6%)",
      family = FONT, size = sz_sub, hjust = 0.5)

  p
}

fig <- draw_flow()

# ── PDF出力 ──────────────────────────────────────────────────
pdf(
  file   = file.path(OUT_DIR, "SupplFigure1_FlowDiagram.pdf"),
  width  = W_IN,
  height = H_IN
)
print(fig)
dev.off()

# ── TIFF出力（600 dpi, LZW）─────────────────────────────────
tiff(
  filename    = file.path(OUT_DIR, "SupplFigure1_FlowDiagram.tiff"),
  width       = W_IN,
  height      = H_IN,
  units       = "in",
  res         = 600,
  compression = "lzw"
)
print(fig)
dev.off()

cat("\u2713 SupplFigure1_FlowDiagram.pdf / .tiff \u3092\u51fa\u529b\u3057\u307e\u3057\u305f\n")
cat("  \u51fa\u529b\u5148:", OUT_DIR, "\n")
