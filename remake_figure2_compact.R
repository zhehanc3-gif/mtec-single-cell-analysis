suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(gridExtra)
})

outdir <- path.expand(
  "~/mtec_analysis/controls/analysis_outs"
)

# 读取之前已经生成的绘图数据
input_file <- file.path(
  outdir,
  "Figure2_pseudotime_plot_data.csv"
)

plot_df <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Figure 2 只保留主 lineage 的三个阶段
plot_df <- plot_df[
  plot_df$stage %in% c(
    "TAC_TEC",
    "Aire_positive",
    "Late_Aire"
  ),
  ,
  drop = FALSE
]

plot_df$stage <- factor(
  plot_df$stage,
  levels = c(
    "TAC_TEC",
    "Aire_positive",
    "Late_Aire"
  )
)

# 作者使用的颜色
stage_colors <- c(
  "TAC_TEC" = "#0066CC",
  "Aire_positive" = "#660099",
  "Late_Aire" = "#FF0000"
)

stage_labels <- c(
  "TAC_TEC" = "TAC-TEC",
  "Aire_positive" = "Aire-positive",
  "Late_Aire" = "Late-Aire"
)

# 与论文接近的纵轴范围
y_limits <- list(
  Aire = c(0, 3.6),
  Fezf2 = c(0, 4.0),
  aire_genes = c(0, 0.105),
  fezf2_genes = c(0, 0.34),
  tra_fantom = c(0, 0.105),
  protein_coding = c(0.05, 0.23)
)

font_family <- "sans"

plot_titles <- list(
  Aire = "Aire",
  Fezf2 = "Fezf2",
  aire_genes = "Aire-dependent genes",
  fezf2_genes = "Fezf2-dependent genes",
  tra_fantom = "TSA genes",
  protein_coding = "Protein coding genes"
)

make_panel <- function(metric, show_legend = FALSE) {

  df <- plot_df
  df$y_value <- df[[metric]]

  ggplot(
    df,
    aes(
      x = pseudotime,
      y = y_value,
      color = stage
    )
  ) +
    # 比旧图更大、更不透明，视觉上会紧凑很多
    geom_point(
      size = 1.15,
      alpha = 0.95
    ) +
    # 与作者函数相同：默认 loess
    geom_smooth(
      aes(group = 1),
      se = FALSE,
      color = "black",
      size = 1
    ) +
    scale_color_manual(
      values = stage_colors,
      breaks = names(stage_colors),
      labels = stage_labels,
      drop = FALSE,
      name = NULL
    ) +
    coord_cartesian(
      xlim = c(0, 16),
      ylim = y_limits[[metric]],
      expand = FALSE
    ) +
    scale_x_continuous(
      breaks = c(0, 5, 10, 15),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    labs(
      title = plot_titles[[metric]],
      x = NULL,
      y = NULL
    ) +
    theme_classic(base_size = 13, base_family = font_family) +
    theme(
      plot.title = element_text(
        size = 17,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 2)
      ),
      axis.text = element_text(
        size = 12,
        color = "black"
      ),
      axis.ticks = element_line(
        size = 0.65,
        color = "black"
      ),
      axis.line = element_line(
        size = 0.8,
        color = "black"
      ),
      plot.margin = margin(
        t = 3,
        r = 7,
        b = 3,
        l = 5
      ),
      legend.position = if (show_legend) "top" else "none",
      legend.text = element_text(size = 14),
      legend.key.width = unit(1.1, "cm")
    )
}

# 按论文中的顺序
p_aire <- make_panel("Aire")
p_fezf2 <- make_panel("Fezf2")

p_aire_dep <- make_panel("aire_genes")
p_fezf2_dep <- make_panel("fezf2_genes")

p_tsa <- make_panel("tra_fantom")
p_protein <- make_panel("protein_coding")

# 单独提取顶部图例
legend_plot <- make_panel(
  "Aire",
  show_legend = TRUE
)

legend_grob <- ggplotGrob(legend_plot)

legend_index <- which(
  vapply(
    legend_grob$grobs,
    function(x) x$name,
    character(1)
  ) == "guide-box"
)

if (length(legend_index) != 1) {
  stop("无法提取图例。")
}

legend_grob <- legend_grob$grobs[[legend_index]]

# 左侧两组公共纵轴标题
y_label_top <- textGrob(
  "Normalized log\nexpression",
  rot = 90,
  gp = gpar(fontsize = 15, fontfamily = font_family)
)

y_label_bottom <- textGrob(
  "Average normalized log\nexpression",
  rot = 90,
  gp = gpar(fontsize = 15, fontfamily = font_family)
)

# 使用 layout_matrix，让下面的纵轴标题跨两行
panel_grid <- arrangeGrob(
  grobs = list(
    y_label_top,
    p_aire,
    p_fezf2,
    y_label_bottom,
    p_aire_dep,
    p_fezf2_dep,
    p_tsa,
    p_protein
  ),
  layout_matrix = rbind(
    c(1, 2, 3),
    c(4, 5, 6),
    c(4, 7, 8)
  ),
  widths = c(0.13, 1, 1),
  heights = c(1, 1, 1)
)

x_label <- textGrob(
  "Pseudotime",
  gp = gpar(fontsize = 16, fontfamily = font_family)
)

panel_letter <- textGrob(
  "e",
  x = unit(0, "npc"),
  just = "left",
  gp = gpar(fontsize = 19, fontface = "bold", fontfamily = font_family)
)

top_row <- arrangeGrob(
  panel_letter,
  legend_grob,
  ncol = 2,
  widths = c(0.08, 0.92)
)

final_figure <- arrangeGrob(
  top_row,
  panel_grid,
  x_label,
  ncol = 1,
  heights = c(0.09, 0.86, 0.05)
)

# 保存 PDF
pdf(
  file.path(
    outdir,
    "Figure2_pseudotime_compact.pdf"
  ),
  width = 8.2,
  height = 8.7,
  useDingbats = FALSE
)

grid.draw(final_figure)
dev.off()

# 保存 PNG
png(
  file.path(
    outdir,
    "Figure2_pseudotime_compact.png"
  ),
  width = 2460,
  height = 2610,
  res = 300
)

grid.draw(final_figure)
dev.off()

cat("\nCOMPACT FIGURE COMPLETED\n")
cat(
  file.path(
    outdir,
    "Figure2_pseudotime_compact.png"
  ),
  "\n"
)
