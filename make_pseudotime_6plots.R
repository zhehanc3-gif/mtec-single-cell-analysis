suppressPackageStartupMessages({
  library(Seurat)
  library(slingshot)
  library(ggplot2)
  library(gridExtra)
})

outdir <- path.expand(
  "~/mtec_analysis/controls/analysis_outs"
)

# 读取带 gene-set score 的 Seurat 对象
load(
  file.path(outdir, "controls_with_gene_scores.rda")
)

# 读取 Slingshot 对象
load(
  file.path(outdir, "controls_merged_slingshot.rda")
)

expr <- mtec_plot@data
meta <- mtec_plot@meta.data

# 确保 metadata 与表达矩阵顺序一致
meta <- meta[
  colnames(expr),
  ,
  drop = FALSE
]

# 提取 curve1 pseudotime
pt <- as.data.frame(
  slingshot::slingPseudotime(sce.umap)
)

if (!"curve1" %in% colnames(pt)) {
  stop(
    "curve1 not found. Available columns: ",
    paste(colnames(pt), collapse = ", ")
  )
}

pt$cell <- rownames(pt)

# 只保留 curve1 上且 pseudotime 在 0-16 的细胞
pt <- pt[
  !is.na(pt$curve1) &
  pt$curve1 >= 0 &
  pt$curve1 <= 16,
  ,
  drop = FALSE
]

cells <- intersect(
  pt$cell,
  colnames(expr)
)

pt <- pt[
  match(cells, pt$cell),
  ,
  drop = FALSE
]

cat("Cells plotted on curve1:", length(cells), "\n")
cat(
  "Pseudotime range:",
  paste(range(pt$curve1), collapse = " to "),
  "\n"
)

# 建立基础绘图数据
plot_df <- data.frame(
  cell = cells,
  pseudotime = pt$curve1,
  stage = meta[cells, "stage"],
  stringsAsFactors = FALSE
)

plot_df$stage <- factor(
  plot_df$stage,
  levels = c(
    "cTEC",
    "Ccl21a_high",
    "TAC_TEC",
    "Aire_positive",
    "Late_Aire",
    "Tuft"
  )
)

# 读取四个 gene-set score
plot_df$protein_coding <- meta[
  cells,
  "protein_coding"
]

plot_df$tra_fantom <- meta[
  cells,
  "tra_fantom"
]

plot_df$aire_genes <- meta[
  cells,
  "aire_genes"
]

plot_df$fezf2_genes <- meta[
  cells,
  "fezf2_genes"
]

# 读取 Aire 和 Fezf2 单基因表达
plot_df$Aire <- as.numeric(
  expr["Aire", cells]
)

plot_df$Fezf2 <- as.numeric(
  expr["Fezf2", cells]
)

# 保存完整绘图数据，方便以后核查
write.csv(
  plot_df,
  file.path(outdir, "Figure2_pseudotime_plot_data.csv"),
  row.names = FALSE,
  quote = FALSE
)

stage_colors <- c(
  "cTEC" = "#CC6600",
  "Ccl21a_high" = "#009933",
  "TAC_TEC" = "#0066CC",
  "Aire_positive" = "#660099",
  "Late_Aire" = "#FF0000",
  "Tuft" = "#990000"
)

plot_info <- list(
  protein_coding = list(
    label = "Protein-coding genes",
    file = "Figure2_pseudotime_I_protein_coding.pdf"
  ),
  tra_fantom = list(
    label = "Tissue-restricted genes",
    file = "Figure2_pseudotime_II_TRA.pdf"
  ),
  Aire = list(
    label = "Aire",
    file = "Figure2_pseudotime_III_Aire.pdf"
  ),
  aire_genes = list(
    label = "Aire-dependent genes",
    file = "Figure2_pseudotime_IV_Aire_dependent.pdf"
  ),
  Fezf2 = list(
    label = "Fezf2",
    file = "Figure2_pseudotime_V_Fezf2.pdf"
  ),
  fezf2_genes = list(
    label = "Fezf2-dependent genes",
    file = "Figure2_pseudotime_VI_Fezf2_dependent.pdf"
  )
)

plots_for_grid <- list()

for (metric in names(plot_info)) {

  current <- plot_df
  current$y_value <- current[[metric]]

  p <- ggplot(
    current,
    aes(
      x = pseudotime,
      y = y_value,
      color = stage
    )
  ) +
    geom_point(
      size = 0.7,
      alpha = 0.65
    ) +
    geom_smooth(
      aes(group = 1),
      se = FALSE,
      color = "black"
    ) +
    scale_color_manual(
      values = stage_colors,
      drop = FALSE
    ) +
    xlim(0, 16) +
    labs(
      x = "Pseudotime",
      y = plot_info[[metric]]$label,
      color = "Stage"
    ) +
    theme_classic(
      base_size = 14
    )

  ggsave(
    filename = file.path(
      outdir,
      plot_info[[metric]]$file
    ),
    plot = p,
    width = 7,
    height = 3
  )

  plots_for_grid[[metric]] <- p +
    theme(
      legend.position = "none"
    )

  cat("Completed:", metric, "\n")
}

# 合并为六联 PDF
pdf(
  file.path(
    outdir,
    "Figure2_pseudotime_6panel.pdf"
  ),
  width = 14,
  height = 9
)

gridExtra::grid.arrange(
  grobs = plots_for_grid,
  ncol = 2,
  nrow = 3
)

dev.off()

# 合并为六联 PNG
png(
  file.path(
    outdir,
    "Figure2_pseudotime_6panel.png"
  ),
  width = 4200,
  height = 2700,
  res = 300
)

gridExtra::grid.arrange(
  grobs = plots_for_grid,
  ncol = 2,
  nrow = 3
)

dev.off()

cat("\n================================\n")
cat("PSEUDOTIME FIGURES COMPLETED\n")
cat(
  "Combined PDF:",
  file.path(outdir, "Figure2_pseudotime_6panel.pdf"),
  "\n"
)
cat(
  "Combined PNG:",
  file.path(outdir, "Figure2_pseudotime_6panel.png"),
  "\n"
)
cat("================================\n")
