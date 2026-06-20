suppressPackageStartupMessages({
  library(Seurat)
  library(slingshot)
  library(ggplot2)
})

outdir <- path.expand(
  "~/mtec_analysis/controls/analysis_outs"
)

load(file.path(outdir, "controls_merged_slingshot_seurat.rda"))
load(file.path(outdir, "controls_merged_slingshot.rda"))

umap <- as.data.frame(
  mtec_wt_slingshot@dr$umap@cell.embeddings
)

colnames(umap)[1:2] <- c("UMAP1", "UMAP2")

umap$stage <- as.character(
  mtec_wt_slingshot@meta.data[
    rownames(umap),
    "stage"
  ]
)

# 使用论文中的名称
umap$stage[umap$stage == "Early_Aire"] <- "TAC_TEC"
umap$stage[umap$stage == "Cortico_medullary"] <- "cTEC"

umap$stage <- factor(
  umap$stage,
  levels = c(
    "cTEC",
    "Ccl21a_high",
    "TAC_TEC",
    "Aire_positive",
    "Late_Aire",
    "Tuft"
  )
)

stage_cols <- c(
  "cTEC" = "#CC6600",
  "Ccl21a_high" = "#009933",
  "TAC_TEC" = "#0066CC",
  "Aire_positive" = "#660099",
  "Late_Aire" = "#FF0000",
  "Tuft" = "#990000"
)

# Figure 2D 只使用第一条 lineage
curve <- slingCurves(sce.umap)[[1]]

curve_df <- as.data.frame(
  curve$s[curve$ord, 1:2]
)

colnames(curve_df) <- c("UMAP1", "UMAP2")

p2d <- ggplot(
  umap,
  aes(UMAP1, UMAP2, color = stage)
) +
  geom_point(size = 0.65, alpha = 0.95) +
  geom_path(
    data = curve_df,
    aes(UMAP1, UMAP2),
    inherit.aes = FALSE,
    linewidth = 1.1,
    color = "black"
  ) +
  scale_color_manual(
    values = stage_cols,
    drop = FALSE
  ) +
  coord_fixed() +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    color = NULL
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "right"
  )

ggsave(
  file.path(outdir, "Figure2D_reproduced.pdf"),
  p2d,
  width = 7,
  height = 6
)

ggsave(
  file.path(outdir, "Figure2D_reproduced.png"),
  p2d,
  width = 7,
  height = 6,
  dpi = 300
)

cat("Figure 2D completed.\n")
