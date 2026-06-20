suppressPackageStartupMessages({
  library(Seurat)
  library(mTEC.10x.pipeline)
  library(slingshot)
  library(RColorBrewer)
})

input_file <- path.expand(
  "~/mtec_analysis/controls/analysis_outs/controls_merged_empty.rda"
)

output_dir <- path.expand(
  "~/mtec_analysis/controls/analysis_outs"
)

loaded_name <- load(input_file)
mtec_wt <- get(loaded_name[1])

cat("Loaded object:", loaded_name[1], "\n")
cat("Cells:", ncol(mtec_wt@data), "\n")
cat("Genes:", nrow(mtec_wt@data), "\n")

# 作者原始 cluster → stage 映射
stage_list <- c(
  "0" = "Aire_positive",
  "1" = "Late_Aire",
  "2" = "Early_Aire",
  "3" = "Aire_positive",
  "4" = "Ccl21a_high",
  "5" = "Cortico_medullary",
  "6" = "Late_Aire",
  "7" = "Tuft",
  "8" = "unknown",
  "9" = "unknown"
)

clusters <- as.character(mtec_wt@ident)
mtec_wt@meta.data$stage <- unname(stage_list[clusters])

if (any(is.na(mtec_wt@meta.data$stage))) {
  stop("Some clusters could not be mapped to stages.")
}

stage_levels <- c(
  "Cortico_medullary",
  "Ccl21a_high",
  "Early_Aire",
  "Aire_positive",
  "Late_Aire",
  "Tuft",
  "unknown"
)

mtec_wt@meta.data$stage <- factor(
  mtec_wt@meta.data$stage,
  levels = stage_levels
)

cat("\nStage counts before removing unknown:\n")
print(table(mtec_wt@meta.data$stage))

# 保存作者后续作图需要的带 stage 标签对象
seurat_output <- file.path(
  output_dir,
  "seurat_controls_merged.rda"
)

save(mtec_wt, file = seurat_output)

# Slingshot 使用副本，删除 unknown
mtec_wt_slingshot <- mtec_wt
mtec_wt_slingshot <- Seurat::SetAllIdent(
  mtec_wt_slingshot,
  id = "stage"
)

mtec_wt_slingshot <- Seurat::SubsetData(
  mtec_wt_slingshot,
  ident.remove = "unknown"
)

slingshot_umap <- mtec_wt_slingshot@dr$umap@cell.embeddings
cluster_labels <- mtec_wt_slingshot@meta.data[
  rownames(slingshot_umap),
  "stage"
]

cat("\nCells used by Slingshot:",
    nrow(slingshot_umap), "\n")

set.seed(1)

sce.umap <- slingshot(
  slingshot_umap,
  clusterLabels = cluster_labels,
  start.clus = "Early_Aire",
  extend = "n"
)

slingshot_output <- file.path(
  output_dir,
  "controls_merged_slingshot.rda"
)

save(sce.umap, file = slingshot_output)

# 额外保存带 stage 的 Slingshot 输入对象
save(
  mtec_wt_slingshot,
  file = file.path(
    output_dir,
    "controls_merged_slingshot_seurat.rda"
  )
)

# 导出 pseudotime 数值
pseudotime <- as.data.frame(
  slingPseudotime(sce.umap)
)

write.csv(
  pseudotime,
  file.path(output_dir, "controls_pseudotime.csv"),
  quote = FALSE
)

# 诊断图
stage_color <- c(
  "Cortico_medullary" = "#CC6600",
  "Ccl21a_high" = "#009933",
  "Early_Aire" = "#0066CC",
  "Aire_positive" = "#660099",
  "Late_Aire" = "#FF0000",
  "Tuft" = "#990000"
)

pdf(
  file.path(output_dir, "controls_slingshot_check.pdf"),
  width = 7,
  height = 6
)

plot(
  slingshot_umap,
  col = stage_color[as.character(cluster_labels)],
  pch = 16,
  cex = 0.55,
  asp = 1,
  xlab = "UMAP1",
  ylab = "UMAP2"
)

lines(
  sce.umap,
  lwd = 2,
  type = "curves"
)

legend(
  "topright",
  legend = names(stage_color),
  col = stage_color,
  pch = 16,
  cex = 0.7,
  bty = "n"
)

dev.off()

cat("\n============================\n")
cat("SLINGSHOT COMPLETED\n")
cat("Seurat object:", seurat_output, "\n")
cat("Slingshot object:", slingshot_output, "\n")
cat("Pseudotime CSV:",
    file.path(output_dir, "controls_pseudotime.csv"), "\n")
cat("============================\n")
