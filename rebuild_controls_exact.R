suppressPackageStartupMessages({
  library(Seurat)
  library(mTEC.10x.pipeline)
})

root <- path.expand("~/mtec_analysis")
project_name <- "RankL_ablation"

make_empty_object <- function(sample_name) {
  matrix_dir <- file.path(
    root,
    sample_name,
    paste0(sample_name, "_count"),
    "outs",
    "filtered_gene_bc_matrices",
    "mm10"
  )

  if (!dir.exists(matrix_dir)) {
    stop("Cannot find matrix directory: ", matrix_dir)
  }

  output_dir <- file.path(root, sample_name, "analysis_outs")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  output_file <- file.path(
    output_dir,
    paste0("seurat_", sample_name, "_empty.rda")
  )

  cat("\nCreating Seurat object for", sample_name, "\n")
  cat("Input:", matrix_dir, "\n")

  mtec_data <- Read10X(data.dir = matrix_dir)

  mtec <- CreateSeuratObject(
    raw.data = mtec_data,
    min.cells = 3,
    min.genes = 200,
    project = project_name
  )

  cat(
    "Initial object:",
    nrow(mtec@raw.data), "genes x",
    ncol(mtec@raw.data), "cells\n"
  )

  save(mtec, file = output_file)
  cat("Saved:", output_file, "\n")

  return(output_file)
}

load_and_process <- function(file_path, sample_name) {
  object_name <- load(file_path)
  seurat_obj <- get(object_name[1])

  cat("\nProcessing", sample_name, "\n")

  seurat_obj <- add_perc_mito(seurat_obj)
  seurat_obj <- process_cells(seurat_obj)
  seurat_obj@meta.data$exp <- sample_name

  cat(
    "After single-sample processing:",
    nrow(seurat_obj@data), "genes x",
    ncol(seurat_obj@data), "cells\n"
  )

  return(seurat_obj)
}

# Create the two single-sample Seurat objects
beg_file <- make_empty_object("isoControlBeg")
end_file <- make_empty_object("isoControlEnd")

beg_obj <- load_and_process(beg_file, "isoControlBeg")
end_obj <- load_and_process(end_file, "isoControlEnd")

# Merge controls exactly as in the authors' control script
cat("\nMerging the two controls\n")

mtec_wt <- Seurat::MergeSeurat(
  object1 = beg_obj,
  object2 = end_obj,
  add.cell.id1 = "isoControlBeg",
  add.cell.id2 = "isoControlEnd"
)

controls_out <- file.path(root, "controls", "analysis_outs")
dir.create(controls_out, recursive = TRUE, showWarnings = FALSE)

quality_pdf <- file.path(
  controls_out,
  "controls_merged_quality.pdf"
)

pdf(quality_pdf)

mtec_wt <- add_perc_mito(mtec_wt)
qc_plot(mtec_wt)

cat("\nProcessing merged control object\n")
mtec_wt <- process_cells(mtec_wt)

cat("\nRunning PCA significance plots\n")
PC_plots(
  mtec_wt,
  jackstraw = TRUE,
  test_pcs = 1:20
)

cat("\nRunning authors' clustering with PCs 1:13\n")
mtec_wt <- group_cells(
  mtec_wt,
  dims_use = 1:13
)


# 先保存已经完成分群的对象，避免后续画图失败导致结果丢失
output_file <- file.path(
  controls_out,
  "controls_merged_empty.rda"
)

save(mtec_wt, file = output_file)
cat("Saved clustered object:", output_file, "
")

# 不再调用作者包中参数名写错的 tSNE_PCA()
# 直接调用 plotDimRed()，并明确传入 seurat_object
try(
  print(
    plotDimRed(
      seurat_object = mtec_wt,
      col_by = "exp",
      plot_type = "umap"
    )
  ),
  silent = TRUE
)

try(
  print(
    plotDimRed(
      seurat_object = mtec_wt,
      col_by = "cluster",
      plot_type = "umap"
    )
  ),
  silent = TRUE
)

try(
  print(
    plotDimRed(
      seurat_object = mtec_wt,
      col_by = "cluster",
      plot_type = "pca"
    )
  ),
  silent = TRUE
)

try(
  print(
    plotDimRed(
      seurat_object = mtec_wt,
      col_by = "cluster",
      plot_type = "tsne"
    )
  ),
  silent = TRUE
)

dev.off()


cat("\n====================================\n")
cat("CONTROL CLUSTERING COMPLETED\n")
cat("Cells:", ncol(mtec_wt@data), "\n")
cat("Genes:", nrow(mtec_wt@data), "\n")
cat("\nCluster counts:\n")
print(table(mtec_wt@ident))
cat("\nSaved:", output_file, "\n")
cat("QC PDF:", quality_pdf, "\n")
cat("====================================\n")
