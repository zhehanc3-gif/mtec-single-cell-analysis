suppressPackageStartupMessages({
  library(Seurat)
  library(DESeq2)
  library(SummarizedExperiment)
  library(Single.mTEC.Transcriptomes)
  library(GEOquery)
  library(Biobase)
  library(limma)
})

root <- path.expand("~/mtec_analysis")
outdir <- file.path(root, "controls", "analysis_outs")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. 读取作者提供的参考数据
# ============================================================

data("tras", package = "Single.mTEC.Transcriptomes")
data("aireDependentSansom",
     package = "Single.mTEC.Transcriptomes")
data("fantom", package = "Single.mTEC.Transcriptomes")
data("geneNames", package = "Single.mTEC.Transcriptomes")
data("biotypes", package = "Single.mTEC.Transcriptomes")

# 读取已经完成分群的 control Seurat 对象
loaded <- load(
  file.path(outdir, "seurat_controls_merged.rda")
)

mtec_wt <- get(loaded[1])
seurat_genes <- rownames(mtec_wt@data)

cat("Genes present in control Seurat object:",
    length(seurat_genes), "\n")

# geneNames:
# names   = gene symbols
# values  = Ensembl/FANTOM identifiers
mapping_ids <- sub(
  ",.*$",
  "",
  as.character(geneNames)
)

id_to_symbol <- setNames(
  names(geneNames),
  mapping_ids
)

# 同时兼容输入本身已是 symbol 或仍是 Ensembl ID
to_symbols <- function(x) {
  x <- as.character(unlist(x))
  x <- x[!is.na(x) & nzchar(x)]

  direct <- x[x %in% seurat_genes]

  clean_ids <- sub(",.*$", "", x)
  mapped <- unname(id_to_symbol[clean_ids])
  mapped <- mapped[!is.na(mapped) & nzchar(mapped)]

  unique(c(direct, mapped))
}

# ============================================================
# 2. Protein-coding genes
# ============================================================

protein_ids <- names(biotype)[
  as.character(biotype) == "protein_coding"
]

protein_coding <- to_symbols(protein_ids)
protein_coding <- intersect(
  protein_coding,
  seurat_genes
)

cat("protein_coding:", length(protein_coding), "\n")

# ============================================================
# 3. Aire-dependent genes
# ============================================================

aire_genes <- to_symbols(aireDependentSansom)
aire_genes <- intersect(
  aire_genes,
  protein_coding
)

cat("aire_genes:", length(aire_genes), "\n")

# ============================================================
# 4. FANTOM tissue-restricted genes
#    作者标准：
#    normalized expression > 5
#    expressed in fewer than 5 tissues
# ============================================================

fantom_counts <- DESeq2::counts(
  dxdFANTOM,
  normalized = TRUE
)

fantom_tissues <- as.character(
  SummarizedExperiment::colData(dxdFANTOM)$tissue
)

if (any(is.na(fantom_tissues))) {
  stop("FANTOM tissue metadata contains NA values.")
}

# 每个组织内的样本取平均
means_by_tissue <- sapply(
  split(seq_len(ncol(fantom_counts)), fantom_tissues),
  function(index) {
    rowMeans(
      fantom_counts[, index, drop = FALSE]
    )
  }
)

if (is.null(dim(means_by_tissue))) {
  stop("Failed to construct FANTOM tissue matrix.")
}

# FANTOM 行名可能包含逗号，作者只保留第一个 ID
fantom_ids <- vapply(
  strsplit(
    rownames(means_by_tissue),
    ",",
    fixed = TRUE
  ),
  `[`,
  character(1),
  1
)

# 同一 ID 的重复行取平均
id_split <- split(
  seq_len(nrow(means_by_tissue)),
  fantom_ids
)

means_by_id <- t(
  vapply(
    id_split,
    function(index) {
      colMeans(
        means_by_tissue[index, , drop = FALSE]
      )
    },
    numeric(ncol(means_by_tissue))
  )
)

colnames(means_by_id) <- colnames(means_by_tissue)

# 把 FANTOM/Ensembl ID 转成 gene symbol
matched <- match(
  mapping_ids,
  rownames(means_by_id)
)

valid <- !is.na(matched)

means_by_symbol <- means_by_id[
  matched[valid],
  ,
  drop = FALSE
]

rownames(means_by_symbol) <- names(geneNames)[valid]

# 如果多个 ID 对应同一个 symbol，再取一次平均
symbol_split <- split(
  seq_len(nrow(means_by_symbol)),
  rownames(means_by_symbol)
)

means_by_symbol <- t(
  vapply(
    symbol_split,
    function(index) {
      colMeans(
        means_by_symbol[index, , drop = FALSE]
      )
    },
    numeric(ncol(means_by_symbol))
  )
)

colnames(means_by_symbol) <- colnames(means_by_tissue)

number_of_tissues <- rowSums(
  means_by_symbol > 5,
  na.rm = TRUE
)

tra_fantom <- names(
  number_of_tissues[number_of_tissues < 5]
)

tra_fantom <- intersect(
  tra_fantom,
  protein_coding
)

tra_fantom <- intersect(
  tra_fantom,
  seurat_genes
)

cat("FANTOM tissues:", ncol(means_by_symbol), "\n")
cat("tra_fantom:", length(tra_fantom), "\n")

# 保存 FANTOM tissue counts，方便核查
write.csv(
  data.frame(
    gene = names(number_of_tissues),
    number_of_tissues = as.integer(number_of_tissues),
    stringsAsFactors = FALSE
  ),
  file.path(outdir, "FANTOM_tissue_counts.csv"),
  row.names = FALSE,
  quote = FALSE
)

# ============================================================
# 5. Fezf2-dependent genes
#    严格采用作者对 GSE69105 的 GEO2R/limma 规则：
#    G1-G0、logFC < -1、adjusted P < 0.05
# ============================================================

geo_dir <- file.path(root, "GEO_cache")
dir.create(geo_dir, recursive = TRUE, showWarnings = FALSE)

cat("\nDownloading/loading GSE69105...\n")

gset_list <- GEOquery::getGEO(
  "GSE69105",
  GSEMatrix = TRUE,
  AnnotGPL = TRUE,
  destdir = geo_dir
)

if (length(gset_list) > 1) {
  gset_names <- attr(gset_list, "names")
  idx <- grep("GPL1261", gset_names)

  if (length(idx) == 0) {
    idx <- 1
  } else {
    idx <- idx[1]
  }
} else {
  idx <- 1
}

gset <- gset_list[[idx]]

Biobase::fvarLabels(gset) <- make.names(
  Biobase::fvarLabels(gset)
)

# 作者原始分组：前 5 个 G0，后 5 个 G1
gsms <- "0000011111"

if (ncol(gset) != nchar(gsms)) {
  stop(
    "GSE69105 sample number is ",
    ncol(gset),
    ", but the authors' grouping expects 10 samples."
  )
}

sml <- substring(
  gsms,
  first = seq_len(nchar(gsms)),
  last = seq_len(nchar(gsms))
)

ex <- Biobase::exprs(gset)

qx <- as.numeric(
  quantile(
    ex,
    c(0, 0.25, 0.5, 0.75, 0.99, 1),
    na.rm = TRUE
  )
)

log_required <-
  (qx[5] > 100) ||
  (
    qx[6] - qx[1] > 50 &&
    qx[2] > 0
  ) ||
  (
    qx[2] > 0 &&
    qx[2] < 1 &&
    qx[4] > 1 &&
    qx[4] < 2
  )

if (log_required) {
  ex[ex <= 0] <- NA
  Biobase::exprs(gset) <- log2(ex)
}

sml <- paste0("G", sml)
group_factor <- factor(sml)

gset$description <- group_factor

design <- model.matrix(
  ~ description + 0,
  data = Biobase::pData(gset)
)

colnames(design) <- levels(group_factor)

fit <- limma::lmFit(gset, design)

contrast_matrix <- limma::makeContrasts(
  G1 - G0,
  levels = design
)

fit2 <- limma::contrasts.fit(
  fit,
  contrast_matrix
)

fit2 <- limma::eBayes(
  fit2,
  proportion = 0.01
)

tt <- limma::topTable(
  fit2,
  adjust.method = "fdr",
  sort.by = "B",
  number = 30000
)

# 找出 GEO platform annotation 中的 gene symbol 列
symbol_candidates <- c(
  "Gene.symbol",
  "Gene.Symbol",
  "GENE_SYMBOL",
  "Symbol",
  "gene_assignment"
)

symbol_column <- symbol_candidates[
  symbol_candidates %in% colnames(tt)
][1]

if (is.na(symbol_column)) {
  stop(
    "Could not find gene-symbol column. Available columns: ",
    paste(colnames(tt), collapse = ", ")
  )
}

tt_filtered <- tt[
  !is.na(tt$logFC) &
  !is.na(tt$adj.P.Val) &
  tt$logFC < -1 &
  tt$adj.P.Val < 0.05,
  ,
  drop = FALSE
]

fezf2_symbols <- as.character(
  tt_filtered[[symbol_column]]
)

fezf2_symbols <- unlist(
  strsplit(
    fezf2_symbols,
    "///",
    fixed = TRUE
  )
)

fezf2_symbols <- trimws(fezf2_symbols)
fezf2_symbols <- fezf2_symbols[
  !is.na(fezf2_symbols) &
  nzchar(fezf2_symbols) &
  fezf2_symbols != "---"
]

fezf2_genes <- unique(fezf2_symbols)
fezf2_genes <- intersect(
  fezf2_genes,
  protein_coding
)

fezf2_genes <- intersect(
  fezf2_genes,
  seurat_genes
)

cat("fezf2_genes:", length(fezf2_genes), "\n")

write.csv(
  tt,
  file.path(outdir, "GSE69105_all_limma_results.csv"),
  row.names = TRUE
)

write.csv(
  tt_filtered,
  file.path(outdir, "GSE69105_Fezf2_filtered.csv"),
  row.names = TRUE
)

# ============================================================
# 6. 保存 Figure 2E 所需 gene_lists
# ============================================================

gene_lists <- list(
  protein_coding = sort(unique(protein_coding)),
  tra_fantom = sort(unique(tra_fantom)),
  aire_genes = sort(unique(aire_genes)),
  fezf2_genes = sort(unique(fezf2_genes))
)

output_file <- file.path(
  outdir,
  "gene_lists.rda"
)

save(
  gene_lists,
  file = output_file
)

size_table <- data.frame(
  gene_set = names(gene_lists),
  n_genes = vapply(
    gene_lists,
    length,
    integer(1)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  size_table,
  file.path(outdir, "gene_list_sizes.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("\n==============================\n")
cat("GENE LISTS COMPLETED\n")
print(size_table)
cat("Saved:", output_file, "\n")
cat("==============================\n")

if (any(size_table$n_genes == 0)) {
  stop("At least one required gene list is empty.")
}
