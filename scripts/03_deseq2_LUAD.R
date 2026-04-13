#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(apeglm)
})

# -----------------------------
# Config
# -----------------------------
COUNTS_PATH  <- "data/processed/luad_counts.tsv"
COLDATA_PATH <- "metadata/luad_coldata.tsv"

OUT_ALL_SHR  <- "results/tables/deseq2_all.tsv"
OUT_ALL_RAW  <- "results/tables/deseq2_all_unshrunk.tsv"
OUT_SIG      <- "results/tables/deseq2_sig.tsv"

MA_PNG    <- "results/figures/ma_plot.png"
PCA_PNG   <- "results/figures/pca.png"          # match README checklist
VOL_PNG   <- "results/figures/volcano.png"
HEAT_PNG  <- "results/figures/heatmap_top50.png"

LOW_COUNT_MIN_TOTAL <- 10      # keep genes with >= this total count across all samples
TOP_HEATMAP_N <- 50

# -----------------------------
# Helpers
# -----------------------------
assert_file <- function(path) {
  if (!file.exists(path)) stop(sprintf("Missing required file: %s", path))
}

dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1) Input checks
# -----------------------------
assert_file(COUNTS_PATH)
assert_file(COLDATA_PATH)

message("Reading counts: ", COUNTS_PATH)
counts_dt <- fread(COUNTS_PATH)

if (ncol(counts_dt) < 3) stop("Counts matrix has too few columns (need gene_id + >=2 samples).")

gene_col <- names(counts_dt)[1]
setnames(counts_dt, gene_col, "gene_id")

if (anyDuplicated(counts_dt$gene_id)) {
  dup_n <- sum(duplicated(counts_dt$gene_id))
  stop(sprintf("Counts matrix has duplicated gene_id (%d duplicates). Fix before DESeq2.", dup_n))
}

# counts matrix: genes x samples
gene_ids <- counts_dt$gene_id
counts_mat <- as.matrix(counts_dt[, -1, with = FALSE])
rownames(counts_mat) <- gene_ids

message("Reading coldata: ", COLDATA_PATH)
coldata <- fread(COLDATA_PATH)
stopifnot(all(c("sample_id", "condition") %in% names(coldata)))

# enforce reference level
coldata[, condition := factor(condition, levels = c("Normal", "Tumor"))]

# -----------------------------
# 2) Align samples (strict)
# -----------------------------
samples_in_counts <- colnames(counts_mat)
samples_in_coldata <- coldata$sample_id

common <- intersect(samples_in_coldata, samples_in_counts)

message("Samples in coldata: ", nrow(coldata))
message("Samples in counts: ", length(samples_in_counts))
message("Common samples: ", length(common))

if (length(common) < 2) stop("Too few overlapping samples between coldata and counts.")

missing_in_counts <- setdiff(samples_in_coldata, samples_in_counts)
missing_in_coldata <- setdiff(samples_in_counts, samples_in_coldata)

if (length(missing_in_counts) > 0) {
  message("Warning: samples present in coldata but missing in counts (showing up to 5):")
  print(head(missing_in_counts, 5))
}
if (length(missing_in_coldata) > 0) {
  message("Warning: samples present in counts but missing in coldata (showing up to 5):")
  print(head(missing_in_coldata, 5))
}

# keep coldata order, subset to common
coldata2 <- as.data.frame(coldata[sample_id %in% common])
rownames(coldata2) <- coldata2$sample_id
counts_mat2 <- counts_mat[, rownames(coldata2), drop = FALSE]

stopifnot(identical(colnames(counts_mat2), rownames(coldata2)))

# -----------------------------
# 3) Sanity + low-count filter
# -----------------------------
# force numeric -> integer-ish
storage.mode(counts_mat2) <- "numeric"

if (any(is.na(counts_mat2))) stop("Counts matrix contains NA values.")
if (any(counts_mat2 < 0)) stop("Counts matrix contains negative values (should not happen).")

# low-count filter
keep <- rowSums(counts_mat2) >= LOW_COUNT_MIN_TOTAL
counts_mat2 <- counts_mat2[keep, , drop = FALSE]
message("Genes kept after low-count filter: ", nrow(counts_mat2))

# ensure integer mode for DESeq2
counts_mat2 <- round(counts_mat2)
storage.mode(counts_mat2) <- "integer"

# -----------------------------
# 4) Run DESeq2
# -----------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat2,
  colData   = coldata2,
  design    = ~ condition
)

dds <- DESeq(dds)

# raw results
res_raw <- results(dds, contrast = c("condition", "Tumor", "Normal"))

# shrink LFC (better for ranking/volcano readability)
res_shr <- lfcShrink(dds, coef = "condition_Tumor_vs_Normal", type = "apeglm")

# -----------------------------
# 5) Write tables
# -----------------------------
res_raw_df <- as.data.frame(res_raw)
res_raw_df$gene_id <- rownames(res_raw_df)
res_raw_df <- res_raw_df[, c("gene_id", setdiff(colnames(res_raw_df), "gene_id"))]
fwrite(as.data.table(res_raw_df), OUT_ALL_RAW, sep = "\t")
message("Wrote: ", OUT_ALL_RAW)

res_df <- as.data.frame(res_shr)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[, c("gene_id", setdiff(colnames(res_df), "gene_id"))]
fwrite(as.data.table(res_df), OUT_ALL_SHR, sep = "\t")
message("Wrote: ", OUT_ALL_SHR)

sig_df <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 1, ]
sig_df <- sig_df[order(sig_df$padj), ]
fwrite(as.data.table(sig_df), OUT_SIG, sep = "\t")
message("Wrote: ", OUT_SIG)
message("Significant genes (padj<0.05 & |log2FC|>=1): ", nrow(sig_df))

# -----------------------------
# 6) Plots
# -----------------------------
# MA
png(MA_PNG, width = 1200, height = 900, res = 150)
plotMA(res_shr, ylim = c(-6, 6), main = "DESeq2 MA Plot (LUAD Tumor vs Normal)")
dev.off()
message("Saved: ", MA_PNG)

# PCA (VST)
vsd <- vst(dds, blind = TRUE)
pca_df <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_df, "percentVar"))

p <- ggplot(pca_df, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3, alpha = 0.85) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA (VST) - LUAD Tumor vs Normal") +
  theme_minimal()

ggsave(PCA_PNG, plot = p, width = 8, height = 6, dpi = 150)
message("Saved: ", PCA_PNG)

# Volcano (uses adjusted p-values; avoid infinite)
vol <- res_df
vol$neglog10padj <- -log10(vol$padj)
vol$neglog10padj[is.infinite(vol$neglog10padj)] <- NA

vol$class <- "Not Sig"
vol$class[!is.na(vol$padj) & vol$padj < 0.05 & vol$log2FoldChange >= 1]  <- "Up"
vol$class[!is.na(vol$padj) & vol$padj < 0.05 & vol$log2FoldChange <= -1] <- "Down"

vp <- ggplot(vol, aes(log2FoldChange, neglog10padj, color = class)) +
  geom_point(alpha = 0.6, size = 1) +
  ggtitle("Volcano Plot (LUAD Tumor vs Normal)") +
  xlab("log2 Fold Change (Tumor/Normal)") +
  ylab("-log10 adjusted p-value") +
  theme_minimal()

ggsave(VOL_PNG, plot = vp, width = 8, height = 6, dpi = 150)
message("Saved: ", VOL_PNG)

# Heatmap of top N by padj (use RAW for padj ordering, plot VST expression)
res_raw_ordered <- res_raw[order(res_raw$padj), ]
top_genes <- rownames(res_raw_ordered)[!is.na(res_raw_ordered$padj)]
top_genes <- head(top_genes, TOP_HEATMAP_N)

mat <- assay(vsd)[top_genes, , drop = FALSE]
mat_z <- t(scale(t(mat)))  # row z-score

ann <- data.frame(condition = coldata2$condition)
rownames(ann) <- rownames(coldata2)

png(HEAT_PNG, width = 1400, height = 1000, res = 150)
pheatmap(
  mat_z,
  annotation_col = ann,
  show_colnames = FALSE,
  fontsize_row = 7,
  main = paste0("Top ", length(top_genes), " DE Genes (VST Z-score)")
)
dev.off()
message("Saved: ", HEAT_PNG)

message("Done.")
