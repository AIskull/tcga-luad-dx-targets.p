#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(ggplot2)
})

FIG_DIR <- "results/figures"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

COUNTS_PATH  <- "data/processed/luad_counts.tsv"
COLDATA_PATH <- "metadata/luad_coldata.tsv"

PANEL_RANKSTABLE <- "results/tables/diagnostic_panel_rankstable_top20.tsv"
HIGHLIGHT_PATH   <- "results/tables/highlight_genes.tsv"
DESEQ_SIG_PATH   <- "results/tables/deseq2_sig.tsv"

message("Reading counts: ", COUNTS_PATH)
counts_dt <- fread(COUNTS_PATH)
gene_col <- names(counts_dt)[1]
gene_ids <- counts_dt[[1]]

counts_mat <- as.matrix(counts_dt[, -1, with = FALSE])
rownames(counts_mat) <- gene_ids

message("Reading coldata: ", COLDATA_PATH)
coldata <- fread(COLDATA_PATH)
stopifnot(all(c("sample_id", "condition") %in% names(coldata)))
coldata <- as.data.frame(coldata)
rownames(coldata) <- coldata$sample_id

# Align sample order
common <- intersect(colnames(counts_mat), rownames(coldata))
message("Common samples: ", length(common))
stopifnot(length(common) > 0)

counts_mat <- counts_mat[, common, drop = FALSE]
coldata <- coldata[common, , drop = FALSE]

# Build DESeq object only for normalization (fast)
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = coldata,
  design    = ~ condition
)
dds <- estimateSizeFactors(dds)
norm <- counts(dds, normalized = TRUE)
log_norm <- log2(norm + 1)

# Choose genes to plot (priority: rank-stable panel -> highlight -> top sig)
genes_to_plot <- character(0)

if (file.exists(PANEL_RANKSTABLE)) {
  message("Found rank-stable panel: ", PANEL_RANKSTABLE)
  panel <- fread(PANEL_RANKSTABLE)
  if ("gene" %in% names(panel)) {
    genes_to_plot <- unique(panel$gene)
  } else if ("gene_id" %in% names(panel)) {
    genes_to_plot <- unique(panel$gene_id)
  }
}

if (length(genes_to_plot) == 0 && file.exists(HIGHLIGHT_PATH)) {
  message("No panel genes found; using highlight genes: ", HIGHLIGHT_PATH)
  hg <- fread(HIGHLIGHT_PATH)
  if ("gene_id" %in% names(hg)) genes_to_plot <- unique(hg$gene_id)
}

if (length(genes_to_plot) == 0 && file.exists(DESEQ_SIG_PATH)) {
  message("No panel/highlight found; using top DESeq2 sig genes: ", DESEQ_SIG_PATH)
  sig <- fread(DESEQ_SIG_PATH)
  sig <- sig[!is.na(padj)]
  sig <- sig[order(padj)]
  genes_to_plot <- unique(sig$gene_id[1:min(10, nrow(sig))])
}

# Keep only genes that exist in matrix
genes_to_plot <- genes_to_plot[genes_to_plot %in% rownames(log_norm)]
if (length(genes_to_plot) == 0) {
  stop("No requested genes exist in the counts matrix. Check gene IDs in your tables vs counts rownames.")
}

message("Genes to plot (n=", length(genes_to_plot), "): ", paste(genes_to_plot, collapse = ", "))

plot_one <- function(gene) {
  df <- data.frame(
    sample    = colnames(log_norm),
    expr      = as.numeric(log_norm[gene, ]),
    condition = coldata$condition,
    stringsAsFactors = FALSE
  )

  p <- ggplot(df, aes(x = condition, y = expr)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.35, size = 1.2) +
    labs(
      title = paste0("LUAD Tumor vs Normal: ", gene),
      x = NULL,
      y = "log2(normalized count + 1)"
    ) +
    theme_bw(base_size = 12)

  out <- file.path(FIG_DIR, paste0("boxplot_", gene, ".png"))
  ggsave(out, p, width = 6.5, height = 4.5, dpi = 150)
  message("Saved: ", out)
}

for (g in genes_to_plot) {
  tryCatch(plot_one(g), error = function(e) message("Skipping ", g, ": ", e$message))
}

message("Done.")
