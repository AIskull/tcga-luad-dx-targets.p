#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# ---------------------------
# Config
# ---------------------------
SCRIPT_ARGS <- commandArgs(trailingOnly = FALSE)
SCRIPT_FILE <- sub("^--file=", "", SCRIPT_ARGS[grep("^--file=", SCRIPT_ARGS, value = TRUE)][1])

SCRIPT_DIR <- if (!is.na(SCRIPT_FILE) && nzchar(SCRIPT_FILE)) {
  dirname(normalizePath(SCRIPT_FILE))
} else {
  normalizePath(getwd())
}

PROJECT_ROOT <- Sys.getenv("TCGA_LUAD_ROOT")
if (!nzchar(PROJECT_ROOT)) {
  PROJECT_ROOT <- dirname(SCRIPT_DIR)
}
PROJECT_ROOT <- normalizePath(PROJECT_ROOT, mustWork = FALSE)

cat(sprintf("Using project root: %s\n", PROJECT_ROOT))

COL_PATH <- file.path(PROJECT_ROOT, "metadata/luad_coldata.tsv")

TUMOR_GZ  <- file.path(PROJECT_ROOT, "data/raw/GSE62944_RAW/GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FeatureCounts.txt.gz")
NORMAL_GZ <- file.path(PROJECT_ROOT, "data/raw/GSE62944_RAW/GSM1697009_06_01_15_TCGA_24.normal_Rsubread_FeatureCounts.txt.gz")

OUT_PATH <- file.path(PROJECT_ROOT, "data/processed/luad_counts.tsv")

# ---------------------------
# Helpers
# ---------------------------
assert_file <- function(path) {
  if (!file.exists(path)) stop(sprintf("Missing required file: %s", path))
}

assert_cmd <- function(cmd) {
  ok <- nzchar(Sys.which(cmd))
  if (!ok) stop(sprintf("Required system command not found in PATH: %s", cmd))
}

read_counts_gz <- function(path_gz) {
  # Read via zcat to avoid unzipping huge files to disk
  dt <- fread(cmd = paste("zcat", shQuote(path_gz)))

  if (ncol(dt) < 2) {
    stop(sprintf("Counts file looks wrong (too few columns): %s", path_gz))
  }

  # First column is gene id/symbol
  gene_col <- names(dt)[1]
  setnames(dt, gene_col, "gene_id")

  # Guard against duplicated gene_id early (merge would explode silently)
  if (anyDuplicated(dt$gene_id)) {
    dup_n <- sum(duplicated(dt$gene_id))
    stop(sprintf("Counts file has duplicated gene_id (%d duplicates): %s", dup_n, path_gz))
  }

  dt
}

# ---------------------------
# 1) Input checks
# ---------------------------
assert_file(COL_PATH)
assert_file(TUMOR_GZ)
assert_file(NORMAL_GZ)
assert_cmd("zcat")

# ---------------------------
# 2) Load coldata
# ---------------------------
coldata <- fread(COL_PATH)
stopifnot(all(c("sample_id", "condition") %in% names(coldata)))

# Defensive: ensure unique sample IDs in coldata
if (anyDuplicated(coldata$sample_id)) {
  dup_n <- sum(duplicated(coldata$sample_id))
  stop(sprintf("coldata has duplicated sample_id (%d duplicates): %s", dup_n, COL_PATH))
}

# Force consistent ordering: Normal then Tumor (matches DESeq2 reference)
coldata[, condition := factor(condition, levels = c("Normal", "Tumor"))]
setorder(coldata, condition, sample_id)

tumor_ids  <- coldata[condition == "Tumor",  sample_id]
normal_ids <- coldata[condition == "Normal", sample_id]

cat("Tumor samples requested:", length(tumor_ids), "\n")
cat("Normal samples requested:", length(normal_ids), "\n")

# ---------------------------
# 3) Read counts
# ---------------------------
cat("Reading tumor counts...\n")
tumor_dt <- read_counts_gz(TUMOR_GZ)

cat("Reading normal counts...\n")
normal_dt <- read_counts_gz(NORMAL_GZ)

# ---------------------------
# 4) Subset to requested samples that exist
# ---------------------------
tumor_keep  <- intersect(tumor_ids,  names(tumor_dt))
normal_keep <- intersect(normal_ids, names(normal_dt))

missing_tumor  <- setdiff(tumor_ids,  tumor_keep)
missing_normal <- setdiff(normal_ids, normal_keep)

cat("Tumor samples found in matrix:", length(tumor_keep), "\n")
cat("Normal samples found in matrix:", length(normal_keep), "\n")

if (length(missing_tumor) > 0) {
  cat("Warning: Missing tumor samples (showing up to 5):\n")
  print(head(missing_tumor, 5))
}
if (length(missing_normal) > 0) {
  cat("Warning: Missing normal samples (showing up to 5):\n")
  print(head(missing_normal, 5))
}

# Keep only gene_id + kept sample columns
tumor_sub  <- tumor_dt[,  c("gene_id", tumor_keep),  with = FALSE]
normal_sub <- normal_dt[, c("gene_id", normal_keep), with = FALSE]

# ---------------------------
# 5) Merge tumor + normal by gene_id
# ---------------------------
cat("Merging tumor + normal by gene_id...\n")
merged <- merge(tumor_sub, normal_sub, by = "gene_id", all = FALSE)

# Sanity
stopifnot(!anyDuplicated(merged$gene_id))

# ---------------------------
# 6) Re-order sample columns to match coldata order exactly
# ---------------------------
final_sample_order <- coldata$sample_id
final_sample_order <- final_sample_order[final_sample_order %in% names(merged)]

merged <- merged[, c("gene_id", final_sample_order), with = FALSE]

cat("Final matrix genes:", nrow(merged), "\n")
cat("Final matrix samples:", ncol(merged) - 1, "\n")

# ---------------------------
# 7) Write output
# ---------------------------
dir.create(dirname(OUT_PATH), showWarnings = FALSE, recursive = TRUE)
fwrite(merged, OUT_PATH, sep = "\t")
cat("Wrote:", OUT_PATH, "\n")
