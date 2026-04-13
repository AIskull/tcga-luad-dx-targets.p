#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------
# Config
# -----------------------------
ALL_PATH <- "results/tables/deseq2_all.tsv"
SIG_PATH <- "results/tables/deseq2_sig.tsv"

OUT_DIR <- "results/tables"
FIG_DIR <- "results/figures"  # not used here, but kept for consistency
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1.0

# "Highlight" thresholds (for boxplot candidates)
HIGHLIGHT_N            <- 25
HIGHLIGHT_MIN_BASEMEAN <- 50
HIGHLIGHT_EXTREME_LFC  <- 8

# -----------------------------
# Helpers
# -----------------------------
safe_fread <- function(path) {
  if (!file.exists(path)) {
    stop("Missing file: ", path)
  }
  fread(path, sep = "\t", data.table = TRUE)
}

pick_panel_file <- function() {
  candidates <- c(
    "results/tables/diagnostic_panel_rankstable_top20.tsv",
    "results/tables/diagnostic_panel_stable_top20.tsv",
    "results/tables/diagnostic_panel_top20_strict.tsv",
    "results/tables/diagnostic_panel_top20.tsv"
  )
  for (f in candidates) {
    if (file.exists(f)) return(f)
  }
  return(NA_character_)
}

write_tsv <- function(dt, path) {
  fwrite(dt, path, sep = "\t", quote = FALSE)
  cat("Wrote:", path, "\n")
}

# -----------------------------
# Main
# -----------------------------
cat("Reading:", ALL_PATH, "\n")
all <- safe_fread(ALL_PATH)
all <- all[!is.na(padj)]

# Ensure direction column exists
all[, direction := ifelse(log2FoldChange > 0, "Up_in_Tumor", "Down_in_Tumor")]

cat("Reading:", SIG_PATH, "\n")
sig <- safe_fread(SIG_PATH)
sig <- sig[!is.na(padj)]
sig[, direction := ifelse(log2FoldChange > 0, "Up_in_Tumor", "Down_in_Tumor")]

cat("Building summary numbers...\n")

# Padj-only stats from ALL results
n_all_padj <- all[padj < PADJ_CUTOFF, .N]
n_all_up   <- all[padj < PADJ_CUTOFF & log2FoldChange > 0, .N]
n_all_down <- all[padj < PADJ_CUTOFF & log2FoldChange < 0, .N]

# Strict stats (padj + absLFC) — use SIG (usually already strict-filtered, but we compute anyway)
n_strict <- sig[padj < PADJ_CUTOFF & abs(log2FoldChange) >= LFC_CUTOFF, .N]
n_strict_up <- sig[padj < PADJ_CUTOFF & log2FoldChange >= LFC_CUTOFF, .N]
n_strict_down <- sig[padj < PADJ_CUTOFF & log2FoldChange <= -LFC_CUTOFF, .N]

summary_numbers <- rbindlist(list(
  data.table(metric = "n_all_padj<0.05", value = n_all_padj),
  data.table(metric = "n_all_up_padj<0.05", value = n_all_up),
  data.table(metric = "n_all_down_padj<0.05", value = n_all_down),

  data.table(metric = paste0("n_sig_strict_padj<", PADJ_CUTOFF, "_absLFC>=", LFC_CUTOFF), value = n_strict),
  data.table(metric = paste0("n_sig_up_strict_padj<", PADJ_CUTOFF, "_LFC>=", LFC_CUTOFF), value = n_strict_up),
  data.table(metric = paste0("n_sig_down_strict_padj<", PADJ_CUTOFF, "_LFC<=", -LFC_CUTOFF), value = n_strict_down)
))

write_tsv(summary_numbers, file.path(OUT_DIR, "summary_numbers.tsv"))

cat("Selecting top 10 up/down by padj (from SIG)...\n")
# Top hits: from SIG
top_up <- sig[log2FoldChange > 0][order(padj)][1:10]
top_down <- sig[log2FoldChange < 0][order(padj)][1:10]

write_tsv(top_up, file.path(OUT_DIR, "top10_up.tsv"))
write_tsv(top_down, file.path(OUT_DIR, "top10_down.tsv"))

cat("Selecting highlight genes (boxplot candidates)...\n")
# Highlight genes = strict + decent expression + either very small padj or extreme effect
highlight <- sig[
  baseMean >= HIGHLIGHT_MIN_BASEMEAN &
    (abs(log2FoldChange) >= HIGHLIGHT_EXTREME_LFC | padj <= 1e-20)
][order(padj)][1:HIGHLIGHT_N]

# Add reason tag (simple, readable)
highlight[, highlight_reason := fifelse(
  abs(log2FoldChange) >= HIGHLIGHT_EXTREME_LFC,
  "strict + extreme_effect",
  "strict + very_small_padj"
)]

write_tsv(highlight, file.path(OUT_DIR, "highlight_genes.tsv"))

# Panel-derived convenience tables (optional)
panel_file <- pick_panel_file()
if (!is.na(panel_file)) {
  cat("Found diagnostic panel:", panel_file, "\n")
  panel <- safe_fread(panel_file)

  # Expect columns: gene, coef (or similar)
  # Normalize column names just in case
  cn <- names(panel)
  if (!("gene" %in% cn) && ("gene_id" %in% cn)) setnames(panel, "gene_id", "gene")
  if (!("coef" %in% cn) && ("coefficient" %in% cn)) setnames(panel, "coefficient", "coef")

  if (all(c("gene", "coef") %in% names(panel))) {
    panel_up <- panel[order(-coef)][1:min(10, .N)]
    panel_down <- panel[order(coef)][1:min(10, .N)]
    write_tsv(panel_up, file.path(OUT_DIR, "top_up_from_panel.tsv"))
    write_tsv(panel_down, file.path(OUT_DIR, "top_down_from_panel.tsv"))
  } else {
    cat("Panel file exists but doesn't have expected columns {gene, coef}. Skipping panel summaries.\n")
  }
} else {
  cat("No diagnostic panel found. Skipping panel-derived tables.\n")
}

cat("Done.\n")
