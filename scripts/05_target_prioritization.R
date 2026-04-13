#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

infile  <- "results/tables/deseq2_sig.tsv"
out_top <- "results/tables/target_shortlist_top50.tsv"
out_flag <- "results/tables/targets_flagged_actionable.tsv"

dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

cat("Reading:", infile, "\n")
sig <- fread(infile)
sig <- sig[!is.na(padj)]

# Basic hygiene
# Remove rows missing key stats
sig <- sig[!is.na(log2FoldChange) & !is.na(baseMean) & !is.na(pvalue)]

# Optional: avoid prioritizing genes with extremely low expression
MIN_BASEMEAN <- 10
sig <- sig[baseMean >= MIN_BASEMEAN]

# Scoring components
sig[, log2FC_abs := abs(log2FoldChange)]
sig[, padj_score := -log10(padj + 1e-300)]
sig[, expr_score := log10(baseMean + 1)]

# TargetScore (tune weights later)
# - effect size matters
# - confidence matters
# - expression matters a bit
sig[, TargetScore := (1.0 * log2FC_abs) + (0.5 * padj_score) + (0.2 * expr_score)]

# Direction tags (useful for interpretation)
sig[, direction := fifelse(log2FoldChange >= 0, "Up_in_Tumor", "Down_in_Tumor")]

# "Actionable-ish" genes starter list (expand any time)
actionable <- c(
  "EGFR","KRAS","ERBB2","MET","ALK","RET","ROS1","BRAF","NTRK1","NTRK2","NTRK3",
  "PIK3CA","PTEN","MAP2K1","MAPK1","RAF1",
  "PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","TIGIT",
  "VEGFA","KDR","FLT1",
  "FGFR1","FGFR2","FGFR3",
  "KIT","JAK1","JAK2","STAT3",
  "MTOR"
)

sig[, is_actionable_known := gene_id %chin% actionable]

# Prioritize upregulated targets first (often more “druggable” as inhibition hypotheses)
# Keep downregulated too, but they’ll rank separately (replacement/activation hypotheses)
sig_up   <- sig[direction == "Up_in_Tumor"][order(-TargetScore)]
sig_down <- sig[direction == "Down_in_Tumor"][order(-TargetScore)]

top50 <- rbind(
  sig_up[1:min(35, .N)],
  sig_down[1:min(15, .N)],
  fill = TRUE
)

# Write outputs
fwrite(top50, out_top, sep = "\t")
fwrite(sig[is_actionable_known == TRUE][order(-TargetScore)], out_flag, sep = "\t")

cat("Wrote:", out_top, "\n")
cat("Wrote:", out_flag, "\n")
cat("Top50 breakdown:",
    "Up=", nrow(top50[direction=="Up_in_Tumor"]),
    "Down=", nrow(top50[direction=="Down_in_Tumor"]), "\n")
cat("Flagged known actionable:", sig[is_actionable_known==TRUE, .N], "\n")
cat("Done.\n")
