#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

INPUT_MAP  <- "metadata/sample_map_all.tsv"
OUTPUT_COL <- "metadata/luad_coldata.tsv"

# ---------------------------
# 1) Basic file checks
# ---------------------------
if (!file.exists(INPUT_MAP)) {
  stop(sprintf(
    "Missing input file: %s\nPut the GEO sample mapping file there and re-run.",
    INPUT_MAP
  ))
}

# ---------------------------
# 2) Read mapping file
#    (GEO file often has NO header)
# ---------------------------
map <- fread(INPUT_MAP, header = FALSE)

if (ncol(map) < 2) {
  stop(sprintf(
    "Mapping file %s has %d columns; expected at least 2 (sample_id, cancer_type).",
    INPUT_MAP, ncol(map)
  ))
}

# Keep only first two columns in case extra columns exist
map <- map[, .(V1, V2)]
setnames(map, c("sample_id", "cancer_type"))

# ---------------------------
# 3) Derive Tumor/Normal from TCGA barcode
#    -01A- = Tumor, -11A- = Normal
# ---------------------------
map[, sample_type := fifelse(
  grepl("-11A-", sample_id), "Normal",
  fifelse(grepl("-01A-", sample_id), "Tumor", NA_character_)
)]

# ---------------------------
# 4) Filter to LUAD tumor + normal
# ---------------------------
luad <- map[cancer_type == "LUAD"]

cat("Total rows in map:", nrow(map), "\n")
cat("Rows labeled LUAD:", nrow(luad), "\n")

# Drop rows that aren't -01A- or -11A-
na_before <- nrow(luad)
luad <- luad[!is.na(sample_type)]
cat("Dropped LUAD rows with unknown sample_type:", na_before - nrow(luad), "\n")

# Keep only Tumor/Normal (defensive)
luad <- luad[sample_type %in% c("Tumor", "Normal")]

# Deduplicate sample IDs (defensive)
dup_n <- sum(duplicated(luad$sample_id))
if (dup_n > 0) {
  cat("Warning: duplicated sample_id rows found:", dup_n, "→ keeping first occurrence.\n")
  setorder(luad, sample_id)
  luad <- luad[!duplicated(sample_id)]
}

# ---------------------------
# 5) Build DESeq2 colData
# ---------------------------
coldata <- data.table(
  sample_id = luad$sample_id,
  condition = factor(luad$sample_type, levels = c("Normal", "Tumor"))
)

# Deterministic order
setorder(coldata, condition, sample_id)

cat("LUAD samples kept:", nrow(coldata), "\n")
print(table(coldata$condition))

# ---------------------------
# 6) Write output
# ---------------------------
dir.create(dirname(OUTPUT_COL), showWarnings = FALSE, recursive = TRUE)
fwrite(coldata, OUTPUT_COL, sep = "\t")
cat("Wrote:", OUTPUT_COL, "\n")

