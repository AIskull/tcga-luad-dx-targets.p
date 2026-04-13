#!/usr/bin/env bash
set -euo pipefail

# Resolve project root from the script location unless overridden.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${TCGA_LUAD_ROOT:-$SCRIPT_DIR/..}"
cd "$PROJECT_ROOT"

echo "== Running TCGA LUAD pipeline =="
echo "Project root: $PROJECT_ROOT"
mkdir -p results/tables results/figures

# ---- required local files (not tracked in git) ----
RAW_FAMILY_TAR="$PROJECT_ROOT/data/raw/GSE62944_family.tar"
SAMPLE_MAP="$PROJECT_ROOT/metadata/sample_map.tsv"

COUNTS="$PROJECT_ROOT/data/processed/luad_counts.tsv"
COLDATA="$PROJECT_ROOT/metadata/luad_coldata.tsv"

# If processed inputs exist, skip the manual-download message.
if [[ -f "$COUNTS" && -f "$COLDATA" ]]; then
  echo ""
  echo "Found existing local inputs:"
  echo "  - $COUNTS"
  echo "  - $COLDATA"
  echo "Skipping GEO download instructions."
  echo ""
else
  echo ""
  echo "Processed inputs not found."
  echo "Will attempt to download/build required files."
  echo ""
fi

# Step 0: Download raw GEO files only if missing
if [[ ! -f "$RAW_FAMILY_TAR" || ! -f "$SAMPLE_MAP" ]]; then
  echo "Raw GEO files missing; running downloader..."
  bash scripts/00_download_gse62944.sh
else
  echo "Raw GEO files already present; skipping download."
fi

# Build processed files if missing
if [[ ! -f "$COLDATA" ]]; then
  Rscript scripts/01_build_coldata.R
fi

if [[ ! -f "$COUNTS" ]]; then
  Rscript scripts/02_subset_counts.R
fi

# Safety: if still missing, stop early with a clear error
if [[ ! -f "$COUNTS" || ! -f "$COLDATA" ]]; then
  echo ""
  echo "ERROR: Required files still missing after build steps:"
  [[ ! -f "$COUNTS" ]] && echo "  - $COUNTS"
  [[ ! -f "$COLDATA" ]] && echo "  - $COLDATA"
  echo "Check your GEO downloads and paths."
  exit 1
fi

# Main analysis
Rscript scripts/03_deseq2_LUAD.R
Rscript scripts/05_target_prioritization.R
python scripts/06_panel_rank_stability.py
Rscript scripts/07_make_summary_tables.R
Rscript scripts/08_boxplots_key_genes.R
python scripts/09_external_validation_cptac_luad.py

echo "== Done. Outputs in results/ =="
