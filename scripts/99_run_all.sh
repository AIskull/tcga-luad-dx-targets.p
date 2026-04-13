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
RAW_TUMOR="$PROJECT_ROOT/data/raw/GSE62944_RAW/GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FeatureCounts.txt.gz"
RAW_NORMAL="$PROJECT_ROOT/data/raw/GSE62944_RAW/GSM1697009_06_01_15_TCGA_24.normal_Rsubread_FeatureCounts.txt.gz"
SAMPLE_MAP="$PROJECT_ROOT/metadata/sample_map_all.tsv"
DATA_SOURCES_FILE="$PROJECT_ROOT/data_sources.tsv"

COUNTS="$PROJECT_ROOT/data/processed/luad_counts.tsv"
COLDATA="$PROJECT_ROOT/metadata/luad_coldata.tsv"

print_status() {
  local label="$1"
  local path="$2"
  if [[ -f "$path" ]]; then
    echo "  OK $label: $path"
    return 0
  fi
  echo "  MISSING $label: $path"
  return 1
}

verify_optional_checksums() {
  local manifest="$1"
  if [[ ! -f "$manifest" ]]; then
    echo "No checksum manifest found: $manifest"
    echo "Running with checksum validation disabled."
    return 0
  fi

  echo "Verifying optional checksums from: $manifest"
  local failed=0
  while IFS=$'\t' read -r dataset accession filename local_path download_url sha256; do
    [[ -z "$dataset" || "${dataset:0:1}" == "#" ]] && continue
    [[ "$dataset" == "dataset" ]] && continue
    [[ -z "$local_path" ]] && continue
    if [[ -z "$sha256" ]]; then
      echo "WARNING: missing checksum for: $local_path"
      continue
    fi
    local full_path="$PROJECT_ROOT/$local_path"
    if [[ ! -f "$full_path" ]]; then
      echo "WARNING: checksum entry without local file: $full_path"
      continue
    fi
    local actual
    actual="$(sha256sum "$full_path" | awk '{print $1}')"
    if [[ "$actual" != "$sha256" ]]; then
      echo "ERROR: checksum mismatch for $full_path"
      echo "  expected: $sha256"
      echo "  actual:   $actual"
      failed=1
    else
      echo "  checksum OK: $local_path"
    fi
  done < "$manifest"

  if [[ "$failed" -ne 0 ]]; then
    echo "Checksum validation failed."
    echo "Please re-download affected files and retry."
    exit 1
  fi
  echo "Checksum verification complete."
}

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

# Step 0: Validate required input files
echo "Checking required inputs:"
if ! print_status "raw tumor matrix" "$RAW_TUMOR"; then
  echo "Run scripts/00_download_gse62944.sh to fetch missing files."
  bash scripts/00_download_gse62944.sh
  if ! print_status "raw tumor matrix" "$RAW_TUMOR"; then
    exit 1
  fi
fi

if ! print_status "raw normal matrix" "$RAW_NORMAL"; then
  echo "Run scripts/00_download_gse62944.sh to fetch missing files."
  bash scripts/00_download_gse62944.sh
  if ! print_status "raw normal matrix" "$RAW_NORMAL"; then
    exit 1
  fi
fi

if ! print_status "sample map" "$SAMPLE_MAP"; then
  echo "ERROR: required mapping file is missing: $SAMPLE_MAP"
  echo "This repository ships with metadata/sample_map_all.tsv; verify it was not removed."
  exit 1
fi

if [[ -f "$DATA_SOURCES_FILE" ]]; then
  verify_optional_checksums "$DATA_SOURCES_FILE"
else
  echo "No data_sources.tsv manifest found; skipping optional checksum checks."
fi

# Build processed files if missing
if [[ ! -f "$COLDATA" ]]; then
  echo "Building metadata with: scripts/01_build_coldata.R"
  Rscript scripts/01_build_coldata.R
fi

if [[ ! -f "$COUNTS" ]]; then
  echo "Building processed counts with: scripts/02_subset_counts.R"
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
