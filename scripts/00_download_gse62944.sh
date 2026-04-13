#!/usr/bin/env bash
set -euo pipefail

# ================================
# GSE62944 data acquisition helper
# ================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${TCGA_LUAD_ROOT:-$SCRIPT_DIR/..}"
DATA_SOURCES_FILE="${PROJECT_ROOT}/data_sources.tsv"
RAW_DIR="$PROJECT_ROOT/data/raw/GSE62944_RAW"

mkdir -p "$RAW_DIR"

print_missing() {
  echo "  Missing: $1"
  echo "  └─ expected at: $2"
  echo
}

echo "======================================"
echo "TCGA LUAD RNA-seq project"
echo "Data source: GEO accession GSE62944"
echo "======================================"
echo
echo "What this repo expects (tracked results are safe to publish):"
echo
echo "  - Processed repository outputs are tracked."
echo "  - Raw GEO files stay local and are intentionally not tracked."
echo "  - CPTAC validation is downloaded through the script dependency at runtime."
echo

if [[ -f "$DATA_SOURCES_FILE" ]]; then
  echo "Public source list: $DATA_SOURCES_FILE"
  echo
  echo "Required GEO files:"
  while IFS=$'\t' read -r dataset accession filename local_path download_url sha256; do
    [[ -z "$dataset" || "${dataset:0:1}" == "#" ]] && continue
    [[ "$dataset" == "dataset" ]] && continue
    [[ -z "$local_path" ]] && continue

    echo "  - ${filename} ($dataset)"
    if [[ -f "${PROJECT_ROOT}/${local_path}" ]]; then
      echo "    status: present"
    else
      echo "    status: missing"
    fi
    echo "    source: ${download_url:-not listed}"
    if [[ -n "${sha256:-}" ]]; then
      echo "    checksum: sha256 provided in manifest"
    else
      echo "    checksum: not provided (manifest line has no hash)"
    fi
    echo
  done < "$DATA_SOURCES_FILE"
else
  echo "No data source manifest found:"
  echo "  $DATA_SOURCES_FILE"
  echo
  echo "Expected files:"
  print_missing "TCGA LUAD tumor count matrix" "$RAW_DIR/GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FeatureCounts.txt.gz"
  print_missing "TCGA LUAD normal count matrix" "$RAW_DIR/GSM1697009_06_01_15_TCGA_24.normal_Rsubread_FeatureCounts.txt.gz"
fi

echo "If you have the files available, place them under:"
echo "  $RAW_DIR"
echo
echo "If you would like an optional checksum check, add a sha256 value"
echo "in the manifest's sha256 column for any row."
echo

echo "To continue the pipeline once downloads are in place:"
echo "  Rscript scripts/01_build_coldata.R"
echo

if [[ "${1:-}" == "--download" ]]; then
  echo "Attempting to download missing files..."
  if ! command -v curl >/dev/null 2>&1; then
    echo "curl is required for --download. Install curl and rerun."
    echo "Manual alternative:"
    echo "  1) open https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62944"
    echo "  2) download the two FeatureCounts supplementary files"
    exit 1
  fi

  if [[ ! -f "$DATA_SOURCES_FILE" ]]; then
    echo "No manifest found; cannot auto-download."
    echo "Create $DATA_SOURCES_FILE with columns: dataset accession filename local_path download_url sha256"
    exit 1
  fi

  while IFS=$'\t' read -r dataset accession filename local_path download_url sha256; do
    [[ -z "$dataset" || "${dataset:0:1}" == "#" ]] && continue
    [[ "$dataset" == "dataset" ]] && continue
    [[ -z "$local_path" || -z "${download_url:-}" ]] && continue

    OUT_PATH="${PROJECT_ROOT}/${local_path}"
    mkdir -p "$(dirname "$OUT_PATH")"
    if [[ -f "$OUT_PATH" ]]; then
      echo "Skipping (already present): $OUT_PATH"
      continue
    fi
    echo "Downloading: $download_url"
    curl -L --fail --retry 2 -o "$OUT_PATH" "$download_url"
  done < "$DATA_SOURCES_FILE"

  echo "Download pass complete."
  echo "Optional: verify checksums in scripts/99_run_all.sh"
fi

echo
echo "======================================"
