#!/usr/bin/env bash
set -euo pipefail

# ================================
# GSE62944 data acquisition helper
# ================================

mkdir -p data/raw metadata

echo "======================================"
echo "TCGA LUAD RNA-seq project"
echo "Data source: GEO GSE62944"
echo "======================================"
echo
echo "This project uses processed TCGA counts"
echo "from GEO accession GSE62944."
echo
echo "Please manually download the following"
echo "supplementary files from the GEO page:"
echo
echo "  1) Gene-level raw count matrix"
echo "     → place into: data/raw/"
echo
echo "  2) Sample annotation / mapping file"
echo "     (sample ID → cancer type + tumor/normal)"
echo "     → place into: metadata/"
echo
echo "Expected directories:"
echo "  data/raw/"
echo "  metadata/"
echo
echo "IMPORTANT:"
echo "  - Do NOT commit data/raw/ to GitHub"
echo "  - Large GEO files should remain local"
echo
echo "Once files are in place, proceed with:"
echo "  Rscript scripts/01_build_coldata.R"
echo
echo "======================================"

