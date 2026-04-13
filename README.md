# TCGA LUAD RNA-seq Analysis: Diagnostic Biomarkers & Therapeutic Target Discovery

## Overview

Lung adenocarcinoma (LUAD) is the most common subtype of lung cancer and a leading cause of cancer-related mortality worldwide. Despite advances in targeted therapies and immunotherapy, early detection and molecular stratification remain critical challenges.

RNA sequencing (RNA-seq) enables systematic characterization of transcriptional differences between tumor and normal tissue. These differences can be leveraged to:

- Identify **diagnostic biomarkers**
- Build **classification models** distinguishing tumor from normal samples
- Prioritize **therapeutic targets** based on expression shifts and biological relevance

This repository implements an **end-to-end, fully reproducible RNA-seq analysis pipeline** using TCGA LUAD data, with both **diagnostic modeling** and **drug-discovery–oriented analyses**, followed by **external validation** on an independent cohort.

---

## Data Sources

### TCGA (Training Cohort)
- **Dataset:** GEO **GSE62944**
- **Origin:** The Cancer Genome Atlas (TCGA)
- **Data type:** Gene-level RNA-seq count matrix
- **Samples:**
  - 497 LUAD tumor samples
  - 58 matched normal lung samples
- **Genes analyzed:** ~22,500 after filtering
- **Required local raw files (public):**
  - `data/raw/GSE62944_RAW/GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FeatureCounts.txt.gz`
  - `data/raw/GSE62944_RAW/GSM1697009_06_01_15_TCGA_24.normal_Rsubread_FeatureCounts.txt.gz`
- **Mapping file used for LUAD filtering:** `metadata/sample_map_all.tsv` (tracked in this repo)

### CPTAC (External Validation Cohort)
- **Dataset:** CPTAC LUAD
- **Samples:**
  - Tumor: 110
  - Normal-adjacent tissue (NAT): 211
- **Platform:** RNA-seq (WashU FPKM)

> Raw data files are not committed to GitHub. See **Reproducibility** for instructions.

---

## Pipeline Overview

The analysis follows a modular, script-based workflow:

1. **Data acquisition**
   - Download TCGA counts and sample annotations

2. **Metadata construction**
   - LUAD sample identification
   - Tumor vs Normal labeling

3. **Count matrix processing**
   - LUAD-only subsetting
   - Tumor/normal merging

4. **Differential expression analysis**
   - DESeq2 modeling
   - Log2 fold-change shrinkage (apeglm)
   - Multiple testing correction

5. **Exploratory analysis**
   - MA plot
   - PCA (variance-stabilized counts)
   - Volcano plot
   - Heatmap of top DE genes

6. **Diagnostic classifier development**
   - Penalized logistic regression
   - Rank-stability feature selection
   - ROC/AUC evaluation

7. **Therapeutic target prioritization**
   - Effect size and confidence-weighted scoring
   - Annotation of known actionable cancer targets

8. **Interpretability & reporting**
   - Summary tables
   - Expression boxplots for key biomarkers

---

## Key Results

### Differential Expression Summary

| Metric | Count |
|------|------|
| Significant DEGs (padj < 0.05) | 15,682 |
| Strict DEGs (padj < 0.05 & \|log2FC\| ≥ 1) | 6,287 |
| Upregulated in tumor | 4,217 |
| Downregulated in tumor | 2,070 |

These results demonstrate extensive transcriptional reprogramming in LUAD tumors.

---

## Exploratory Visualizations

- **MA plot**: global expression shifts  
- **PCA (VST)**: strong tumor–normal separation  
- **Volcano plot**: effect size vs statistical significance  
- **Heatmap (top 50 DE genes)**: consistent tumor-specific expression patterns  

Collectively, these plots confirm robust biological signal and cohort separability.

---

## Diagnostic Gene Panel

A **rank-stability cross-validation framework** was used to identify genes that remain consistently informative across folds.

- **Model:** Penalized logistic regression
- **Features:** Top 20 rank-stable genes
- **Internal performance (TCGA):**
  - Mean ROC AUC ≈ **0.9998**

**Representative diagnostic genes:**
- *SEMA5B*
- *FOXM1*
- *SPP1*
- *AFAP1-AS1*

Expression boxplots demonstrate clear, biologically interpretable tumor–normal shifts.

---

## External Validation: CPTAC LUAD

To assess cross-cohort generalization, the TCGA-trained diagnostic model was evaluated on an **independent CPTAC LUAD dataset**.

### Validation Strategy
- Model trained **exclusively on TCGA**
- No retraining or fine-tuning on CPTAC
- Fixed rank-stable gene panel transferred across cohorts

### External Performance
- **ROC AUC:** 0.739
- **Genes transferred:** 19 / 20

**Default threshold (0.5):**
- Sensitivity: 81.8%
- Specificity: 57.3%
- Accuracy: 65.7%
- Balanced accuracy: 69.6%

**Optimized threshold (Youden’s J):**
- Sensitivity: 100%
- Specificity: 47.9%
- Balanced accuracy: 73.9%

### Interpretation
Despite differences in cohort composition and sequencing pipelines, the diagnostic signal generalizes across datasets. Threshold-dependent trade-offs suggest distinct use cases, ranging from **screening-oriented detection** to **confirmatory diagnostics**.

---

## Therapeutic Target Prioritization

Genes were ranked using a composite **TargetScore** incorporating:

- Absolute log2 fold change
- Statistical confidence (−log10 adjusted p-value)
- Expression magnitude
- Direction of dysregulation

**Flagged cancer-relevant targets include:**
- EGFR
- ERBB2
- MET
- RET
- NTRK2 / NTRK3
- TIGIT
- LAG3

> These candidates represent **hypothesis-generating targets**, not validated therapeutic recommendations.

---

## Reproducibility

### Public-safe visibility
- Raw GEO files are intentionally **not committed** to the repository.
- `metadata/sample_map_all.tsv`, scripts, and published outputs are tracked for reproducibility.
- CPTAC validation is dependency-driven (`cptac` package), and is not committed as raw data.

### Data setup (required before running)
Use the local manifest and downloader:

```bash
cd /path/to/tcga-luad-dx-targets-showcase
bash scripts/00_download_gse62944.sh
```

The downloader prints exact source links and tells you whether each file is present.
It does not auto-download by default. To attempt automatic retrieval, use:

```bash
bash scripts/00_download_gse62944.sh --download
```

### Optional checksum verification
If `data_sources.tsv` includes SHA256 values, this pipeline verifies them when present.

```bash
# Generate a standard hash-check file (when hash values are present)
awk -F '\t' 'BEGIN {OFS=" "} /^[^#]/ && NF >= 6 && $6 != "" {print $6, $4}' data_sources.tsv > data_sources.sha256
sha256sum -c data_sources.sha256
```

If hashes are missing in the manifest, the run will continue with a warning.

### Quick run (recommended)
```bash
conda activate luad_tcga
cd /path/to/tcga-luad-dx-targets-showcase
export TCGA_LUAD_ROOT="$PWD"
bash scripts/99_run_all.sh
```

`scripts/99_run_all.sh` will:
- Run required preprocessing scripts from the repo root.
- Attempt to download/build missing local inputs when possible.
- Fail fast with clear messages when required files are still missing.

If you don’t want to set `TCGA_LUAD_ROOT`, the scripts attempt to infer the project root from their location.

For manual step-by-step execution:
```bash
Rscript scripts/00_download_gse62944.sh
Rscript scripts/01_build_coldata.R
Rscript scripts/02_subset_counts.R
Rscript scripts/03_deseq2_LUAD.R
python scripts/06_panel_rank_stability.py
Rscript scripts/05_target_prioritization.R
Rscript scripts/07_make_summary_tables.R
Rscript scripts/08_boxplots_key_genes.R
python scripts/09_external_validation_cptac_luad.py
```

### Expected outputs

Tables

deseq2_all.tsv

deseq2_sig.tsv

diagnostic_panel_rankstable_top20.tsv

classifier_metrics_rankstable.tsv

target_shortlist_top50.tsv

targets_flagged_actionable.tsv

external_metrics_cptac_luad.tsv

Figures

ma_plot.png

pca_vst.png

volcano.png

heatmap_top50.png

roc_curve.png

roc_curve_external_cptac_luad.png

boxplot_*.png

### Public posting checklist
- Ensure no personal secrets, API tokens, or credentials appear in environment files.
- Ensure no untracked large raw/processed inputs are committed.
- Confirm your branch list includes only intended work branches before sharing.
- Push only intended branches/tags to GitHub.

Limitations

Bulk RNA-seq cannot resolve cell-type–specific effects

No survival or outcome modeling included

Functional validation not performed

Next Steps

Molecular subtyping within LUAD

Survival analysis (Kaplan–Meier, Cox models)

Pathway enrichment (GO / KEGG / Reactome)

Multi-omics integration (mutation, CNV)

Extension to additional cancers (brain, breast, prostate)

Author

Kevin Arredondo
B.S. Molecular Cell Biology
Bioinformatics & Translational Cancer Research
