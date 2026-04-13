# LUAD RNA-seq Differential Expression and External Validation

## Overview

This repository documents an independent computational biology project focused on lung adenocarcinoma (LUAD) tumor-vs-normal RNA-seq analysis using public datasets.

The workflow covers data acquisition, metadata construction, count processing, differential expression analysis, classifier development, and external validation. The goal is to show a reproducible end-to-end analysis in R and Python using public LUAD transcriptomic resources.

> Research use only. This repository is intended for research communication and reproducibility, not for clinical diagnosis, treatment selection, or patient-care decisions.

---

## Skills Demonstrated

- RNA-seq data processing and cohort filtering
- Metadata curation for tumor-vs-normal analysis
- Differential expression analysis with DESeq2
- Logistic regression classifier development and evaluation
- Cross-cohort validation with an independent LUAD dataset
- Reproducible workflow development in R and Python
- Figure and summary-table generation for scientific reporting

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
- **Access path in this repo:** transcriptomics loaded through the Python `cptac` package for external validation only; the workflow does not redistribute CPTAC raw files or clinical tables

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
Despite differences in cohort composition and sequencing pipelines, the classifier retains measurable cross-cohort signal. The reduced external performance relative to the TCGA training cohort highlights the difficulty of transferring transcriptomic models across datasets with different technical and biological characteristics.

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

## Limitations

- The external validation analysis is transcriptomics-only and does not incorporate orthogonal molecular or clinical features.
- Performance differs between the TCGA training cohort and the CPTAC external cohort, so results should be interpreted as research findings rather than deployment-ready models.
- Bulk RNA-seq does not resolve cell-type-specific effects.
- Functional validation of highlighted genes was not performed in this repository.

---

## Reproducibility

### Public-safe visibility
- Raw GEO files are intentionally **not committed** to the repository.
- `metadata/sample_map_all.tsv`, scripts, and published outputs are tracked for reproducibility.
- CPTAC validation is dependency-driven (`cptac` package), and is not committed as raw data.

### Data access and attribution
This repository is a public-facing analysis and reproducibility snapshot. Users should obtain source data from the official providers and comply with upstream access, attribution, and reuse terms.

- **GEO / GSE62944:** training inputs in this project are derived from public GEO accession `GSE62944`. GEO states that NCBI places no restrictions on use or distribution of GEO data, but submitters may assert patent, copyright, or other intellectual property rights in submitted data.
- **TCGA / GDC:** where open-access TCGA or GDC resources are involved, reuse requires proper accreditation. Controlled-access resources require separate authorization and are not redistributed in this repository.
- **CPTAC LUAD:** the external validation workflow uses transcriptomics accessed through the Python `cptac` package. This repository does not redistribute CPTAC raw files, CPTAC clinical tables, controlled-access data, or access credentials.

When reusing this repository, please cite or acknowledge the original data resources and relevant upstream publications, including GEO `GSE62944`, TCGA / GDC, and CPTAC.

### Data setup (required before running)
Use the local manifest and downloader:

```bash
cd /path/to/tcga-luad-dx-targets.p
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
cd /path/to/tcga-luad-dx-targets.p
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

## Public Posting Checklist
- Ensure no personal secrets, API tokens, or credentials appear in environment files.
- Ensure no untracked large raw/processed inputs are committed.
- Confirm your branch list includes only intended work branches before sharing.
- Push only intended branches/tags to GitHub.

## Next Steps

Molecular subtyping within LUAD

Survival analysis (Kaplan–Meier, Cox models)

Pathway enrichment (GO / KEGG / Reactome)

Multi-omics integration (mutation, CNV)

Extension to additional cancers (brain, breast, prostate)

Author

Kevin Arredondo
B.S. Molecular Cell Biology
