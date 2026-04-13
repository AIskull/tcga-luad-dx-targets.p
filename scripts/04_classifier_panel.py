#!/usr/bin/env python3
"""
Step 8) Diagnostic classifier + minimal gene panel (LUAD Tumor vs Normal)

Inputs:
  - data/processed/luad_counts.tsv
  - metadata/luad_coldata.tsv
  - results/tables/deseq2_sig.tsv

Outputs:
  - results/tables/diagnostic_panel_top20.tsv
  - results/tables/classifier_metrics.tsv
  - results/tables/confusion_matrix.tsv
  - results/figures/roc_curve.png
"""

import os
import numpy as np
import pandas as pd

from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.metrics import roc_auc_score, roc_curve, confusion_matrix
import matplotlib.pyplot as plt


# -----------------------------
# Config
# -----------------------------
COUNTS_PATH = "data/processed/luad_counts.tsv"
COLDATA_PATH = "metadata/luad_coldata.tsv"
DESEQ_SIG_PATH = "results/tables/deseq2_sig.tsv"

OUT_PANEL = "results/tables/diagnostic_panel_top20.tsv"
OUT_METRICS = "results/tables/classifier_metrics.tsv"
OUT_CM = "results/tables/confusion_matrix.tsv"
OUT_ROC = "results/figures/roc_curve.png"

TOP_N_FEATURES = 200          # pick top N DE genes (by padj then |log2FC|)
N_SPLITS = 5                  # CV folds
RANDOM_STATE = 42
THRESHOLD = 0.5               # for confusion matrix


# -----------------------------
# Helpers
# -----------------------------
def must_exist(path: str) -> None:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing required file: {path}")


def main():
    # Ensure output dirs exist
    os.makedirs("results/tables", exist_ok=True)
    os.makedirs("results/figures", exist_ok=True)

    # Check inputs exist
    must_exist(COUNTS_PATH)
    must_exist(COLDATA_PATH)
    must_exist(DESEQ_SIG_PATH)

    print(f"Reading counts: {COUNTS_PATH}")
    counts = pd.read_csv(COUNTS_PATH, sep="\t")
    gene_col = counts.columns[0]
    counts = counts.set_index(gene_col)  # rows=genes, cols=samples

    print(f"Reading coldata: {COLDATA_PATH}")
    coldata = pd.read_csv(COLDATA_PATH, sep="\t")
    if "sample_id" not in coldata.columns or "condition" not in coldata.columns:
        raise ValueError("coldata must contain columns: sample_id, condition")
    coldata = coldata.set_index("sample_id")

    # y: Tumor=1, Normal=0
    y = (coldata["condition"].astype(str) == "Tumor").astype(int)

    # Align samples robustly (intersection + same order)
    counts_samples = list(counts.columns)
    coldata_samples = list(coldata.index)

    common = [s for s in coldata_samples if s in counts_samples]
    print(f"Samples in coldata: {len(coldata_samples)}")
    print(f"Samples in counts: {len(counts_samples)}")
    print(f"Common samples: {len(common)}")

    if len(common) < 10:
        raise ValueError("Too few common samples after alignment — check sample IDs.")

    # X: samples x genes
    X = counts[common].T
    y = y.loc[common]

    # log2(count+1)
    X = np.log2(X + 1)

    # Load DESeq2 significant genes and pick features
    print(f"Reading DESeq2 sig: {DESEQ_SIG_PATH}")
    sig = pd.read_csv(DESEQ_SIG_PATH, sep="\t")

    needed_cols = {"gene_id", "padj", "log2FoldChange"}
    if not needed_cols.issubset(set(sig.columns)):
        raise ValueError(f"deseq2_sig.tsv must contain columns: {sorted(needed_cols)}")

    sig = sig.dropna(subset=["padj", "log2FoldChange"]).copy()
    # sort by padj, then bigger absolute fold-change
    sig["absLFC"] = sig["log2FoldChange"].abs()
    sig = sig.sort_values(["padj", "absLFC"], ascending=[True, False])

    top_genes = sig["gene_id"].head(TOP_N_FEATURES).astype(str).tolist()

    # Keep only genes that actually exist in X
    top_genes = [g for g in top_genes if g in X.columns]
    if len(top_genes) < 10:
        raise ValueError("After filtering to genes present in counts, too few features remain.")

    print(f"Using TOP_N_FEATURES requested: {TOP_N_FEATURES}")
    print(f"Features present in matrix: {len(top_genes)}")

    X_sel = X[top_genes]

    # Model pipeline
    pipe = Pipeline([
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(max_iter=5000, solver="liblinear"))
    ])

    cv = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=RANDOM_STATE)

    aucs = []
    all_probs = np.zeros(len(y), dtype=float)

    print("Running cross-validation...")
    for fold, (train_idx, test_idx) in enumerate(cv.split(X_sel, y), start=1):
        X_train, X_test = X_sel.iloc[train_idx], X_sel.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

        pipe.fit(X_train, y_train)
        probs = pipe.predict_proba(X_test)[:, 1]
        all_probs[test_idx] = probs

        fold_auc = roc_auc_score(y_test, probs)
        aucs.append(fold_auc)
        print(f"  Fold {fold}/{N_SPLITS} AUC: {fold_auc:.4f}")

    mean_auc = float(np.mean(aucs))
    std_auc = float(np.std(aucs))

    # Confusion matrix using CV probabilities threshold
    y_pred = (all_probs >= THRESHOLD).astype(int)
    cm = confusion_matrix(y, y_pred)
    cm_df = pd.DataFrame(
        cm,
        index=["True_Normal(0)", "True_Tumor(1)"],
        columns=["Pred_Normal(0)", "Pred_Tumor(1)"]
    )
    cm_df.to_csv(OUT_CM, sep="\t")

    # Fit on ALL data to get coefficients for a "panel"
    print("Fitting final model on all data for coefficients...")
    pipe.fit(X_sel, y)
    coef = pipe.named_steps["clf"].coef_[0]

    coef_df = (
        pd.DataFrame({"gene": top_genes, "coef": coef})
        .sort_values("coef", ascending=False)
        .reset_index(drop=True)
    )

    # Panel = top 10 positive + top 10 negative
    panel = pd.concat([coef_df.head(10), coef_df.tail(10)], axis=0)
    panel.to_csv(OUT_PANEL, sep="\t", index=False)

    # Metrics table
    metrics_df = pd.DataFrame({
        "metric": ["roc_auc_cv_mean", "roc_auc_cv_std", "n_samples", "n_features_used", "threshold_for_cm"],
        "value": [mean_auc, std_auc, int(len(y)), int(len(top_genes)), float(THRESHOLD)]
    })
    metrics_df.to_csv(OUT_METRICS, sep="\t", index=False)

    # ROC curve plot (CV probs)
    fpr, tpr, _ = roc_curve(y, all_probs)
    plt.figure()
    plt.plot(fpr, tpr)
    plt.plot([0, 1], [0, 1], linestyle="--")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"LUAD Tumor vs Normal ROC (CV mean AUC={mean_auc:.3f})")
    plt.savefig(OUT_ROC, dpi=150, bbox_inches="tight")

    print(f"Wrote: {OUT_PANEL}")
    print(f"Wrote: {OUT_METRICS}")
    print(f"Wrote: {OUT_CM}")
    print(f"Saved: {OUT_ROC}")
    print("Done.")


if __name__ == "__main__":
    main()
