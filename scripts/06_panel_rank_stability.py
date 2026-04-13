#!/usr/bin/env python3
"""
Rank-stable diagnostic panel builder (heuristic panel selection).

Pipeline:
- Load counts + coldata
- log2(count+1)
- Use DESeq2 significant genes (results/tables/deseq2_sig.tsv) to define a TOP_N feature pool
- Cross-validate logistic regression on that pool
- Inside each fold: rank genes by |coef| and store their ranks
- Aggregate mean rank across folds -> pick FINAL_PANEL genes (most rank-stable)
- Fit final model on ALL samples using FINAL_PANEL genes -> write panel coefficients
- Save out-of-fold ROC, CV AUC stats, and confusion matrix

Outputs:
  results/tables/diagnostic_panel_rankstable_top20.tsv
  results/tables/classifier_metrics_rankstable.tsv
  results/tables/confusion_matrix_rankstable.tsv
  results/tables/rank_stability_table_rankstable.tsv
  results/figures/roc_curve_rankstable.png

Important (README honesty):
- This script uses DESeq2_sig.tsv created on the full dataset to choose the TOP_N pool,
  so it is NOT leakage-free evaluation. Use 04_classifier_panel_cv_strict.py for strict CV.
"""

import os
import sys
import numpy as np
import pandas as pd
from collections import defaultdict

from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.metrics import roc_auc_score, roc_curve, confusion_matrix

import matplotlib
matplotlib.use("Agg")  # headless-safe (prevents Wayland/Qt errors)
import matplotlib.pyplot as plt


# ---------------- Config ----------------
TOP_N = 200          # feature pool size from DESeq2_sig (top by padj)
FINAL_PANEL = 20     # final panel size
N_SPLITS = 5
RANDOM_STATE = 42
THRESHOLD = 0.5

COUNTS_PATH = "data/processed/luad_counts.tsv"
COLDATA_PATH = "metadata/luad_coldata.tsv"
DESEQ2_SIG_PATH = "results/tables/deseq2_sig.tsv"

OUT_TABLES = "results/tables"
OUT_FIGS = "results/figures"

OUT_PANEL = f"{OUT_TABLES}/diagnostic_panel_rankstable_top{FINAL_PANEL}.tsv"
OUT_METRICS = f"{OUT_TABLES}/classifier_metrics_rankstable.tsv"
OUT_CM = f"{OUT_TABLES}/confusion_matrix_rankstable.tsv"
OUT_RANKTAB = f"{OUT_TABLES}/rank_stability_table_rankstable.tsv"
OUT_ROC = f"{OUT_FIGS}/roc_curve_rankstable.png"


# ---------------- Helpers ----------------
def log(msg: str):
    print(msg, flush=True)

def ensure_dirs():
    os.makedirs(OUT_TABLES, exist_ok=True)
    os.makedirs(OUT_FIGS, exist_ok=True)

def assert_file(path: str):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing required file: {path}")

def build_model():
    return Pipeline([
        ("scaler", StandardScaler(with_mean=True, with_std=True)),
        ("clf", LogisticRegression(
            max_iter=5000,
            solver="lbfgs",
            penalty="l2",
            class_weight="balanced",
            random_state=RANDOM_STATE
        ))
    ])

def load_counts_coldata():
    assert_file(COUNTS_PATH)
    assert_file(COLDATA_PATH)

    log(f"Reading counts: {COUNTS_PATH}")
    counts = pd.read_csv(COUNTS_PATH, sep="\t")
    gene_col = counts.columns[0]
    counts = counts.set_index(gene_col)

    if counts.index.duplicated().any():
        dup = int(counts.index.duplicated().sum())
        raise ValueError(f"Counts has duplicated gene IDs ({dup}). Fix before running.")

    log(f"Reading coldata: {COLDATA_PATH}")
    coldata = pd.read_csv(COLDATA_PATH, sep="\t")
    if "sample_id" not in coldata.columns or "condition" not in coldata.columns:
        raise ValueError("coldata must contain columns: sample_id, condition")

    coldata = coldata.set_index("sample_id")

    # align sample IDs
    common = coldata.index.intersection(counts.columns)
    log(f"Samples in coldata: {coldata.shape[0]}")
    log(f"Samples in counts: {counts.shape[1]}")
    log(f"Common samples: {len(common)}")
    if len(common) < 2:
        raise ValueError("Too few overlapping sample IDs between counts and coldata.")

    coldata = coldata.loc[common].copy()
    counts = counts[common].copy()

    y = (coldata["condition"].astype(str) == "Tumor").astype(int)
    X = counts.T

    # log2(count+1)
    if (X.values < 0).any():
        raise ValueError("Counts contain negative values. Counts must be >=0.")
    X = np.log2(X + 1.0)

    return X, y, coldata

def load_top_genes_from_deseq2_sig(top_n: int):
    assert_file(DESEQ2_SIG_PATH)
    log(f"Reading DESeq2 sig: {DESEQ2_SIG_PATH}")

    sig = pd.read_csv(DESEQ2_SIG_PATH, sep="\t")
    if "gene_id" not in sig.columns or "padj" not in sig.columns:
        raise ValueError("deseq2_sig.tsv must include columns: gene_id, padj")

    sig = sig.dropna(subset=["padj"]).copy()
    # stable ordering: primarily padj, then absolute effect
    if "log2FoldChange" in sig.columns:
        sig["absLFC"] = sig["log2FoldChange"].abs()
        sig = sig.sort_values(["padj", "absLFC"], ascending=[True, False])
    else:
        sig = sig.sort_values(["padj"], ascending=[True])

    genes = sig["gene_id"].astype(str).head(top_n).tolist()
    if len(genes) == 0:
        raise ValueError("No genes available in deseq2_sig.tsv after filtering padj.")

    return genes

def save_roc(y, probs, auc_mean, path):
    fpr, tpr, _ = roc_curve(y, probs)
    plt.figure()
    plt.plot(fpr, tpr)
    plt.plot([0, 1], [0, 1], linestyle="--")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"Rank-stable Panel ROC (CV mean AUC={auc_mean:.3f})")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()


# ---------------- Main ----------------
def main():
    ensure_dirs()

    X_all, y, _ = load_counts_coldata()
    top_genes = load_top_genes_from_deseq2_sig(TOP_N)

    # Keep only genes present in expression matrix
    present = [g for g in top_genes if g in X_all.columns]
    missing = [g for g in top_genes if g not in X_all.columns]
    log(f"TOP_N requested: {TOP_N}")
    log(f"Features present in matrix: {len(present)}")
    if len(missing) > 0:
        log(f"Note: {len(missing)} genes from DESeq2_sig not found in counts matrix (ignored).")

    if len(present) < FINAL_PANEL:
        raise ValueError(f"Not enough features present ({len(present)}) to build FINAL_PANEL={FINAL_PANEL}.")

    X = X_all[present].copy()

    cv = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=RANDOM_STATE)

    rank_store = defaultdict(list)
    probs_oof = np.zeros(len(y), dtype=float)
    aucs = []

    log("Running rank-stability CV...")
    for fold, (tr, te) in enumerate(cv.split(X, y), start=1):
        model = build_model()
        model.fit(X.iloc[tr], y.iloc[tr])

        probs = model.predict_proba(X.iloc[te])[:, 1]
        probs_oof[te] = probs

        auc = roc_auc_score(y.iloc[te], probs)
        aucs.append(auc)
        log(f"  Fold {fold}/{N_SPLITS} AUC: {auc:.4f}")

        coef = model.named_steps["clf"].coef_[0]
        ranked = pd.Series(np.abs(coef), index=X.columns).sort_values(ascending=False)

        for r, gene in enumerate(ranked.index, start=1):
            rank_store[gene].append(r)

    # Aggregate ranks across folds
    rank_df = pd.DataFrame({
        "gene": list(rank_store.keys()),
        "mean_rank": [float(np.mean(v)) for v in rank_store.values()],
        "std_rank": [float(np.std(v)) for v in rank_store.values()],
        "folds_seen": [len(v) for v in rank_store.values()],
    }).sort_values(["mean_rank", "std_rank"], ascending=[True, True])

    rank_df.to_csv(OUT_RANKTAB, sep="\t", index=False)
    log(f"Wrote: {OUT_RANKTAB}")

    final_genes = rank_df.head(FINAL_PANEL)["gene"].tolist()

    # Fit final model on ALL data using fixed panel
    X_final = X_all[final_genes]
    final_model = build_model()
    final_model.fit(X_final, y)

    coef = final_model.named_steps["clf"].coef_[0]
    panel = pd.DataFrame({"gene": final_genes, "coef": coef}).sort_values("coef", ascending=False)
    panel.to_csv(OUT_PANEL, sep="\t", index=False)
    log(f"Wrote: {OUT_PANEL}")

    # Metrics + confusion matrix (OOF)
    auc_mean = float(np.mean(aucs))
    auc_std = float(np.std(aucs))

    pd.DataFrame({
        "metric": ["roc_auc_cv_mean", "roc_auc_cv_std", "n_samples", "n_features_pool", "n_features_final", "threshold_for_cm"],
        "value": [auc_mean, auc_std, float(len(y)), float(len(present)), float(FINAL_PANEL), float(THRESHOLD)]
    }).to_csv(OUT_METRICS, sep="\t", index=False)
    log(f"Wrote: {OUT_METRICS}")

    preds = (probs_oof >= THRESHOLD).astype(int)
    cm = confusion_matrix(y, preds, labels=[0, 1])
    cm_df = pd.DataFrame(cm, index=["True_Normal", "True_Tumor"], columns=["Pred_Normal", "Pred_Tumor"])
    cm_df.to_csv(OUT_CM, sep="\t")
    log(f"Wrote: {OUT_CM}")

    save_roc(y, probs_oof, auc_mean, OUT_ROC)
    log(f"Saved: {OUT_ROC}")

    log("Done. Rank-stable diagnostic panel created.")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        log(f"ERROR: {e}")
        sys.exit(1)

