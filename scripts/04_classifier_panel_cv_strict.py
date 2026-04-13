#!/usr/bin/env python3
"""
Leakage-free classifier (strict CV) with fold-specific feature selection.

What this does (defensible / grad-ready):
- Loads counts + labels
- Uses log2(count+1) expression
- For each CV fold:
    * Selects top N genes using ONLY the training split (Welch t-test ranking)
    * Trains logistic regression on training split
    * Evaluates on held-out split
- Aggregates out-of-fold probabilities -> ROC/AUC + confusion matrix
- Fits a final model on ALL samples using a "global" top-N gene set (for a final panel)

Writes:
  results/tables/classifier_metrics_strict.tsv
  results/tables/confusion_matrix_strict.tsv
  results/tables/panel_stability_strict.tsv
  results/tables/diagnostic_panel_top20_strict.tsv
  results/tables/global_feature_ranking_strict.tsv
  results/figures/roc_curve_strict.png
"""

import os
import sys
import numpy as np
import pandas as pd

from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.metrics import roc_auc_score, roc_curve, confusion_matrix

import matplotlib
matplotlib.use("Agg")  # headless-safe
import matplotlib.pyplot as plt


# ---------------- Config ----------------
COUNTS_PATH = "data/processed/luad_counts.tsv"
COLDATA_PATH = "metadata/luad_coldata.tsv"

OUT_TABLES = "results/tables"
OUT_FIGS = "results/figures"

TOP_N_FEATURES = 200
N_SPLITS = 5
RANDOM_STATE = 42

PANEL_POS = 10
PANEL_NEG = 10

THRESHOLD = 0.5


# ---------------- Helpers ----------------
def log(msg: str):
    print(msg, flush=True)

def assert_file(path: str):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing required file: {path}")

def ensure_dirs():
    os.makedirs(OUT_TABLES, exist_ok=True)
    os.makedirs(OUT_FIGS, exist_ok=True)

def safe_log2p1(df: pd.DataFrame) -> pd.DataFrame:
    # counts should be non-negative; handle safely
    if (df.values < 0).any():
        raise ValueError("Counts contain negative values. Counts must be >= 0.")
    return np.log2(df + 1.0)

def load_data():
    assert_file(COUNTS_PATH)
    assert_file(COLDATA_PATH)

    log(f"Reading counts: {COUNTS_PATH}")
    counts = pd.read_csv(COUNTS_PATH, sep="\t")

    if counts.shape[1] < 3:
        raise ValueError("Counts matrix has too few columns (need gene_id + >=2 samples).")

    gene_col = counts.columns[0]
    counts = counts.set_index(gene_col)

    # sanity: duplicate gene ids?
    if counts.index.duplicated().any():
        dup_n = int(counts.index.duplicated().sum())
        raise ValueError(f"Counts matrix has duplicated gene IDs ({dup_n}). Fix before modeling.")

    log(f"Reading coldata: {COLDATA_PATH}")
    coldata = pd.read_csv(COLDATA_PATH, sep="\t")
    if "sample_id" not in coldata.columns or "condition" not in coldata.columns:
        raise ValueError("coldata must contain columns: sample_id, condition")

    coldata = coldata.set_index("sample_id")

    log(f"Samples in coldata: {coldata.shape[0]}")
    log(f"Samples in counts: {counts.shape[1]}")

    common = coldata.index.intersection(counts.columns)
    log(f"Common samples: {len(common)}")

    if len(common) < 2:
        raise ValueError("Too few overlapping sample IDs between counts and coldata.")

    # keep consistent order
    coldata = coldata.loc[common].copy()
    counts = counts[common].copy()

    if "condition" not in coldata.columns:
        raise ValueError("coldata missing 'condition' column after indexing.")

    # binary labels
    y = (coldata["condition"].astype(str) == "Tumor").astype(int)

    # X: samples x genes
    X = counts.T
    X = safe_log2p1(X)

    # final sanity
    if X.isna().any().any():
        raise ValueError("Expression matrix contains NA values after log transform.")

    return X, y, coldata

def build_model():
    # deterministic + stable
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

def welch_t_rank(X_train: pd.DataFrame, y_train: pd.Series) -> pd.DataFrame:
    """
    Rank genes by training-only Welch t-test (fast, no scipy needed).
    We use p-approx only for ranking (NOT significance claims).
    """
    X0 = X_train[y_train == 0]
    X1 = X_train[y_train == 1]

    n0 = int(X0.shape[0])
    n1 = int(X1.shape[0])
    if n0 < 2 or n1 < 2:
        raise ValueError(f"Need >=2 samples per class in training fold (got n0={n0}, n1={n1}).")

    m0 = X0.mean(axis=0).to_numpy()
    m1 = X1.mean(axis=0).to_numpy()
    v0 = X0.var(axis=0, ddof=1).to_numpy()
    v1 = X1.var(axis=0, ddof=1).to_numpy()

    denom = np.sqrt(v0 / n0 + v1 / n1)
    denom[denom == 0] = np.nan
    tstat = (m1 - m0) / denom

    # Normal approx for two-sided p-value: p = 2*(1 - Phi(|t|))
    # Phi(x) = 0.5*(1 + erf(x/sqrt(2)))
    from math import erf, sqrt
    abs_t = np.abs(tstat)
    phi = 0.5 * (1.0 + np.vectorize(erf)(abs_t / sqrt(2.0)))
    p_approx = 2.0 * (1.0 - phi)

    genes = X_train.columns.to_numpy()

    df = pd.DataFrame({
        "gene": genes,
        "tstat": tstat,
        "p_approx": p_approx,
        "mean_normal": m0,
        "mean_tumor": m1,
        "mean_diff": (m1 - m0),
    })

    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=["p_approx"]).sort_values(["p_approx"], ascending=True)

    return df

def select_top_n(rank_df: pd.DataFrame, top_n: int):
    return rank_df["gene"].head(top_n).tolist()

def save_confusion(cm: np.ndarray, path: str):
    cm_df = pd.DataFrame(
        cm,
        index=["True_Normal(0)", "True_Tumor(1)"],
        columns=["Pred_Normal(0)", "Pred_Tumor(1)"]
    )
    cm_df.to_csv(path, sep="\t")

def save_roc(y: pd.Series, probs: np.ndarray, path: str, title: str):
    fpr, tpr, _ = roc_curve(y, probs)
    plt.figure()
    plt.plot(fpr, tpr)
    plt.plot([0, 1], [0, 1], linestyle="--")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(title)
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()


# ---------------- Main ----------------
def main():
    ensure_dirs()

    X, y, _ = load_data()

    log(f"Using TOP_N_FEATURES requested: {TOP_N_FEATURES}")
    log("Running leakage-free cross-validation (feature selection inside folds)...")

    cv = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=RANDOM_STATE)

    aucs = []
    all_probs = np.zeros(len(y), dtype=float)

    # stability: counts of how often each gene was selected
    selection_counts = {}

    for fold, (train_idx, test_idx) in enumerate(cv.split(X, y), start=1):
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

        rank_df = welch_t_rank(X_train, y_train)
        fold_genes = select_top_n(rank_df, TOP_N_FEATURES)

        for g in fold_genes:
            selection_counts[g] = selection_counts.get(g, 0) + 1

        model = build_model()
        model.fit(X_train[fold_genes], y_train)

        probs = model.predict_proba(X_test[fold_genes])[:, 1]
        all_probs[test_idx] = probs

        fold_auc = roc_auc_score(y_test, probs)
        aucs.append(fold_auc)
        log(f"  Fold {fold}/{N_SPLITS} AUC: {fold_auc:.4f} (features={len(fold_genes)})")

    mean_auc = float(np.mean(aucs))
    std_auc = float(np.std(aucs))

    # out-of-fold confusion matrix
    preds = (all_probs >= THRESHOLD).astype(int)
    cm = confusion_matrix(y, preds, labels=[0, 1])

    # --- write metrics
    metrics_path = os.path.join(OUT_TABLES, "classifier_metrics_strict.tsv")
    pd.DataFrame({
        "metric": [
            "roc_auc_cv_mean",
            "roc_auc_cv_std",
            "n_samples",
            "n_splits",
            "top_n_features_per_fold",
            "threshold_for_cm",
        ],
        "value": [
            mean_auc,
            std_auc,
            float(len(y)),
            float(N_SPLITS),
            float(TOP_N_FEATURES),
            float(THRESHOLD),
        ]
    }).to_csv(metrics_path, sep="\t", index=False)
    log(f"Wrote: {metrics_path}")

    cm_path = os.path.join(OUT_TABLES, "confusion_matrix_strict.tsv")
    save_confusion(cm, cm_path)
    log(f"Wrote: {cm_path}")

    roc_path = os.path.join(OUT_FIGS, "roc_curve_strict.png")
    save_roc(
        y,
        all_probs,
        roc_path,
        title=f"LUAD Tumor vs Normal ROC (Leakage-free CV AUC={mean_auc:.3f}±{std_auc:.3f})"
    )
    log(f"Saved: {roc_path}")

    # --- stability table
    stability_path = os.path.join(OUT_TABLES, "panel_stability_strict.tsv")
    stab = pd.DataFrame({
        "gene": list(selection_counts.keys()),
        "folds_selected": list(selection_counts.values()),
    })
    stab["selection_frequency"] = stab["folds_selected"] / float(N_SPLITS)
    stab = stab.sort_values(["folds_selected", "gene"], ascending=[False, True])
    stab.to_csv(stability_path, sep="\t", index=False)
    log(f"Wrote: {stability_path}")

    # --- final "global" top-N + coefficients (panel)
    log("Fitting final model on ALL data for coefficients (global top-N selection)...")
    global_rank = welch_t_rank(X, y).copy()
    global_rank_path = os.path.join(OUT_TABLES, "global_feature_ranking_strict.tsv")
    global_rank.to_csv(global_rank_path, sep="\t", index=False)
    log(f"Wrote: {global_rank_path}")

    global_genes = select_top_n(global_rank, TOP_N_FEATURES)

    model_final = build_model()
    model_final.fit(X[global_genes], y)

    coef = model_final.named_steps["clf"].coef_[0]
    coef_df = pd.DataFrame({"gene": global_genes, "coef": coef}).sort_values("coef", ascending=False)

    panel = pd.concat([coef_df.head(PANEL_POS), coef_df.tail(PANEL_NEG)], axis=0)
    panel_path = os.path.join(OUT_TABLES, "diagnostic_panel_top20_strict.tsv")
    panel.to_csv(panel_path, sep="\t", index=False)
    log(f"Wrote: {panel_path}")

    log("Done.")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        log(f"ERROR: {e}")
        sys.exit(1)

