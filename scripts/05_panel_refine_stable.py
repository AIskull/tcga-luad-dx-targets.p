#!/usr/bin/env python3
"""
Refine a STABLE minimal diagnostic panel from leakage-free CV results.

Inputs:
- data/processed/luad_counts.tsv
- metadata/luad_coldata.tsv
- results/tables/panel_stability_strict.tsv  (from 04_classifier_panel_cv_strict.py)

What it does:
1) Filters genes by fold-stability (selected in >= MIN_FOLDS of CV folds)
2) Fits a model using ONLY those stable genes (still strict CV evaluation)
3) Builds a top panel (top positive + top negative coefficients) from final fit
4) Saves stable metrics + ROC

Outputs:
- results/tables/diagnostic_panel_stable_top20.tsv
- results/tables/classifier_metrics_stable.tsv
- results/tables/confusion_matrix_stable.tsv
- results/figures/roc_curve_stable.png
"""

import os
import numpy as np
import pandas as pd

from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.metrics import roc_auc_score, roc_curve, confusion_matrix
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------- Config ----------
COUNTS_PATH = "data/processed/luad_counts.tsv"
COLDATA_PATH = "metadata/luad_coldata.tsv"
STABILITY_PATH = "results/tables/panel_stability_strict.tsv"

OUT_DIR_TABLES = "results/tables"
OUT_DIR_FIGS = "results/figures"

N_SPLITS = 5
RANDOM_STATE = 42

# Stability filter: require a gene be selected in at least this many folds
MIN_FOLDS = 2   # 2/5 folds is the best available stability here

# After stability filter, cap features to this many (most stable first)
MAX_STABLE_FEATURES = 200

# Panel size
PANEL_POS = 10
PANEL_NEG = 10

THRESHOLD = 0.5

# ---------- Helpers ----------
def ensure_dirs():
    os.makedirs(OUT_DIR_TABLES, exist_ok=True)
    os.makedirs(OUT_DIR_FIGS, exist_ok=True)

def log(msg: str):
    print(msg, flush=True)

def build_model():
    return Pipeline([
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(max_iter=5000, solver="lbfgs"))
    ])

def load_counts_and_labels():
    log(f"Reading counts: {COUNTS_PATH}")
    counts = pd.read_csv(COUNTS_PATH, sep="\t")
    gene_col = counts.columns[0]
    counts = counts.set_index(gene_col)

    log(f"Reading coldata: {COLDATA_PATH}")
    coldata = pd.read_csv(COLDATA_PATH, sep="\t").set_index("sample_id")

    common = coldata.index.intersection(counts.columns)
    log(f"Common samples: {len(common)}")
    if len(common) == 0:
        raise ValueError("No overlapping sample IDs between counts and coldata.")

    counts = counts[common]
    coldata = coldata.loc[common]

    y = (coldata["condition"] == "Tumor").astype(int)

    X = counts.T
    X = np.log2(X + 1)  # log transform

    return X, y, coldata

def load_stable_genes():
    log(f"Reading stability: {STABILITY_PATH}")
    stab = pd.read_csv(STABILITY_PATH, sep="\t")

    # Filter by stability threshold
    stable = stab[stab["folds_selected"] >= MIN_FOLDS].copy()
    stable = stable.sort_values(["folds_selected", "gene"], ascending=[False, True])

    genes = stable["gene"].tolist()
    if len(genes) == 0:
        raise ValueError(
            f"No genes passed stability filter (MIN_FOLDS={MIN_FOLDS}). "
            "Try lowering MIN_FOLDS to 3."
        )

    # Cap number of stable features if huge
    genes = genes[:MAX_STABLE_FEATURES]
    log(f"Stable genes kept: {len(genes)} (MIN_FOLDS={MIN_FOLDS}, cap={MAX_STABLE_FEATURES})")
    return genes

def main():
    ensure_dirs()

    X, y, _ = load_counts_and_labels()
    stable_genes = load_stable_genes()

    # Keep only stable genes present in matrix
    stable_genes = [g for g in stable_genes if g in X.columns]
    if len(stable_genes) == 0:
        raise ValueError("None of the stable genes were found in X columns.")

    Xs = X[stable_genes]

    log("Running strict CV using ONLY stable genes...")
    cv = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=RANDOM_STATE)

    aucs = []
    all_probs = np.zeros(len(y), dtype=float)

    for fold, (train_idx, test_idx) in enumerate(cv.split(Xs, y), start=1):
        X_train, X_test = Xs.iloc[train_idx], Xs.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

        model = build_model()
        model.fit(X_train, y_train)
        probs = model.predict_proba(X_test)[:, 1]
        all_probs[test_idx] = probs

        fold_auc = roc_auc_score(y_test, probs)
        aucs.append(fold_auc)
        log(f"  Fold {fold}/{N_SPLITS} AUC: {fold_auc:.4f} (features={Xs.shape[1]})")

    mean_auc = float(np.mean(aucs))
    std_auc = float(np.std(aucs))

    # Confusion matrix (out-of-fold)
    preds = (all_probs >= THRESHOLD).astype(int)
    cm = confusion_matrix(y, preds, labels=[0, 1])

    # Save metrics
    metrics_path = os.path.join(OUT_DIR_TABLES, "classifier_metrics_stable.tsv")
    pd.DataFrame({
        "metric": [
            "roc_auc_cv_mean",
            "roc_auc_cv_std",
            "n_samples",
            "n_splits",
            "n_features_used",
            "min_folds_for_stability",
            "threshold_for_cm",
        ],
        "value": [
            mean_auc,
            std_auc,
            float(len(y)),
            float(N_SPLITS),
            float(Xs.shape[1]),
            float(MIN_FOLDS),
            float(THRESHOLD),
        ]
    }).to_csv(metrics_path, sep="\t", index=False)
    log(f"Wrote: {metrics_path}")

    # Confusion matrix file
    cm_path = os.path.join(OUT_DIR_TABLES, "confusion_matrix_stable.tsv")
    cm_df = pd.DataFrame(
        cm,
        index=["True_Normal(0)", "True_Tumor(1)"],
        columns=["Pred_Normal(0)", "Pred_Tumor(1)"]
    )
    cm_df.to_csv(cm_path, sep="\t")
    log(f"Wrote: {cm_path}")

    # ROC curve
    fpr, tpr, _ = roc_curve(y, all_probs)
    roc_path = os.path.join(OUT_DIR_FIGS, "roc_curve_stable.png")
    plt.figure()
    plt.plot(fpr, tpr)
    plt.plot([0, 1], [0, 1], linestyle="--")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"ROC (Stable genes only) AUC={mean_auc:.3f}±{std_auc:.3f} | n={Xs.shape[1]}")
    plt.savefig(roc_path, dpi=150, bbox_inches="tight")
    plt.close()
    log(f"Saved: {roc_path}")

    # Fit final model on all data to get coefficients for a stable panel
    log("Fitting final model on all data (stable genes) for coefficients...")
    final_model = build_model()
    final_model.fit(Xs, y)

    coef = final_model.named_steps["clf"].coef_[0]
    coef_df = pd.DataFrame({"gene": stable_genes, "coef": coef}).sort_values("coef", ascending=False)

    panel = pd.concat([coef_df.head(PANEL_POS), coef_df.tail(PANEL_NEG)], axis=0)
    panel_path = os.path.join(OUT_DIR_TABLES, "diagnostic_panel_stable_top20.tsv")
    panel.to_csv(panel_path, sep="\t", index=False)
    log(f"Wrote: {panel_path}")

    log("Done.")

if __name__ == "__main__":
    main()
