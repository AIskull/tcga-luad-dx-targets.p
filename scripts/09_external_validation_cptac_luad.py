#!/usr/bin/env python3
"""
External validation on CPTAC LUAD (tumor vs NAT/normal).

Trains on TCGA LUAD (your processed counts) using a fixed gene panel
then tests on CPTAC LUAD transcriptomics.

Key design choices:
- NO clinical calls to CPTAC (avoids generator/len() bug in your cptac install).
- CPTAC labels are built by loading transcriptomics split by tissue_type:
    NAT -> 0, tumor -> 1
- Handles CPTAC duplicate gene columns (multi-transcript -> same gene symbol).
- Reports BOTH:
    (a) metrics + confusion at fixed threshold 0.5
    (b) optimized threshold via Youden’s J + confusion at that threshold

Inputs (TCGA):
  - data/processed/luad_counts.tsv
  - metadata/luad_coldata.tsv
  - results/tables/diagnostic_panel_rankstable_top20.tsv

Outputs:
  - results/tables/external_metrics_cptac_luad.tsv
  - results/tables/external_confusion_cptac_luad.tsv
  - results/tables/external_confusion_cptac_luad_optimal.tsv
  - results/figures/roc_curve_external_cptac_luad.png
"""

from __future__ import annotations

import os
import numpy as np
import pandas as pd
from pathlib import Path

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve, confusion_matrix

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# -------------------------
# Paths / config
# -------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = Path(os.environ.get("TCGA_LUAD_ROOT", SCRIPT_DIR.parent)).resolve()

TCGA_COUNTS = str(PROJECT_ROOT / "data/processed/luad_counts.tsv")
TCGA_COLDATA = str(PROJECT_ROOT / "metadata/luad_coldata.tsv")
PANEL_PATH = str(PROJECT_ROOT / "results/tables/diagnostic_panel_rankstable_top20.tsv")

OUT_METRICS = str(PROJECT_ROOT / "results/tables/external_metrics_cptac_luad.tsv")
OUT_CM = str(PROJECT_ROOT / "results/tables/external_confusion_cptac_luad.tsv")
OUT_CM_OPT = str(PROJECT_ROOT / "results/tables/external_confusion_cptac_luad_optimal.tsv")
OUT_ROC = str(PROJECT_ROOT / "results/figures/roc_curve_external_cptac_luad.png")

DEFAULT_THRESHOLD = 0.5


# -------------------------
# Utilities
# -------------------------
def must_exist(path: str) -> None:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing required file: {path}")


def drop_duplicate_columns(df: pd.DataFrame, label: str = "") -> pd.DataFrame:
    """Keep first occurrence of any duplicated column names."""
    if df.columns.duplicated().any():
        dupes = df.columns[df.columns.duplicated()].unique().tolist()
        prefix = f"{label}: " if label else ""
        print(f"WARNING: {prefix}dropping duplicated columns (keeping first): {dupes}")
        df = df.loc[:, ~df.columns.duplicated(keep="first")].copy()
    return df


def compute_basic_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    """Return common classification metrics at a fixed threshold."""
    # Confusion components
    # y_true/y_pred are 0/1
    tn = float(((y_true == 0) & (y_pred == 0)).sum())
    fp = float(((y_true == 0) & (y_pred == 1)).sum())
    fn = float(((y_true == 1) & (y_pred == 0)).sum())
    tp = float(((y_true == 1) & (y_pred == 1)).sum())

    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else float("nan")  # TPR / recall
    specificity = tn / (tn + fp) if (tn + fp) > 0 else float("nan")  # TNR
    precision = tp / (tp + fp) if (tp + fp) > 0 else float("nan")    # PPV
    accuracy = (tp + tn) / (tp + tn + fp + fn) if (tp + tn + fp + fn) > 0 else float("nan")
    balanced_accuracy = (sensitivity + specificity) / 2 if np.isfinite(sensitivity) and np.isfinite(specificity) else float("nan")

    return {
        "tn": tn, "fp": fp, "fn": fn, "tp": tp,
        "sensitivity": float(sensitivity),
        "specificity": float(specificity),
        "precision": float(precision),
        "accuracy": float(accuracy),
        "balanced_accuracy": float(balanced_accuracy),
    }


def pick_threshold_youdenJ(y_true: np.ndarray, probs: np.ndarray) -> tuple[float, float, float]:
    """
    Choose threshold maximizing Youden’s J = TPR - FPR.
    Returns (threshold, sensitivity_at_thr, specificity_at_thr)
    """
    fpr, tpr, thresholds = roc_curve(y_true, probs)
    j = tpr - fpr
    best_idx = int(np.argmax(j))
    thr = float(thresholds[best_idx])
    sens = float(tpr[best_idx])
    spec = float(1.0 - fpr[best_idx])
    return thr, sens, spec


# -------------------------
# TCGA loader
# -------------------------
def load_tcga_train(panel_genes: list[str]) -> tuple[pd.DataFrame, pd.Series]:
    print(f"Reading TCGA counts: {TCGA_COUNTS}")
    counts = pd.read_csv(TCGA_COUNTS, sep="\t")
    counts = counts.set_index(counts.columns[0])  # rows=genes, cols=samples

    print(f"Reading TCGA coldata: {TCGA_COLDATA}")
    col = pd.read_csv(TCGA_COLDATA, sep="\t").set_index("sample_id")
    y = (col["condition"].astype(str) == "Tumor").astype(int)

    common = col.index.intersection(counts.columns)
    if len(common) == 0:
        raise ValueError("No overlapping samples between TCGA counts and coldata.")

    X = counts[common].T  # samples x genes
    y = y.loc[common]

    X = np.log2(X + 1)

    genes_present = [g for g in panel_genes if g in X.columns]
    if len(genes_present) < 5:
        raise ValueError(f"Too few panel genes present in TCGA matrix: {len(genes_present)}")

    X = X[genes_present]
    X = drop_duplicate_columns(X, label="TCGA")
    print(f"TCGA train matrix: n={X.shape[0]} samples, p={X.shape[1]} panel genes")
    return X, y


# -------------------------
# CPTAC helpers (NO clinical calls)
# -------------------------
def _flatten_gene_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    CPTAC omics sometimes return MultiIndex columns (e.g., ('TP53','ENSG...') or with a 'Name' level).
    We want columns to be gene symbols if possible.
    """
    if isinstance(df.columns, pd.MultiIndex):
        if "Name" in df.columns.names:
            lvl = df.columns.names.index("Name")
            df = df.copy()
            df.columns = df.columns.get_level_values(lvl).astype(str)
        else:
            df = df.copy()
            df.columns = df.columns.get_level_values(0).astype(str)
    else:
        df = df.copy()
        df.columns = df.columns.astype(str)
    return df


def _try_get_tx(ds, tissue_type: str | None = None) -> pd.DataFrame:
    """
    Robust transcriptomics loader.

    IMPORTANT: We do NOT call ds.get_clinical() or ds.get_dataframe("clinical")
    because your cptac install throws: TypeError: generator has no len()

    We only call get_transcriptomics() with optional tissue_type and source fallbacks.
    """
    errors: list[str] = []
    kwargs = {}
    if tissue_type is not None:
        kwargs["tissue_type"] = tissue_type

    # Attempt 1: plain call
    try:
        tx = ds.get_transcriptomics(**kwargs)
        return _flatten_gene_columns(tx)
    except Exception as e:
        errors.append(f"default get_transcriptomics({kwargs}) failed: {repr(e)}")

    # Attempt 2: try common sources explicitly
    for src in ["washu", "broad", "harmonized"]:
        try:
            tx = ds.get_transcriptomics(source=src, **kwargs)
            print(f"Loaded CPTAC transcriptomics with source='{src}', tissue_type={tissue_type}")
            return _flatten_gene_columns(tx)
        except Exception as e:
            errors.append(f"get_transcriptomics(source='{src}', {kwargs}) failed: {repr(e)}")

    raise RuntimeError("Could not load CPTAC transcriptomics.\n" + "\n".join(errors))


def load_cptac_luad(panel_genes: list[str]) -> tuple[pd.DataFrame, pd.Series] | tuple[None, None]:
    """
    Load CPTAC LUAD transcriptomics WITHOUT ANY CLINICAL CALLS.

    Build labels by loading transcriptomics split by tissue type:
      - NAT/normal: tissue_type="nat"  -> y=0
      - Tumor:      tissue_type="tumor"-> y=1
    """
    try:
        import cptac
    except Exception:
        print(
            "\nERROR: 'cptac' is not installed in this environment.\n"
            "Fix:\n"
            "  conda activate luad_tcga\n"
            "  pip install git+https://github.com/PayneLab/cptac.git\n"
            "Then rerun:\n"
            "  python scripts/09_external_validation_cptac_luad.py\n"
        )
        return None, None

    ds = cptac.Luad()

    nat_candidates = ["nat", "normal", "adjacent_normal"]
    tumor_candidates = ["tumor", "primary", "primary_tumor"]

    nat_tx = None
    nat_errs: list[str] = []
    for tt in nat_candidates:
        try:
            nat_tx = _try_get_tx(ds, tissue_type=tt)
            print(f"Loaded NAT transcriptomics via tissue_type='{tt}'")
            break
        except Exception as e:
            nat_errs.append(f"{tt}: {repr(e)}")

    tumor_tx = None
    tumor_errs: list[str] = []
    for tt in tumor_candidates:
        try:
            tumor_tx = _try_get_tx(ds, tissue_type=tt)
            print(f"Loaded tumor transcriptomics via tissue_type='{tt}'")
            break
        except Exception as e:
            tumor_errs.append(f"{tt}: {repr(e)}")

    if nat_tx is None or tumor_tx is None:
        msg = ["Could not load CPTAC transcriptomics split by tissue_type."]
        if nat_tx is None:
            msg.append("NAT attempts:\n" + "\n".join(nat_errs))
        if tumor_tx is None:
            msg.append("Tumor attempts:\n" + "\n".join(tumor_errs))
        raise RuntimeError("\n\n".join(msg))

    print(f"CPTAC NAT tx shape (raw): {nat_tx.shape}")
    print(f"CPTAC tumor tx shape (raw): {tumor_tx.shape}")

    # numeric coercion -> log2(x+1)
    nat_tx = nat_tx.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    tumor_tx = tumor_tx.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    nat_tx = np.log2(nat_tx + 1)
    tumor_tx = np.log2(tumor_tx + 1)

    # CPTAC may have duplicate gene columns after flattening
    nat_tx = drop_duplicate_columns(nat_tx, label="CPTAC_NAT")
    tumor_tx = drop_duplicate_columns(tumor_tx, label="CPTAC_TUMOR")

    panel_genes = [str(g) for g in panel_genes]
    genes_present = [g for g in panel_genes if (g in nat_tx.columns and g in tumor_tx.columns)]
    if len(genes_present) < 5:
        raise ValueError(
            f"Too few panel genes present in CPTAC matrix: {len(genes_present)}\n"
            f"Example CPTAC columns (first 30): {list(nat_tx.columns)[:30]}"
        )

    nat_tx = nat_tx[genes_present].copy()
    tumor_tx = tumor_tx[genes_present].copy()

    # Make indices unique so NAT and Tumor rows never collide
    nat_tx.index = nat_tx.index.astype(str) + "_NAT"
    tumor_tx.index = tumor_tx.index.astype(str) + "_TUMOR"

    X = pd.concat([nat_tx, tumor_tx], axis=0)
    y = pd.Series(
        data=[0] * nat_tx.shape[0] + [1] * tumor_tx.shape[0],
        index=X.index,
        name="label"
    )

    print("CPTAC label counts:", y.value_counts().to_dict())
    print(f"CPTAC test matrix: n={X.shape[0]} samples, p={X.shape[1]} panel genes")
    return X, y


# -------------------------
# Main
# -------------------------
def main() -> None:
    print(f"Using project root: {PROJECT_ROOT}")
    os.makedirs(PROJECT_ROOT / "results/tables", exist_ok=True)
    os.makedirs(PROJECT_ROOT / "results/figures", exist_ok=True)

    must_exist(TCGA_COUNTS)
    must_exist(TCGA_COLDATA)
    must_exist(PANEL_PATH)

    panel = pd.read_csv(PANEL_PATH, sep="\t")
    if "gene" not in panel.columns:
        raise ValueError("Panel file must contain a 'gene' column.")
    panel_genes = panel["gene"].astype(str).tolist()
    print(f"Loaded panel genes: {len(panel_genes)}")

    # Load
    X_train, y_train = load_tcga_train(panel_genes)
    X_test, y_test = load_cptac_luad(panel_genes)
    if X_test is None:
        return

    # Align genes (intersection)
    common_genes = [g for g in X_train.columns if g in X_test.columns]
    if len(common_genes) < 5:
        raise ValueError("After alignment, too few genes overlap between TCGA and CPTAC.")

    X_train = X_train[common_genes].copy()
    X_test = X_test[common_genes].copy()

    # Defensive: drop duplicates (should be rare for TCGA, common for CPTAC)
    X_train = drop_duplicate_columns(X_train, label="TCGA_aligned")
    X_test = drop_duplicate_columns(X_test, label="CPTAC_aligned")

    # Force identical column order
    X_test = X_test.reindex(columns=X_train.columns, fill_value=0.0)
    print(f"Using aligned genes: {X_train.shape[1]}")

    # Model
    model = Pipeline([
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(max_iter=5000, solver="lbfgs"))
    ])

    # Use numpy arrays to avoid sklearn feature-name checks
    model.fit(X_train.to_numpy(), y_train.to_numpy())
    probs = model.predict_proba(X_test.to_numpy())[:, 1]
    auc = float(roc_auc_score(y_test.to_numpy(), probs))

    # Thresholding: fixed 0.5
    y_true = y_test.to_numpy().astype(int)
    pred_default = (probs >= DEFAULT_THRESHOLD).astype(int)
    cm_default = confusion_matrix(y_true, pred_default, labels=[0, 1])

    # Thresholding: Youden's J (balanced sensitivity/specificity)
    opt_thr, opt_sens, opt_spec = pick_threshold_youdenJ(y_true, probs)
    pred_opt = (probs >= opt_thr).astype(int)
    cm_opt = confusion_matrix(y_true, pred_opt, labels=[0, 1])

    # Save confusion matrices
    pd.DataFrame(
        cm_default,
        index=["True_Normal(0)", "True_Tumor(1)"],
        columns=["Pred_Normal(0)", "Pred_Tumor(1)"]
    ).to_csv(OUT_CM, sep="\t")

    pd.DataFrame(
        cm_opt,
        index=["True_Normal(0)", "True_Tumor(1)"],
        columns=["Pred_Normal(0)", "Pred_Tumor(1)"]
    ).to_csv(OUT_CM_OPT, sep="\t")

    # Metrics (default + optimal)
    m_def = compute_basic_metrics(y_true, pred_default)
    m_opt = compute_basic_metrics(y_true, pred_opt)

    metrics_rows = [
        ("external_dataset", "CPTAC_LUAD"),
        ("roc_auc", auc),
        ("n_test_samples", float(len(y_true))),
        ("n_panel_genes_used", float(X_train.shape[1])),

        ("default_threshold", float(DEFAULT_THRESHOLD)),
        ("default_sensitivity", m_def["sensitivity"]),
        ("default_specificity", m_def["specificity"]),
        ("default_precision", m_def["precision"]),
        ("default_accuracy", m_def["accuracy"]),
        ("default_balanced_accuracy", m_def["balanced_accuracy"]),

        ("optimal_threshold_youdenJ", float(opt_thr)),
        ("optimal_sensitivity", float(opt_sens)),
        ("optimal_specificity", float(opt_spec)),
        ("optimal_precision", m_opt["precision"]),
        ("optimal_accuracy", m_opt["accuracy"]),
        ("optimal_balanced_accuracy", m_opt["balanced_accuracy"]),
    ]

    pd.DataFrame(metrics_rows, columns=["metric", "value"]).to_csv(OUT_METRICS, sep="\t", index=False)

    # ROC plot (include optimal point)
    fpr, tpr, thresholds = roc_curve(y_true, probs)
    plt.figure()
    plt.plot(fpr, tpr, label=f"ROC (AUC={auc:.3f})")
    plt.plot([0, 1], [0, 1], linestyle="--", label="Chance")

    # mark optimal point
    # find closest threshold index (for plotting marker)
    best_idx = int(np.argmin(np.abs(thresholds - opt_thr)))
    plt.scatter([fpr[best_idx]], [tpr[best_idx]], marker="o", label=f"Youden J thr={opt_thr:.3f}")

    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("External ROC: CPTAC LUAD")
    plt.legend(loc="lower right")
    plt.savefig(OUT_ROC, dpi=150, bbox_inches="tight")
    plt.close()

    # Console summary
    print(f"\nAUC: {auc:.3f}")
    print(f"Default threshold {DEFAULT_THRESHOLD:.2f} -> sensitivity={m_def['sensitivity']:.3f}, specificity={m_def['specificity']:.3f}")
    print(f"Optimal threshold (Youden J) {opt_thr:.3f} -> sensitivity={m_opt['sensitivity']:.3f}, specificity={m_opt['specificity']:.3f}")

    print(f"\nWrote: {OUT_METRICS}")
    print(f"Wrote: {OUT_CM}")
    print(f"Wrote: {OUT_CM_OPT}")
    print(f"Saved: {OUT_ROC}")
    print("Done.")


if __name__ == "__main__":
    main()
