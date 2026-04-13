#!/usr/bin/env python3
"""
(Upgrade #3) Create a single pipeline schematic figure for the README.

Output:
  - results/figures/pipeline_overview.png
"""

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch


OUT = "results/figures/pipeline_overview.png"


def box(ax, x, y, w, h, title, lines):
    rect = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.02,rounding_size=0.02",
        linewidth=1.2,
        facecolor="white"
    )
    ax.add_patch(rect)
    ax.text(x + w*0.03, y + h*0.78, title, fontsize=11, fontweight="bold", va="top")
    ax.text(x + w*0.03, y + h*0.68, "\n".join(lines), fontsize=9, va="top")


def arrow(ax, x1, y1, x2, y2):
    ax.annotate(
        "", xy=(x2, y2), xytext=(x1, y1),
        arrowprops=dict(arrowstyle="->", lw=1.2)
    )


def main():
    os.makedirs("results/figures", exist_ok=True)

    fig = plt.figure(figsize=(10, 6))
    ax = plt.gca()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # Boxes
    box(ax, 0.05, 0.72, 0.40, 0.22, "Data + Labels",
        ["GSE62944 (TCGA processed counts)",
         "LUAD Tumor vs Normal (n=555)",
         "Build coldata + subset counts"])

    box(ax, 0.55, 0.72, 0.40, 0.22, "Differential Expression (DESeq2)",
        ["Filter low counts",
         "Tumor vs Normal LFC + padj",
         "Plots: MA / PCA / Volcano / Heatmap"])

    box(ax, 0.05, 0.40, 0.40, 0.22, "Diagnostic Panel",
        ["Rank-stable 20-gene logistic regression",
         "Leakage-aware CV AUC",
         "ROC curve + confusion matrix"])

    box(ax, 0.55, 0.40, 0.40, 0.22, "Therapeutic Shortlist",
        ["TargetScore = |LFC| + -log10(padj) + expr",
         "Top 50 hypotheses",
         "Flag known actionable genes"])

    box(ax, 0.30, 0.08, 0.40, 0.22, "README Summary Outputs",
        ["Summary tables (top10 up/down, counts)",
         "Key boxplots for panel genes",
         "Figures embedded for quick review"])

    # Arrows
    arrow(ax, 0.45, 0.83, 0.55, 0.83)   # Data -> DE
    arrow(ax, 0.25, 0.72, 0.25, 0.62)   # Data -> Panel
    arrow(ax, 0.75, 0.72, 0.75, 0.62)   # DE -> Targets
    arrow(ax, 0.25, 0.40, 0.45, 0.30)   # Panel -> README
    arrow(ax, 0.75, 0.40, 0.55, 0.30)   # Targets -> README
    arrow(ax, 0.75, 0.62, 0.55, 0.62)   # (visual balance)

    plt.tight_layout()
    plt.savefig(OUT, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved: {OUT}")


if __name__ == "__main__":
    main()
