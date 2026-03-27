"""
Figure 4c: per–time-point methylation fraction distributions as bar histograms (not KDE).
Data: MOESM8.xlsx sheet "Figure 4c" — columns SW-BCP-ALL-375 / 725 / 376 (T1 / T2 / T3).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[1]
XLSX = ROOT / "evoflux" / "41586_2025_9374_MOESM8_ESM.xlsx"
OUT_DIR = ROOT / "figures" / "reproduced"

# Match Nature-style panel labels (years from figure caption).
PANELS: list[tuple[str, str, str]] = [
    ("SW-BCP-ALL-375", "T1, diagnosis, 3.5 years", "B-ALL"),
    ("SW-BCP-ALL-725", "T2, remission, 3.8 years", "Remission"),
    ("SW-BCP-ALL-376", "T3, relapse, 5.8 years", "B-ALL"),
]

COLOR_BALL = "#DAA520"  # golden yellow
COLOR_REMISSION = "#D9C9A8"  # light tan / beige


def _load_fractions() -> dict[str, np.ndarray]:
    df = pd.read_excel(XLSX, sheet_name="Figure 4c", header=0)
    out: dict[str, np.ndarray] = {}
    for col, _, _ in PANELS:
        if col not in df.columns:
            raise KeyError(f"Missing column {col!r} in Figure 4c. Got: {df.columns.tolist()}")
        s = pd.to_numeric(df[col], errors="coerce").dropna().to_numpy(dtype=float)
        s = s[np.isfinite(s)]
        s = np.clip(s, 0.0, 1.0)
        out[col] = s
    return out


def build() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "axes.linewidth": 0.8,
        }
    )

    data = _load_fractions()
    n_bins = 40
    bin_edges = np.linspace(0.0, 1.0, n_bins + 1)
    bin_width = bin_edges[1] - bin_edges[0]
    centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0

    fig, axes = plt.subplots(1, 3, figsize=(9.2, 2.9), sharey=True, sharex=True)

    for ax, (col, title, group) in zip(axes, PANELS):
        vals = data[col]
        dens, _ = np.histogram(vals, bins=bin_edges, density=True)
        color = COLOR_BALL if group == "B-ALL" else COLOR_REMISSION
        ax.bar(
            centers,
            dens,
            width=bin_width * 0.92,
            align="center",
            color=color,
            edgecolor="0.25",
            linewidth=0.35,
        )
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 5.0)
        ax.set_title(title, fontsize=9)
        ax.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
        ax.tick_params(axis="both", labelsize=8)

    axes[0].set_ylabel("Probability density", fontsize=9)
    axes[1].set_xlabel("Fraction methylated", fontsize=9)

    leg = fig.legend(
        handles=[
            Patch(facecolor=COLOR_BALL, edgecolor="0.25", linewidth=0.35, label="B-ALL"),
            Patch(facecolor=COLOR_REMISSION, edgecolor="0.25", linewidth=0.35, label="Remission"),
        ],
        loc="upper right",
        bbox_to_anchor=(0.99, 0.98),
        frameon=True,
        fontsize=8,
    )
    leg.get_frame().set_linewidth(0.5)

    fig.tight_layout(rect=(0, 0, 0.88, 1))
    for ext in ("png", "pdf"):
        fig.savefig(
            OUT_DIR / f"Fig4c.{ext}",
            dpi=600 if ext == "png" else None,
            bbox_inches="tight",
            facecolor="white",
        )
    plt.close(fig)
    print(f"[完成] {OUT_DIR / 'Fig4c.png'}")
    print(f"[完成] {OUT_DIR / 'Fig4c.pdf'}")


if __name__ == "__main__":
    build()
