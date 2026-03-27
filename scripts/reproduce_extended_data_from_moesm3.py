from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


MOESM3_PATH = Path(r"d:\naturedna\evoflux-reproduce\evoflux\41586_2025_9374_MOESM3_ESM.xlsx")
OUT_DIR = Path(r"d:\naturedna\evoflux-reproduce\figures\Extended_Data_Fig")


def _load_moesm3() -> pd.DataFrame:
    meta = pd.read_excel(MOESM3_PATH, sheet_name="supplementary_table_2", header=2)
    qc = pd.read_excel(MOESM3_PATH, sheet_name="supplementary_table_3", header=2)
    inf = pd.read_excel(MOESM3_PATH, sheet_name="supplementary_table_12", header=2)

    # Normalize participant id column name
    meta = meta.rename(columns={"participant_id_anonymous": "PARTICIPANT_ID_ANONYMOUS"})
    qc = qc.rename(columns={"participant_id_anonymous": "PARTICIPANT_ID_ANONYMOUS"})

    df = (
        inf.merge(
            meta[
                [
                    "PARTICIPANT_ID_ANONYMOUS",
                    "sample_id",
                    "cell_type",
                    "cell_type_annot_1",
                    "cell_type_annot_2",
                    "cell_type_annot_3",
                    "disease_subtype",
                    "disease_subtype_epigenetics",
                    "purity_tumor_consensus",
                ]
            ],
            on="PARTICIPANT_ID_ANONYMOUS",
            how="left",
        )
        .merge(
            qc[
                [
                    "PARTICIPANT_ID_ANONYMOUS",
                    "age_sampling",
                    "ssnob_normalization_batch",
                ]
            ],
            on="PARTICIPANT_ID_ANONYMOUS",
            how="left",
        )
    )

    # Derived values commonly used in the paper.
    df["switch_mean"] = df[["mu", "nu", "gamma", "zeta"]].mean(axis=1, skipna=True)
    df["switch_mean_x100"] = 100.0 * df["switch_mean"]
    df["Scancer_per_million"] = df["Scancer"] / 1e6
    df["log10_Scancer"] = np.log10(df["Scancer"].astype(float))
    return df


def _order_diseases(df: pd.DataFrame, col: str = "cell_type_annot_1") -> list[str]:
    # A stable, human-friendly order close to the paper narrative.
    preferred = [
        "B-ALL",
        "T-ALL",
        "B-ALL-remission",
        "T-ALL-remission",
        "B-ALL-relapse",
        "CLL",
        "MBL",
        "RT",
        "MCL",
        "DLBCL-NOS",
        "MM",
        "MGUS",
        "Normal_lymphoid_cell",
        "Whole_blood",
        "PBMCs",
    ]
    present = [x for x in preferred if x in set(df[col].dropna().astype(str))]
    # Append any missing categories in alphabetical order.
    missing = sorted(set(df[col].dropna().astype(str)) - set(present))
    return present + missing


def _save_pdf(fig: plt.Figure, name: str) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / name
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def plot_boxplots_by_disease(df: pd.DataFrame) -> list[Path]:
    sns.set_style("whitegrid")
    order = _order_diseases(df, "cell_type_annot_1")

    # Filter to samples with a disease label.
    d = df[df["cell_type_annot_1"].notna()].copy()

    outputs: list[Path] = []

    specs = [
        ("theta", "Growth rate θ (per year)", "ExtendedData_theta_by_disease.pdf", False),
        ("Scancer_per_million", "Effective population size Nₑ (millions)", "ExtendedData_Ne_by_disease.pdf", True),
        ("cancerAge", "Time since MRCA (T − τ) (years)", "ExtendedData_cancerAge_by_disease.pdf", False),
        ("switch_mean_x100", "Mean switching rate × 100", "ExtendedData_switchrate_by_disease.pdf", True),
    ]

    for y, ylab, fname, logy in specs:
        fig, ax = plt.subplots(figsize=(10.5, 3.6))
        sns.boxplot(
            data=d,
            x="cell_type_annot_1",
            y=y,
            order=order,
            ax=ax,
            fliersize=0,
            linewidth=0.7,
            color="#E6E6E6",
        )
        # Overlay a light jitter for distribution sense (kept subtle).
        sns.stripplot(
            data=d.sample(min(len(d), 4000), random_state=0),
            x="cell_type_annot_1",
            y=y,
            order=order,
            ax=ax,
            size=1.2,
            alpha=0.25,
            color="black",
            jitter=0.25,
        )
        ax.set_xlabel("")
        ax.set_ylabel(ylab)
        ax.tick_params(axis="x", rotation=45, labelsize=7)
        ax.tick_params(axis="y", labelsize=8)
        if logy:
            ax.set_yscale("log")
        outputs.append(_save_pdf(fig, fname))

    return outputs


def plot_scatter_theta_vs_Ne(df: pd.DataFrame) -> Path:
    sns.set_style("whitegrid")
    d = df[df["cell_type_annot_1"].notna()].copy()
    d = d[(d["theta"].notna()) & (d["Scancer"].notna())]

    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    sns.scatterplot(
        data=d,
        x="theta",
        y="Scancer",
        hue="cell_type_annot_1",
        palette="tab10",
        s=12,
        alpha=0.75,
        linewidth=0,
        ax=ax,
        legend=False,
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("θ (per year, log scale)")
    ax.set_ylabel("Nₑ (Scancer, log scale)")
    return _save_pdf(fig, "ExtendedData_scatter_theta_vs_Ne.pdf")


def plot_scatter_tau_vs_switch(df: pd.DataFrame) -> Path:
    sns.set_style("whitegrid")
    d = df[df["cell_type_annot_1"].notna()].copy()
    d = d[(d["tau"].notna()) & (d["switch_mean"].notna())]

    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    sns.scatterplot(
        data=d,
        x="tau",
        y="switch_mean",
        hue="cell_type_annot_1",
        palette="tab10",
        s=12,
        alpha=0.75,
        linewidth=0,
        ax=ax,
        legend=False,
    )
    ax.set_xlabel("τ (age at MRCA, years)")
    ax.set_ylabel("Mean switching rate")
    ax.set_yscale("log")
    return _save_pdf(fig, "ExtendedData_scatter_tau_vs_switch.pdf")


def main() -> None:
    df = _load_moesm3()
    outs = []
    outs.append(plot_scatter_theta_vs_Ne(df))
    outs.append(plot_scatter_tau_vs_switch(df))

    print("WROTE")
    for p in outs:
        print(str(p))


if __name__ == "__main__":
    main()

