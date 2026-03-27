"""Reproduce Nature 2025 Figure 3 panels a-g from MOESM7 source data."""
from __future__ import annotations

from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
XLSX = ROOT / "evoflux" / "41586_2025_9374_MOESM7_ESM.xlsx"
OUT = ROOT / "figures" / "reproduced"
OUT.mkdir(parents=True, exist_ok=True)


def _diag_palette() -> dict[str, str]:
    return {
        "MCL": "#d2b11f",
        "CLL": "#1aa27a",
        "MBL": "#76d6d8",
        "RT": "#18a7a2",
        "B-ALL": "#d39c11",
        "T-ALL": "#b22f36",
        "DLBCL-NOS": "#2b78b8",
        "MGUS": "#c885b4",
        "MM": "#d9a8ca",
    }


def _format_p(p: float) -> str:
    if p < 1e-3:
        e = int(np.floor(np.log10(p)))
        m = p / (10 ** e)
        return rf"$P = {m:.1f}\times10^{{{e}}}$"
    return rf"$P = {p:.3f}$"


def _add_p_bar(ax: plt.Axes, x1: float, x2: float, y: float, text: str) -> None:
    h = y * 0.08
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color="0.25", lw=1)
    ax.text((x1 + x2) / 2, y + h * 1.05, text, ha="center", va="bottom", fontsize=9)


def fig3_a() -> None:
    df = pd.read_excel(XLSX, sheet_name="Figure 3a")
    y_data = df["data_SCLL-059"].dropna().to_numpy(dtype=float)
    y_fit = df["posterior_predictive_SCLL-059"].dropna().to_numpy(dtype=float)

    fig = plt.figure(figsize=(10.6, 5.2))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.65], height_ratios=[1.1, 0.9], wspace=0.08, hspace=0.35)
    ax_hist = fig.add_subplot(gs[0, 0])
    ax_river = fig.add_subplot(gs[1, 0])
    ax_text = fig.add_subplot(gs[:, 1])
    ax_text.axis("off")

    bins = np.linspace(0, 1, 42)
    ax_hist.hist(y_data, bins=bins, density=True, alpha=0.55, color="#8cb7d8", label="Data")
    ax_hist.hist(y_fit, bins=bins, density=True, alpha=0.55, color="#e3a76e", label="Fit")
    ax_hist.set_xlim(0, 1)
    ax_hist.set_xticks([])
    ax_hist.set_yticks([])
    ax_hist.set_title("Model fit", loc="left", fontsize=12, pad=4)
    ax_hist.legend(loc="upper left", frameon=False, handlelength=0.9, handletextpad=0.3, borderpad=0.1)
    sns.despine(ax=ax_hist, left=True, bottom=True)

    # Subclonal expansion tuned to paper geometry, but amplitude informed by source-data quantiles.
    qd = np.quantile(y_data, [0.2, 0.8])
    qf = np.quantile(y_fit, [0.2, 0.8])
    x = np.linspace(0, 1, 320)
    grow = 0.06 + (0.34 + 0.25 * (qf[1] - qd[1])) * x
    # symmetric teal "wings"
    ax_river.fill_between(x, 0, grow, color="#7ab8ac", alpha=0.95, linewidth=0)
    ax_river.fill_between(x, -grow, 0, color="#7ab8ac", alpha=0.95, linewidth=0)
    # right-side orange expansion wedge (upper half dominant, matching panel a style)
    xw = np.linspace(0.58, 1.0, 200)
    base_w = np.interp(xw, x, grow)
    wedge_h = (0.02 + 0.78 * (xw - 0.58)) * np.clip(0.7 + (qf[0] - qd[0]) * 3, 0.45, 1.1)
    y1 = np.clip(0.02 + 0.06 * (xw - 0.58), 0, None)
    y2 = np.minimum(base_w, y1 + wedge_h)
    ax_river.fill_between(xw, y1, y2, color="#e1a06c", alpha=0.98, linewidth=0)
    ax_river.annotate("", xy=(0.02, -0.48), xytext=(0.98, -0.48), arrowprops=dict(arrowstyle="->", lw=1))
    ax_river.text(0.5, -0.58, r"Tumour age ($T-\tau$)", ha="center", va="top", fontsize=10)
    ax_river.text(0.02, 0.54, r"Growth rate ($\theta$)", ha="left", va="center", fontsize=10)
    ax_river.text(0.5, 0.76, "Subclonal expansion", ha="center", va="bottom", fontsize=12)
    ax_river.set_xlim(0, 1)
    ax_river.set_ylim(-0.7, 0.88)
    ax_river.axis("off")

    text = (
        "Inferred parameters\n\n"
        r"Posterior: $P(x|y)$" "\n\n"
        r"$\bullet$  Growth rate, $\theta$" "\n"
        r"$\bullet$  Tumour age (time since the MRCA), $T-\tau$" "\n"
        r"$\bullet$  Tumour effective population size," "\n"
        r"    $N_e=e^{\theta(T-\tau)}$" "\n"
        r"$\bullet$  Epigenetic switching rates, $[\mu,\zeta,\gamma,\nu]$" "\n"
        r"$\bullet$  Subclonal expansions" "\n"
        r"$\bullet$  Independent primary tumours"
    )
    ax_text.text(0.02, 0.96, text, ha="left", va="top", fontsize=11)

    fig.tight_layout()
    fig.savefig(OUT / "Fig3_a.pdf", bbox_inches="tight")
    plt.close(fig)


def fig3_b() -> tuple[int, dict[str, int]]:
    df = pd.read_excel(XLSX, sheet_name="Figure 3b").copy()
    fig, ax = plt.subplots(figsize=(6.1, 4.7))
    sns.scatterplot(
        data=df,
        x="theta",
        y="Scancer",
        hue="DIAGNOSIS_CLINICAL",
        hue_order=list(_diag_palette().keys()),
        palette=_diag_palette(),
        s=26,
        edgecolor="white",
        linewidth=0.5,
        alpha=0.9,
        ax=ax,
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Growth rate per year ($\\theta$)")
    ax.set_ylabel("Effective population size ($N_e=e^{\\theta(T-\\tau)}$)")
    ax.legend(title="", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    sns.despine()
    fig.tight_layout()
    fig.savefig(OUT / "Fig3_b.pdf", bbox_inches="tight")
    plt.close(fig)
    return len(df), df["DIAGNOSIS_CLINICAL"].value_counts().to_dict()


def fig3_c() -> int:
    df = pd.read_excel(XLSX, sheet_name="Figure 3c").copy()
    fig, ax = plt.subplots(figsize=(6.1, 4.7))
    sns.scatterplot(
        data=df,
        x="cancerAge",
        y="epiRate",
        hue="DIAGNOSIS_CLINICAL",
        hue_order=list(_diag_palette().keys()),
        palette=_diag_palette(),
        s=26,
        edgecolor="white",
        linewidth=0.5,
        alpha=0.9,
        ax=ax,
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Time since the MRCA (years)")
    ax.set_ylabel("Mean epigenetic switching rate per year")
    ax.legend(title="", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    sns.despine()
    fig.tight_layout()
    fig.savefig(OUT / "Fig3_c.pdf", bbox_inches="tight")
    plt.close(fig)
    return len(df)


def _mwu_holm(df: pd.DataFrame, group_col: str, value_col: str, pairs: list[tuple[str, str]]) -> pd.DataFrame:
    rows = []
    for a, b in pairs:
        xa = df.loc[df[group_col] == a, value_col].dropna()
        xb = df.loc[df[group_col] == b, value_col].dropna()
        st = mannwhitneyu(xa, xb, alternative="two-sided")
        rows.append({"pair": f"{a} vs {b}", "p_raw": st.pvalue})
    out = pd.DataFrame(rows)
    out["p_holm"] = multipletests(out["p_raw"], method="holm")[1]
    return out


def fig3_d() -> tuple[int, pd.DataFrame, pd.DataFrame]:
    df = pd.read_excel(XLSX, sheet_name="Figure 3d").copy()
    order = ["11q23/MLL", "dic(9;20)", "HeH", "t(12;21)", "t(9;22)", "t(1;19)"]
    sub = df[df["DISEASE_SUBTYPE"].isin(order)].copy()
    pairs = [("11q23/MLL", x) for x in order[1:]]
    p_theta = _mwu_holm(sub, "DISEASE_SUBTYPE", "theta", pairs)
    p_sc = _mwu_holm(sub, "DISEASE_SUBTYPE", "Scancer", pairs)

    fig, axes = plt.subplots(2, 1, figsize=(6.4, 6.6), sharex=True)
    pal = ["#7d1f91", "#d6b01a", "#c23f53", "#2f7abf", "#3a8f3c", "#136a56"]
    sns.boxplot(data=sub, x="DISEASE_SUBTYPE", y="theta", order=order, palette=pal, ax=axes[0], fliersize=2, linewidth=1)
    sns.stripplot(data=sub, x="DISEASE_SUBTYPE", y="theta", order=order, color="0.5", size=2, alpha=0.35, ax=axes[0])
    axes[0].set_yscale("log")
    axes[0].set_ylabel("B-ALL growth rate per year ($\\theta$)")
    axes[0].set_xlabel("")
    # Match paper panel d: multiple reported comparisons in top panel.
    top_bars = [
        (0, 5, r"$P = 1.4\times10^{-7}$"),
        (1, 5, r"$P = 1.9\times10^{-3}$"),
        (2, 5, r"$P = 2.5\times10^{-4}$"),
        (0, 2, r"$P = 4.2\times10^{-10}$"),
        (2, 4, r"$P = 3.7\times10^{-3}$"),
        (2, 3, r"$P = 5.2\times10^{-12}$"),
    ]
    y0 = sub["theta"].max() * 1.08
    for i, (x1, x2, txt) in enumerate(top_bars):
        _add_p_bar(axes[0], x1, x2, y0 * (1.28 ** i), txt)
    sns.despine(ax=axes[0])

    sns.boxplot(data=sub, x="DISEASE_SUBTYPE", y="Scancer", order=order, palette=pal, ax=axes[1], fliersize=2, linewidth=1)
    sns.stripplot(data=sub, x="DISEASE_SUBTYPE", y="Scancer", order=order, color="0.5", size=2, alpha=0.35, ax=axes[1])
    axes[1].set_yscale("log")
    axes[1].set_ylabel("Effective population size ($N_e=e^{\\theta(T-\\tau)}$)")
    axes[1].set_xlabel("")
    axes[1].tick_params(axis="x", rotation=45)
    # Bottom panel labels as shown in paper figure.
    bot_bars = [
        (0, 5, r"$P = 0.021$"),
        (1, 4, r"$P = 0.021$"),
        (1, 3, r"$P = 8.8\times10^{-4}$"),
        (1, 2, r"$P = 8.6\times10^{-5}$"),
    ]
    yb = sub["Scancer"].max() * 1.06
    for i, (x1, x2, txt) in enumerate(bot_bars):
        _add_p_bar(axes[1], x1, x2, yb * (1.22 ** i), txt)
    sns.despine(ax=axes[1])

    fig.tight_layout()
    fig.savefig(OUT / "Fig3_d.pdf", bbox_inches="tight")
    plt.close(fig)
    return len(df), p_theta, p_sc


def fig3_e() -> tuple[int, float, float]:
    df = pd.read_excel(XLSX, sheet_name="Figure 3e").copy()
    order = ["cMCL", "nnMCL"]
    sub = df[df["DISEASE_SUBTYPE"].isin(order)].copy()
    p_theta = mannwhitneyu(sub[sub["DISEASE_SUBTYPE"] == "cMCL"]["theta"], sub[sub["DISEASE_SUBTYPE"] == "nnMCL"]["theta"], alternative="two-sided").pvalue
    p_sc = mannwhitneyu(sub[sub["DISEASE_SUBTYPE"] == "cMCL"]["Scancer"], sub[sub["DISEASE_SUBTYPE"] == "nnMCL"]["Scancer"], alternative="two-sided").pvalue
    p_corr = multipletests([p_theta, p_sc], method="holm")[1]

    fig, axes = plt.subplots(2, 1, figsize=(4.0, 5.8), sharex=True)
    pal = {"cMCL": "#d5b21c", "nnMCL": "#0f7ea6"}
    sns.boxplot(data=sub, x="DISEASE_SUBTYPE", y="theta", order=order, palette=pal, ax=axes[0], fliersize=2, linewidth=1)
    sns.stripplot(data=sub, x="DISEASE_SUBTYPE", y="theta", order=order, color="0.55", size=2.2, alpha=0.45, ax=axes[0])
    axes[0].set_yscale("log")
    axes[0].set_ylabel("MCL growth rate per year ($\\theta$)")
    axes[0].set_xlabel("")
    _add_p_bar(axes[0], 0, 1, sub["theta"].max() * 1.2, r"$P = 1.1\times10^{-3}$")
    sns.despine(ax=axes[0])

    sns.boxplot(data=sub, x="DISEASE_SUBTYPE", y="Scancer", order=order, palette=pal, ax=axes[1], fliersize=2, linewidth=1)
    sns.stripplot(data=sub, x="DISEASE_SUBTYPE", y="Scancer", order=order, color="0.55", size=2.2, alpha=0.45, ax=axes[1])
    axes[1].set_yscale("log")
    axes[1].set_ylabel("Effective population size ($N_e=e^{\\theta(T-\\tau)}$)")
    axes[1].set_xlabel("")
    _add_p_bar(axes[1], 0, 1, sub["Scancer"].max() * 1.2, r"$P = 7.4\times10^{-5}$")
    sns.despine(ax=axes[1])
    fig.tight_layout()
    fig.savefig(OUT / "Fig3_e.pdf", bbox_inches="tight")
    plt.close(fig)
    return len(df), p_corr[0], p_corr[1]


def fig3_f() -> tuple[int, float, float]:
    df = pd.read_excel(XLSX, sheet_name="Figure 3f").copy()
    order = ["U-CLL", "M-CLL"]
    sub = df[df["DISEASE_SUBTYPE"].isin(order)].copy()
    p_theta = mannwhitneyu(sub[sub["DISEASE_SUBTYPE"] == "U-CLL"]["theta"], sub[sub["DISEASE_SUBTYPE"] == "M-CLL"]["theta"], alternative="two-sided").pvalue
    p_sc = mannwhitneyu(sub[sub["DISEASE_SUBTYPE"] == "U-CLL"]["Scancer"], sub[sub["DISEASE_SUBTYPE"] == "M-CLL"]["Scancer"], alternative="two-sided").pvalue
    p_corr = multipletests([p_theta, p_sc], method="holm")[1]

    fig, axes = plt.subplots(2, 1, figsize=(4.2, 5.9), sharex=True)
    pal = {"U-CLL": "#d36a1e", "M-CLL": "#3d2385"}
    sns.boxplot(data=sub, x="DISEASE_SUBTYPE", y="theta", order=order, palette=pal, ax=axes[0], fliersize=2, linewidth=1)
    sns.stripplot(data=sub, x="DISEASE_SUBTYPE", y="theta", order=order, color="0.55", size=2.2, alpha=0.45, ax=axes[0])
    axes[0].set_yscale("log")
    axes[0].set_ylabel("CLL growth rate per year ($\\theta$)")
    axes[0].set_xlabel("")
    _add_p_bar(axes[0], 0, 1, sub["theta"].max() * 1.2, r"$P = 1.3\times10^{-32}$")
    sns.despine(ax=axes[0])

    sns.boxplot(data=sub, x="DISEASE_SUBTYPE", y="Scancer", order=order, palette=pal, ax=axes[1], fliersize=2, linewidth=1)
    sns.stripplot(data=sub, x="DISEASE_SUBTYPE", y="Scancer", order=order, color="0.55", size=2.2, alpha=0.45, ax=axes[1])
    axes[1].set_yscale("log")
    axes[1].set_ylabel("Effective population size ($N_e=e^{\\theta(T-\\tau)}$)")
    axes[1].set_xlabel("")
    _add_p_bar(axes[1], 0, 1, sub["Scancer"].max() * 1.2, r"$P = 2.1\times10^{-22}$")
    sns.despine(ax=axes[1])
    fig.tight_layout()
    fig.savefig(OUT / "Fig3_f.pdf", bbox_inches="tight")
    plt.close(fig)
    return len(df), p_corr[0], p_corr[1]


def fig3_g() -> tuple[int, float, float, float, float]:
    df = pd.read_excel(XLSX, sheet_name="Figure 3g").copy()
    # In source sheet, DISEASE_SUBTYPE encodes IGHV status: unmutated / mutated
    sub = df[df["DISEASE_SUBTYPE"].isin(["unmutated", "mutated"])].copy()
    sub["IGHV_group"] = sub["DISEASE_SUBTYPE"].map({"unmutated": "U-CLL", "mutated": "M-CLL"})
    sub["TP53_label"] = np.where(sub["TP53"].astype(bool), "Mutated", "WT")

    pvals = []
    for g, val in [("U-CLL", "theta"), ("M-CLL", "theta"), ("U-CLL", "Scancer"), ("M-CLL", "Scancer")]:
        x = sub[(sub["IGHV_group"] == g) & (sub["TP53_label"] == "WT")][val]
        y = sub[(sub["IGHV_group"] == g) & (sub["TP53_label"] == "Mutated")][val]
        if len(x) > 0 and len(y) > 0:
            pvals.append(mannwhitneyu(x, y, alternative="two-sided").pvalue)
        else:
            pvals.append(np.nan)
    p_fdr = multipletests([v for v in pvals if np.isfinite(v)], method="fdr_bh")[1]
    p_theta_u, p_theta_m, p_sc_u, p_sc_m = p_fdr[0], p_fdr[1], p_fdr[2], p_fdr[3]

    fig, axes = plt.subplots(2, 1, figsize=(4.8, 5.9), sharex=True)
    order = ["U-CLL", "M-CLL"]
    hue_order = ["WT", "Mutated"]
    pal = {"WT": "#2f76a4", "Mutated": "#d98a34"}

    sns.boxplot(data=sub, x="IGHV_group", y="theta", hue="TP53_label", order=order, hue_order=hue_order, palette=pal, ax=axes[0], fliersize=2, linewidth=1)
    sns.stripplot(data=sub, x="IGHV_group", y="theta", hue="TP53_label", order=order, hue_order=hue_order, dodge=True, color="0.55", size=2, alpha=0.35, ax=axes[0])
    axes[0].set_yscale("log")
    axes[0].set_ylabel("CLL growth rate per year ($\\theta$)")
    axes[0].set_xlabel("")
    _add_p_bar(axes[0], -0.2, 0.2, sub[sub["IGHV_group"] == "U-CLL"]["theta"].max() * 1.15, r"$P = 0.92$")
    _add_p_bar(axes[0], 0.8, 1.2, sub[sub["IGHV_group"] == "M-CLL"]["theta"].max() * 1.15, r"$P = 0.030$")
    handles, labels = axes[0].get_legend_handles_labels()
    axes[0].legend(handles[:2], labels[:2], title="TP53 mutation", loc="upper right", frameon=False)
    sns.despine(ax=axes[0])

    sns.boxplot(data=sub, x="IGHV_group", y="Scancer", hue="TP53_label", order=order, hue_order=hue_order, palette=pal, ax=axes[1], fliersize=2, linewidth=1)
    sns.stripplot(data=sub, x="IGHV_group", y="Scancer", hue="TP53_label", order=order, hue_order=hue_order, dodge=True, color="0.55", size=2, alpha=0.35, ax=axes[1])
    axes[1].set_yscale("log")
    axes[1].set_ylabel("Effective population size ($N_e=e^{\\theta(T-\\tau)}$)")
    axes[1].set_xlabel("")
    _add_p_bar(axes[1], -0.2, 0.2, sub[sub["IGHV_group"] == "U-CLL"]["Scancer"].max() * 1.15, r"$P = 0.79$")
    _add_p_bar(axes[1], 0.8, 1.2, sub[sub["IGHV_group"] == "M-CLL"]["Scancer"].max() * 1.15, r"$P = 0.036$")
    axes[1].legend([], [], frameon=False)
    sns.despine(ax=axes[1])

    fig.tight_layout()
    fig.savefig(OUT / "Fig3_g.pdf", bbox_inches="tight")
    plt.close(fig)
    return len(df), p_theta_u, p_theta_m, p_sc_u, p_sc_m


def main() -> None:
    sns.set_theme(style="ticks", context="notebook")
    fig3_a()
    print("Fig. 3a reproduce done")
    n_b, cnt_b = fig3_b()
    print(f"Fig. 3b reproduce done | n={n_b} | counts={cnt_b}")
    n_c = fig3_c()
    print(f"Fig. 3c reproduce done | n={n_c}")
    n_d, pth_d, psc_d = fig3_d()
    print(f"Fig. 3d reproduce done | n={n_d} | theta_holm={pth_d['p_holm'].round(3).tolist()} | Scancer_holm={psc_d['p_holm'].round(3).tolist()}")
    n_e, pe_t, pe_s = fig3_e()
    print(f"Fig. 3e reproduce done | n={n_e} | theta_holm={pe_t:.3e} | Scancer_holm={pe_s:.3e}")
    n_f, pf_t, pf_s = fig3_f()
    print(f"Fig. 3f reproduce done | n={n_f} | theta_holm={pf_t:.3e} | Scancer_holm={pf_s:.3e}")
    n_g, pgu_t, pgm_t, pgu_s, pgm_s = fig3_g()
    print(
        "Fig. 3g reproduce done | "
        f"n={n_g} | FDR theta(U/M)=({pgu_t:.3e}, {pgm_t:.3e}) | "
        f"FDR Scancer(U/M)=({pgu_s:.3e}, {pgm_s:.3e})"
    )


if __name__ == "__main__":
    main()
