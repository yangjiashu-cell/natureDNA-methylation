"""
Reproduce Nature 2025 Fig. 2 panels a–e from MOESM6 source data only.
Outputs: figures/reproduced/Fig2_a.pdf … Fig2_e.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

# --- paths ---
ROOT = Path(__file__).resolve().parents[1]
XLSX = ROOT / "evoflux" / "41586_2025_9374_MOESM6_ESM.xlsx"
OUT_DIR = ROOT / "figures" / "reproduced"


def fig2_a() -> None:
    """EVOFLUx schematic: pre-MRCA τ, post-MRCA θ, switching rates μ, ν, γ, ζ."""
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")

    # timeline
    ax.annotate(
        "",
        xy=(8.2, 3.0),
        xytext=(1.5, 3.0),
        arrowprops=dict(arrowstyle="->", lw=1.5, color="0.3"),
    )
    ax.text(1.0, 3.15, "Time", fontsize=10, color="0.3")
    ax.axvline(4.5, color="0.5", ls="--", lw=1)
    ax.text(4.55, 5.1, "MRCA", fontsize=10, fontweight="bold")

    # pre-MRCA (single-cell / τ)
    pre = mpatches.FancyBboxPatch(
        (1.6, 2.0),
        2.6,
        2.0,
        boxstyle="round,pad=0.03",
        ec="0.2",
        fc="#e8f4fc",
    )
    ax.add_patch(pre)
    ax.text(2.9, 3.55, "Pre-MRCA", ha="center", fontsize=11, fontweight="bold")
    ax.text(2.9, 3.05, r"effective population size $\tau$", ha="center", fontsize=10)

    # post-MRCA (exponential growth θ)
    post = mpatches.FancyBboxPatch(
        (5.0, 2.0),
        3.0,
        2.0,
        boxstyle="round,pad=0.03",
        ec="0.2",
        fc="#fff3e0",
    )
    ax.add_patch(post)
    ax.text(6.5, 3.55, "Post-MRCA", ha="center", fontsize=11, fontweight="bold")
    ax.text(6.5, 3.05, r"growth rate $\theta$ (exponential)", ha="center", fontsize=10)

    # switching rates box
    sw = mpatches.FancyBboxPatch(
        (2.2, 0.35),
        5.6,
        1.25,
        boxstyle="round,pad=0.03",
        ec="0.2",
        fc="#f3e5f5",
    )
    ax.add_patch(sw)
    ax.text(
        5.0,
        1.35,
        r"Methylation switching: $\mu,\ \nu,\ \gamma,\ \zeta$",
        ha="center",
        fontsize=11,
        fontweight="bold",
    )
    ax.text(
        5.0,
        0.72,
        r"$\mu$: unmeth→meth  ·  $\nu$: meth→unmeth  ·  "
        r"$\gamma$: maintenance (meth)  ·  $\zeta$: maintenance (unmeth)",
        ha="center",
        fontsize=8.5,
        color="0.25",
    )

    ax.set_title("EVOFLUx: fluctuating methylation along a growing lineage", fontsize=12, pad=12)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "Fig2_a.pdf", bbox_inches="tight")
    plt.close(fig)


def fig2_b() -> None:
    """Sheet Figure 2b: recent_mrca vs distant_mrca distributions."""
    df = pd.read_excel(XLSX, sheet_name="Figure 2b")
    long = df.melt(var_name="scenario", value_name="methylation")
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.histplot(
        data=long,
        x="methylation",
        hue="scenario",
        bins=50,
        stat="density",
        common_bins=True,
        element="step",
        ax=ax,
    )
    ax.set_xlabel("Simulated fraction methylated")
    ax.set_ylabel("Density")
    ax.set_title("Distant vs recent MRCA (source data)")
    sns.despine()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "Fig2_b.pdf", bbox_inches="tight")
    plt.close(fig)

    r_mean, d_mean = df["recent_mrca"].mean(), df["distant_mrca"].mean()
    ks = stats.ks_2samp(df["recent_mrca"], df["distant_mrca"])
    print(
        f"[Fig2_b] recent mean={r_mean:.6f}, distant mean={d_mean:.6f}, "
        f"KS D={ks.statistic:.6f}, p={ks.pvalue:.2e}"
    )


def fig2_c() -> None:
    """Sheet Figure 2c: high_growth vs low_growth."""
    df = pd.read_excel(XLSX, sheet_name="Figure 2c")
    long = df.melt(var_name="growth", value_name="methylation")
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.histplot(
        data=long,
        x="methylation",
        hue="growth",
        bins=50,
        stat="density",
        common_bins=True,
        element="step",
        ax=ax,
    )
    ax.set_xlabel("Simulated fraction methylated")
    ax.set_ylabel("Density")
    ax.set_title("High vs low post-MRCA growth (source data)")
    sns.despine()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "Fig2_c.pdf", bbox_inches="tight")
    plt.close(fig)

    ks = stats.ks_2samp(df["high_growth"], df["low_growth"])
    print(
        f"[Fig2_c] high mean={df['high_growth'].mean():.6f}, low mean={df['low_growth'].mean():.6f}, "
        f"KS D={ks.statistic:.6f}, p={ks.pvalue:.2e}"
    )


def _load_figure_2d() -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Row 0 of the sheet is a merged label row; row 1 has real headers (timepoint, clone_id, …, T1–T4).
    Only the first 8 rows carry clone prevalence; all 1000 rows carry T1–T4 methylation.
    """
    raw = pd.read_excel(XLSX, sheet_name="Figure 2d", header=1)
    raw = raw.rename(columns={"Unnamed: 3": "_gap"})
    if "_gap" in raw.columns:
        raw = raw.drop(columns=["_gap"])
    for c in raw.columns:
        raw[c] = pd.to_numeric(raw[c], errors="coerce")
    top = raw.dropna(subset=["timepoint", "clone_id"]).copy()
    bottom = raw.dropna(subset=["T1", "T2", "T3", "T4"], how="any").copy()
    return top, bottom


def fig2_d() -> None:
    """
    Sheet Figure 2d — fixed layout:
    Top: symmetric fish plot (mirror about y=0): teal outer wings + orange inner core from MOESM6 clone freqs.
    Bottom: three panels x=T1 methylation, y=T2/T3/T4; |Δmethylation| yellow→purple; diagonal y=x.
    Output: Fig2_d_fixed.pdf
    """
    top, bottom = _load_figure_2d()
    assert len(bottom) == 1000

    t_order = sorted(top["timepoint"].unique())
    assert len(t_order) == 4
    xpos = np.arange(4, dtype=float)
    labels = ["T1", "T2", "T3", "T4"]

    wide = top.pivot(index="timepoint", columns="clone_id", values="clonal_prev")
    wide = wide.reindex(t_order)
    cids = sorted(wide.columns.astype(float))
    c_low, c_high = min(cids), max(cids)
    frac_a = wide[c_low].to_numpy(dtype=float)
    frac_b = wide[c_high].to_numpy(dtype=float)

    fig = plt.figure(figsize=(10.8, 7.4))
    gs_outer = fig.add_gridspec(2, 1, height_ratios=[1.12, 1.38], hspace=0.36)
    ax_top = fig.add_subplot(gs_outer[0, 0])
    gs_inner = gs_outer[1, 0].subgridspec(1, 4, width_ratios=[1, 1, 1, 0.075], wspace=0.34)
    axes_sc = [fig.add_subplot(gs_inner[0, j]) for j in range(3)]
    cax = fig.add_subplot(gs_inner[0, 3])

    # --- top: vertically symmetric fish plot (mirror about y=0), Nature Fig. 2d style ---
    # frac_a = dominant clone at T1 (teal, outer wings); frac_b = emerging clone (orange, inner core)
    x_anchor = np.array([-0.42, 0.0, 1.0, 2.0, 3.0], dtype=float)
    f_o_anchor = np.array([0.0, frac_b[0], frac_b[1], frac_b[2], frac_b[3]], dtype=float)
    x_dense = np.linspace(-0.42, 3.0, 600)
    f_o = np.interp(x_dense, x_anchor, f_o_anchor)
    win = np.hanning(31)
    win /= win.sum()
    f_o = np.convolve(f_o, win, mode="same")
    f_o = np.clip(f_o, 0.0, 1.0)

    # Envelope: total vertical half-span grows from 0 (left pinch) to 1 by ~T1, then flat
    half_h = np.interp(x_dense, [-0.42, -0.18, 0.0, 3.0], [0.0, 0.35, 1.0, 1.0])
    half_h = np.convolve(half_h, win, mode="same")
    half_h = np.clip(half_h, 0.0, 1.0)

    inner = f_o * half_h
    outer = half_h

    teal = "#6eb8ad"
    orange = "#e8a878"
    bg = "#e6e6e6"

    ax_top.set_facecolor(bg)
    pad_x0, pad_x1, pad_y = -0.55, 3.22, 1.18
    ax_top.add_patch(
        mpatches.FancyBboxPatch(
            (pad_x0, -pad_y),
            pad_x1 - pad_x0,
            2 * pad_y,
            boxstyle="round,pad=0.02,rounding_size=0.06",
            facecolor=bg,
            edgecolor="#bbbbbb",
            linewidth=0.8,
            zorder=0,
        )
    )

    ax_top.fill_between(x_dense, inner, outer, color=teal, linewidth=0, zorder=2)
    ax_top.fill_between(x_dense, -outer, -inner, color=teal, linewidth=0, zorder=2)
    ax_top.fill_between(x_dense, 0.0, inner, color=orange, linewidth=0, zorder=2)
    ax_top.fill_between(x_dense, -inner, 0.0, color=orange, linewidth=0, zorder=2)

    y_line = 1.05
    for i, xc in enumerate(xpos):
        ax_top.plot([xc, xc], [-y_line, y_line], ls=(0, (4, 3.5)), lw=1.05, c="0.15", zorder=3)
    circle_colors = ["#4caf50", "#d4c84a", "#8d6e63", "#e57373"]
    y_circ = 1.14
    for i, xc in enumerate(xpos):
        ax_top.scatter(
            [xc],
            [y_circ],
            s=220,
            c=[circle_colors[i]],
            edgecolors="#444444",
            linewidths=0.9,
            zorder=5,
            clip_on=False,
        )
        ax_top.text(xc, y_circ, labels[i], ha="center", va="center", fontsize=9.5, fontweight="bold", color="0.1", zorder=6)

    ax_top.text(-0.48, 1.02, "Clonal fractions", fontsize=11.5, ha="left", va="center", fontweight="bold")
    ax_top.text(1.29, 1.28, "Timepoints", fontsize=12, ha="center", va="center")
    ax_top.set_xlim(pad_x0 + 0.02, pad_x1 - 0.02)
    ax_top.set_ylim(-1.22, 1.38)
    ax_top.axis("off")

    t1 = bottom["T1"].to_numpy(dtype=float)
    cmap_diff = LinearSegmentedColormap.from_list(
        "yl_purple", ["#fff176", "#ab47bc", "#4a148c"], N=256
    )
    norm = plt.Normalize(vmin=0, vmax=0.25)
    for ax, tj in zip(axes_sc, ["T2", "T3", "T4"]):
        tj_arr = bottom[tj].to_numpy(dtype=float)
        adiff = np.abs(tj_arr - t1)
        ax.scatter(
            t1,
            tj_arr,
            c=adiff,
            cmap=cmap_diff,
            norm=norm,
            s=14,
            alpha=0.88,
            linewidths=0,
            rasterized=True,
        )
        ax.plot([0, 1], [0, 1], color="0.35", lw=0.9, ls="--", zorder=0)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("T1 methylation")
        ax.set_ylabel(f"{tj} methylation")
        sns.despine(ax=ax)

    sm = plt.cm.ScalarMappable(cmap=cmap_diff, norm=norm)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax, extend="max")
    cb.set_label(r"$|\Delta$methylation$|$", fontsize=9)
    cb.ax.tick_params(labelsize=8)

    fig.savefig(OUT_DIR / "Fig2_d_fixed.pdf", bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(
        f"[Fig2_d_fixed] symmetric fish top + scatter bottom; clones={len(cids)}; scatter n={len(bottom)}; "
        f"max|T2-T1|={np.abs(bottom['T2']-bottom['T1']).max():.4f}"
    )


def fig2_e() -> None:
    """
    Sheet Figure 2e — match paper layout:
    Title theme: fraction of fCpGs in intermediate methylation peaks (%).
    Left: overlay histogram T1 (neutral, green) vs T3 (subclonal, brown).
    Right: bar chart of those fractions; χ² test on 2×2 table → P ≈ 7.9×10⁻⁶.
    """
    raw = pd.read_excel(XLSX, sheet_name="Figure 2e", header=1)
    raw = raw.rename(columns={"Unnamed: 2": "_g"})
    if "_g" in raw.columns:
        raw = raw.drop(columns=["_g"])
    sum_col = "Number of simulated fCpGs"
    left = raw[["T1", "T3"]].copy()
    for c in ("T1", "T3"):
        left[c] = pd.to_numeric(left[c], errors="coerce")
    left = left.dropna(subset=["T1", "T3"])

    meta = raw.iloc[:2].copy()
    meta[sum_col] = meta[sum_col].astype(str)
    int_counts = pd.to_numeric(meta["Intermediate"], errors="coerce").to_numpy()
    nint_counts = pd.to_numeric(meta["Non-intermediate"], errors="coerce").to_numpy()
    assert len(left) == 1000
    assert int_counts.tolist() == [73, 135] and nint_counts.tolist() == [927, 865]

    cont = np.array([[int_counts[0], nint_counts[0]], [int_counts[1], nint_counts[1]]], dtype=float)
    chi2, p_chi, _, _ = stats.chi2_contingency(cont)
    n = 1000.0
    frac_inter = int_counts / n * 100.0

    hist_df = pd.concat(
        [
            pd.DataFrame({"methylation": left["T1"], "series": "T1 (neutral)"}),
            pd.DataFrame({"methylation": left["T3"], "series": "T3 (subclonal)"}),
        ],
        ignore_index=True,
    )

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.0))
    fig.suptitle("Fraction of fCpGs in intermediate peaks (%)", fontsize=11.5, y=1.03)

    pal = {"T1 (neutral)": "#7CFC00", "T3 (subclonal)": "#CD853F"}
    sns.histplot(
        data=hist_df,
        x="methylation",
        hue="series",
        hue_order=["T1 (neutral)", "T3 (subclonal)"],
        bins=42,
        stat="density",
        common_bins=True,
        element="step",
        fill=True,
        palette=pal,
        linewidth=1.15,
        alpha=0.55,
        ax=axes[0],
        legend=True,
    )
    leg = axes[0].get_legend()
    if leg is not None:
        for lh in getattr(leg, "legend_handles", []):
            lh.set_alpha(0.75)
    axes[0].set_xlim(0, 1)
    axes[0].set_xlabel("Methylation")
    axes[0].set_ylabel("Density")
    sns.despine(ax=axes[0])

    xb = np.arange(2)
    axes[1].bar(
        xb,
        frac_inter,
        width=0.52,
        color=[pal["T1 (neutral)"], pal["T3 (subclonal)"]],
        edgecolor="#3d3d3d",
        linewidth=0.85,
    )
    axes[1].set_xticks(xb)
    axes[1].set_xticklabels(["T1", "T3"])
    axes[1].set_ylabel("Fraction (%)")
    axes[1].set_ylim(0, max(18, frac_inter.max() * 1.15))
    sns.despine(ax=axes[1])
    axes[1].text(
        0.5,
        0.94,
        r"$P = 7.9 \times 10^{-6}$",
        transform=axes[1].transAxes,
        ha="center",
        va="top",
        fontsize=10,
    )

    fig.tight_layout()
    fig.savefig(OUT_DIR / "Fig2_e.pdf", bbox_inches="tight", dpi=300)
    plt.close(fig)

    p_str = format(p_chi, ".6e").replace("e", "e").replace("\u2212", "-")
    print(
        f"[Fig2_e] n={len(left)}; intermediate % T1={frac_inter[0]:.2f}, T3={frac_inter[1]:.2f}; "
        f"chi2={chi2:.4f}, p={p_str} (reported in figure as 7.9e-6)"
    )


def main(de_only: bool = False) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="white", context="notebook")
    if not de_only:
        fig2_a()
        print("Fig. 2a reproduce done")
        fig2_b()
        print("Fig. 2b reproduce done")
        fig2_c()
        print("Fig. 2c reproduce done")
    fig2_d()
    print("Fig. 2d fixed -> Fig2_d_fixed.pdf")
    fig2_e()
    print("Fig. 2e reproduce done")


if __name__ == "__main__":
    import sys

    main(de_only="--de-only" in sys.argv)
