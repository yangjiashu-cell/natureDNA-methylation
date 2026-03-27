from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Patch, Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


MOESM5 = Path(r"D:\naturedna\evoflux-reproduce\evoflux\41586_2025_9374_MOESM5_ESM.xlsx")
OUT_PDF = Path(r"D:\naturedna\evoflux-reproduce\figures\reproduced\Fig1_c.pdf")


PAPER_COLS = {
    "B cell": "#BDBDBD",
    "T cell": "#FFFFFF",
    "PBMCs": "#CFE8F3",
    "Whole blood": "#B8B5A6",
    "T-ALL": "#CD3333",
    "T-ALL remission": "#C7A5A5",
    "T-ALL relapse": "#CD3333",
    "B-ALL": "#E69F00",
    "B-ALL remission": "#D8CAAB",
    "B-ALL relapse": "#E69F00",
    "MCL": "#F0E442",
    "DLBCL-NOS": "#56B4E9",
    "MBL": "#5A9A89",
    "CLL": "#009E73",
    "RT": "#02C590",
    "MGUS": "#BF99AE",
    "MM": "#CC79A7",
}


def _auto_start_row(df: pd.DataFrame) -> int:
    for i in range(min(200, len(df))):
        v = df.iloc[i, 0]
        if isinstance(v, str) and v.startswith("cg"):
            return i
    raise RuntimeError("Could not find CpG start row (no 'cg*' in first column).")


def _read_figure1c() -> tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_excel(MOESM5, sheet_name="Figure 1c", header=None)
    start = _auto_start_row(df)

    sample_ids = df.iloc[0, 1:].astype(str).values
    sample_group = df.iloc[1, 1:].astype(str).values
    purity = pd.to_numeric(df.iloc[2, 1:], errors="coerce").values.astype(float)
    fcpg_discovery = df.iloc[3, 1:].astype(str).str.lower().isin(["true", "1", "t", "yes"]).values
    is_cancer = df.iloc[4, 1:].astype(str).str.lower().isin(["true", "1", "t", "yes"]).values
    platform = df.iloc[5, 1:].astype(str).values

    cpg_ids = df.iloc[start:, 0].astype(str).values
    beta = df.iloc[start:, 1:].astype(float).values

    beta_df = pd.DataFrame(beta, index=cpg_ids, columns=sample_ids)
    meta = pd.DataFrame(
        {
            "Sample group": sample_group,
            "Tumour fraction": purity,
            "fCpG discovery": np.where(fcpg_discovery, "True", "False"),
            "Sample type": np.where(is_cancer, "Cancer", "Normal"),
            "Platform": np.where(pd.Series(platform).str.contains("EPIC", case=False, na=False), "EPIC", "450k"),
        },
        index=sample_ids,
    )
    return beta_df, meta


def _build_col_colors(meta: pd.DataFrame) -> pd.DataFrame:
    # Sample group colors (fallback grey if unseen)
    sg = meta["Sample group"].astype(str)
    sg_colors = [PAPER_COLS.get(x, "#CCCCCC") for x in sg]

    # Tumour fraction continuous (Blues)
    tf = meta["Tumour fraction"].astype(float).clip(0, 1)
    tf_cmap = plt.cm.Blues
    tf_norm = Normalize(vmin=0, vmax=1)
    tf_colors = [tf_cmap(tf_norm(v)) for v in tf]

    # fCpG discovery True/False (teal/yellow)
    disc = meta["fCpG discovery"].astype(str)
    disc_map = {"True": "#00A6A6", "False": "#E6D34A"}
    disc_colors = [disc_map.get(x, "#666666") for x in disc]

    # Sample type Cancer/Normal (red/blue)
    st = meta["Sample type"].astype(str)
    st_map = {"Cancer": "#D73027", "Normal": "#4575B4"}
    st_colors = [st_map.get(x, "#666666") for x in st]

    # Platform 450k/EPIC (purple/green)
    pf = meta["Platform"].astype(str)
    pf_map = {"450k": "#6A3D9A", "EPIC": "#33A02C"}
    pf_colors = [pf_map.get(x, "#666666") for x in pf]

    col_colors = pd.DataFrame(
        {
            "Sample group": sg_colors,
            "Tumour fraction": tf_colors,
            "fCpG discovery": disc_colors,
            "Sample type": st_colors,
            "Platform": pf_colors,
        },
        index=meta.index,
    )
    return col_colors


def _add_legends(fig: plt.Figure) -> None:
    # Sample group legend (left-top)
    # Build handles only for categories present in the sheet
    handles = []
    for k in [
        "B cell",
        "T cell",
        "PBMCs",
        "Whole blood",
        "T-ALL",
        "T-ALL remission",
        "T-ALL relapse",
        "B-ALL",
        "B-ALL remission",
        "B-ALL relapse",
        "DLBCL-NOS",
        "MBL",
        "CLL",
        "RT",
        "MGUS",
        "MM",
        "MCL",
    ]:
        if k in PAPER_COLS:
            handles.append(Patch(facecolor=PAPER_COLS[k], edgecolor="none", label=k))

    fig.legend(
        handles=handles,
        title="Sample group",
        loc="upper left",
        bbox_to_anchor=(0.02, 0.99),
        ncol=6,
        frameon=False,
        fontsize=8,
        title_fontsize=9,
        handlelength=1.2,
        handleheight=0.8,
        columnspacing=0.9,
    )

    # Right legends
    fig.legend(
        handles=[
            Patch(facecolor="#00A6A6", edgecolor="none", label="True"),
            Patch(facecolor="#E6D34A", edgecolor="none", label="False"),
        ],
        title="fCpG discovery",
        loc="upper right",
        bbox_to_anchor=(0.995, 0.86),
        frameon=False,
        fontsize=8,
        title_fontsize=9,
    )
    fig.legend(
        handles=[
            Patch(facecolor="#D73027", edgecolor="none", label="Cancer"),
            Patch(facecolor="#4575B4", edgecolor="none", label="Normal"),
        ],
        title="Sample type",
        loc="upper right",
        bbox_to_anchor=(0.995, 0.70),
        frameon=False,
        fontsize=8,
        title_fontsize=9,
    )
    fig.legend(
        handles=[
            Patch(facecolor="#6A3D9A", edgecolor="none", label="450k"),
            Patch(facecolor="#33A02C", edgecolor="none", label="EPIC"),
        ],
        title="Platform",
        loc="upper right",
        bbox_to_anchor=(0.995, 0.56),
        frameon=False,
        fontsize=8,
        title_fontsize=9,
    )


def _add_insets(g: sns.matrix.ClusterGrid, meta_reordered: pd.DataFrame) -> None:
    ax = g.ax_heatmap
    n_rows, n_cols = g.data2d.shape

    # Identify cancer vs normal columns after clustering
    is_cancer = meta_reordered["Sample type"].values == "Cancer"
    cancer_cols = np.where(is_cancer)[0]
    normal_cols = np.where(~is_cancer)[0]

    # Pick a representative window in each group
    def pick_window(cols: np.ndarray, width: int) -> tuple[int, int]:
        if len(cols) == 0:
            return 0, min(width, n_cols)
        # choose a window centered in the group
        start = int(max(cols[0], cols[len(cols) // 2] - width // 2))
        end = min(start + width, n_cols)
        start = max(0, end - width)
        return start, end

    # Use moderate window sizes for readability
    cx0, cx1 = pick_window(cancer_cols, width=min(180, n_cols))
    nx0, nx1 = pick_window(normal_cols, width=min(120, n_cols))

    # Row window: take a mid-band for "speckled" visualization
    ry0 = int(n_rows * 0.45)
    ry1 = int(min(n_rows, ry0 + 160))

    # Draw rectangles on the main heatmap
    rect_kwargs = dict(fill=False, edgecolor="black", linewidth=0.6)
    ax.add_patch(Rectangle((cx0, ry0), cx1 - cx0, ry1 - ry0, **rect_kwargs))
    ax.add_patch(Rectangle((nx0, ry0), nx1 - nx0, ry1 - ry0, **rect_kwargs))

    # Create inset axes on the right
    tumor_ax = inset_axes(ax, width="25%", height="25%", loc="center right", borderpad=1.2)
    healthy_ax = inset_axes(ax, width="25%", height="25%", loc="lower right", borderpad=1.2)

    data = g.data2d.values
    vmin, vmax = 0.0, 1.0
    cmap = plt.get_cmap("RdBu_r")

    tumor_ax.imshow(data[ry0:ry1, cx0:cx1], aspect="auto", interpolation="nearest", cmap=cmap, vmin=vmin, vmax=vmax)
    healthy_ax.imshow(data[ry0:ry1, nx0:nx1], aspect="auto", interpolation="nearest", cmap=cmap, vmin=vmin, vmax=vmax)
    for iax, title in [(tumor_ax, "Tumour samples"), (healthy_ax, "Healthy samples")]:
        iax.set_xticks([])
        iax.set_yticks([])
        iax.set_title(title, fontsize=9, pad=2)
        for spine in iax.spines.values():
            spine.set_linewidth(0.6)


def main() -> None:
    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    beta_df, meta = _read_figure1c()

    sns.set(context="paper", style="white")

    col_colors = _build_col_colors(meta)

    # Main clustered heatmap
    g = sns.clustermap(
        beta_df,
        method="average",
        metric="euclidean",
        cmap="RdBu_r",
        vmin=0.0,
        vmax=1.0,
        figsize=(12.8, 8.8),
        xticklabels=False,
        yticklabels=False,
        dendrogram_ratio=(0.08, 0.15),
        colors_ratio=(0.04, 0.04),  # thickness for col_colors bar block
        col_colors=col_colors,
        cbar_pos=(0.02, 0.25, 0.015, 0.35),  # (x, y, width, height) in figure coords
        cbar_kws={"label": "Fraction methylated"},
    )

    # Remove any title (paper has none)
    g.fig.suptitle("")

    # Label the row-annotation names on the left, matching the paper
    # Seaborn doesn't label the col_colors rows; we add text manually.
    cc_ax = g.ax_col_colors
    cc_ax.set_yticks([])
    cc_ax.set_xticks([])
    for spine in cc_ax.spines.values():
        spine.set_visible(False)
    # Add labels aligned with each annotation row
    labels = ["Sample group", "Tumour fraction", "fCpG discovery", "Sample type", "Platform"]
    # y positions in axes fraction coordinates (top -> bottom)
    y_fracs = np.linspace(0.9, 0.1, len(labels))
    for lab, y in zip(labels, y_fracs):
        cc_ax.text(
            -0.01,
            y,
            lab,
            transform=cc_ax.transAxes,
            ha="right",
            va="center",
            fontsize=9,
        )

    # Apply ordering for insets: meta needs to be reordered by the clustering
    col_order = g.dendrogram_col.reordered_ind
    ordered_cols = list(beta_df.columns[col_order])
    meta_reordered = meta.loc[ordered_cols]

    _add_legends(g.fig)
    _add_insets(g, meta_reordered)

    g.savefig(OUT_PDF, dpi=300, bbox_inches="tight")
    plt.close(g.fig)
    print(f"WROTE {OUT_PDF}")


if __name__ == "__main__":
    main()

