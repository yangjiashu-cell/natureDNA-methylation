"""
Fig. 1e — vertical heatmap matching Nature panel (MOESM5 + seaborn).
Stars: white, centered inside cells. Colorbar: left. cmap: RdBu_r [-3, 3].
"""
import argparse
import sys

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def main():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--xlsx",
        default=r"D:\naturedna\evoflux-reproduce\evoflux\41586_2025_9374_MOESM5_ESM.xlsx",
    )
    p.add_argument(
        "--out",
        default=r"D:\naturedna\evoflux-reproduce\figures\reproduced\Fig1_e.pdf",
    )
    args = p.parse_args()

    df = pd.read_excel(args.xlsx, sheet_name="Figure 1e", header=0)
    if "genomic_location" in df.columns:
        loc_col, val_col = "genomic_location", "log2_foldchange"
    else:
        loc_col, val_col = df.columns[0], df.columns[1]
    df_e = df.set_index(loc_col)[[val_col]].copy()
    df_e.columns = ["log2FC"]

    # Normalize index labels to paper spelling / order
    rename_idx = {
        "Open Sea": "Open sea",
        "open sea": "Open sea",
        "OPEN SEA": "Open sea",
    }
    df_e = df_e.rename(index=lambda x: rename_idx.get(str(x), str(x)))

    order = ["Island", "Open sea", "Shelf", "Shore"]
    missing = [k for k in order if k not in df_e.index]
    if missing:
        print(f"ERROR: missing rows after rename: {missing}, have: {list(df_e.index)}", file=sys.stderr)
        sys.exit(1)
    df_e = df_e.reindex(order)

    # Significance stars (paper figure); centered in cells, white
    stars = ["***", "***", "**", "***"]

    fig, ax = plt.subplots(figsize=(2.5, 6))
    hm = sns.heatmap(
        df_e,
        annot=False,
        cmap="RdBu_r",
        vmin=-3,
        vmax=3,
        linewidths=0,
        linecolor=None,
        cbar_kws={
            "label": r"$\log_2$(FC) of fCpGs vs non-fCpGs fractions",
            "location": "left",
        },
        ax=ax,
    )

    # White significance stars inside each cell (single column → x=0.5)
    n = len(df_e)
    for i, star in enumerate(stars):
        ax.text(0.5, i + 0.5, star, va="center", ha="center", fontsize=11, color="white", fontweight="bold")

    ax.set_xlabel("")
    ax.set_ylabel("Genomic location", labelpad=8)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    ax.set_xticks([])

    ax.set_title("e", fontsize=14, pad=10, loc="left", fontweight="bold")

    plt.tight_layout()
    fig.savefig(args.out, bbox_inches="tight")
    plt.close()
    print("Wrote", args.out)


if __name__ == "__main__":
    main()
