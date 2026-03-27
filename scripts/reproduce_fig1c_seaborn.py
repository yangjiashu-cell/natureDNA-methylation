import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path


def main() -> None:
    file_path = Path(r"D:\naturedna\evoflux-reproduce\evoflux\41586_2025_9374_MOESM5_ESM.xlsx")
    out_path = Path(r"D:\naturedna\evoflux-reproduce\figures\reproduced\Fig1_c.pdf")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_excel(file_path, sheet_name="Figure 1c", header=None)

    # Auto-detect CpG matrix start row (first row whose first column starts with 'cg')
    start_row = None
    for i in range(15):
        v = df.iloc[i, 0]
        if isinstance(v, str) and v.startswith("cg"):
            start_row = i
            break
    if start_row is None:
        # Fallback: scan a larger window
        for i in range(min(200, len(df))):
            v = df.iloc[i, 0]
            if isinstance(v, str) and v.startswith("cg"):
                start_row = i
                break
    if start_row is None:
        raise RuntimeError("Could not find CpG start row (no 'cg*' in first column).")

    cpg_ids = df.iloc[start_row:, 0].astype(str).values
    sample_names = df.iloc[0, 1:].astype(str).values
    beta_mat = df.iloc[start_row:, 1:].astype(float).values

    beta_df = pd.DataFrame(beta_mat, index=cpg_ids, columns=sample_names)

    sns.set(context="paper", style="white")

    g = sns.clustermap(
        beta_df,
        method="average",
        metric="euclidean",
        cmap="RdBu_r",
        figsize=(22, 14),
        xticklabels=False,
        yticklabels=False,
        dendrogram_ratio=(0.08, 0.15),
        cbar_kws={"label": "β value (DNA methylation)"},
    )

    g.fig.suptitle(
        "Fig. 1c - 978 fCpGs hierarchical clustered heatmap (reproduced)",
        y=0.95,
    )

    g.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(g.fig)

    print(f"WROTE {out_path}")


if __name__ == "__main__":
    main()

