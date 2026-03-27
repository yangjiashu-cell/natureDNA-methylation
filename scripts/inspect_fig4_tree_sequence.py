"""
Inspect MOESM8 Figure 4: sheet names, tree-sequence block, and clone-fraction block.

The Excel file does NOT contain separate sheets named
  "Figure 4a - Top (tree sequence)" / "Figure 4b - Top (tree sequence)".
Those strings appear as cell A0 *inside* sheets "Figure 4a" and "Figure 4b".

Run:
  python inspect_fig4_tree_sequence.py
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

XLSX = Path(__file__).resolve().parents[1] / "evoflux" / "41586_2025_9374_MOESM8_ESM.xlsx"


def main() -> None:
    xl = pd.ExcelFile(XLSX)
    print("=== All sheet names ===")
    print(xl.sheet_names)

    # header=None: MOESM8 row0 is a merged title ("Figure 4a - Top (tree sequence)").
    # Default read_excel(header=0) would wrongly promote row0 to column names.
    xls = pd.read_excel(XLSX, sheet_name=None, header=None)

    for sheet in ("Figure 4a", "Figure 4b"):
        df = xls[sheet]
        print(f"\n=== {sheet} ===")
        print(f"A0 (block title): {df.iloc[0, 0]!r}")
        hdr = df.iloc[1, 0:5].tolist()
        print("Row1 header [node .. mutations]:", hdr)
        tree_rows: list[dict] = []
        i = 2
        while i < len(df) and pd.notna(df.iloc[i, 0]):
            tree_rows.append(
                {
                    "node_label": df.iloc[i, 0],
                    "parent_node_label": df.iloc[i, 1],
                    "parent_node_row": df.iloc[i, 2],
                    "no.of.mutations.assigned": df.iloc[i, 3],
                }
            )
            i += 1
        print(pd.DataFrame(tree_rows).to_string(index=False))
        # Clone-info header is in row 1 cols 5-8 or 5-9
        print("Clone-info title (row0 col5):", df.iloc[0, 5])
        tcols = df.iloc[1, 5 : 5 + (5 if sheet == "Figure 4b" else 4)].tolist()
        print("Time columns:", tcols)


if __name__ == "__main__":
    main()
