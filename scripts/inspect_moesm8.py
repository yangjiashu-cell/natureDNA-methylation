"""Inspect MOESM8 Excel structure for Figure 4."""
from pathlib import Path

import pandas as pd

XLSX = Path(r"D:\naturedna\evoflux-reproduce\evoflux\41586_2025_9374_MOESM8_ESM.xlsx")
xl = pd.ExcelFile(XLSX)
lines = []
lines.append("所有 sheet 名称： " + str(xl.sheet_names))
for sheet in xl.sheet_names:
    df = pd.read_excel(xl, sheet_name=sheet, nrows=3)
    cols = df.columns.tolist()
    lines.append(f"Sheet '{sheet}' 列数={len(cols)} 前15列名： {cols[:15]}")
    if len(cols) > 15:
        lines.append(f"  ... 后5列： {cols[-5:]}")
out = "\n".join(lines)
print(out)
Path(r"D:\naturedna\evoflux-reproduce\figures\reproduced\moesm8_sheet_inspection.txt").write_text(
    out, encoding="utf-8"
)
