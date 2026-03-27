"""Step-1 MOESM6 sheet inspection (UTF-8)."""
from pathlib import Path

import pandas as pd

xl = pd.ExcelFile(r"D:\naturedna\evoflux-reproduce\evoflux\41586_2025_9374_MOESM6_ESM.xlsx")
lines = []
lines.append("所有 sheet 名称： " + str(xl.sheet_names))
for sheet in xl.sheet_names:
    df = pd.read_excel(xl, sheet_name=sheet, nrows=5)
    lines.append(f"Sheet '{sheet}' 前5列名： {df.columns.tolist()}")
out = "\n".join(lines)
print(out)
Path(r"D:\naturedna\evoflux-reproduce\figures\reproduced\moesm6_step1_sheet_inspection.txt").write_text(
    out, encoding="utf-8"
)
