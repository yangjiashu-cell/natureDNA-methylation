"""Inspect MOESM5 Fig.1 source data workbook."""
import pandas as pd

path = r"d:\naturedna\evoflux-reproduce\evoflux\41586_2025_9374_MOESM5_ESM.xlsx"
x = pd.ExcelFile(path)
print("SHEETS:", x.sheet_names)
for s in x.sheet_names:
    df = pd.read_excel(x, s, nrows=5)
    n = len(pd.read_excel(x, s))
    print(f"\n=== {s} === rows={n} cols={len(df.columns)}")
    print(list(df.columns)[:20])
