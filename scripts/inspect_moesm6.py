import pandas as pd

p = r"d:\naturedna\evoflux-reproduce\evoflux\41586_2025_9374_MOESM6_ESM.xlsx"
x = pd.ExcelFile(p)
for s in x.sheet_names:
    df = pd.read_excel(x, s, nrows=8)
    n = len(pd.read_excel(x, s))
    print("===", s, "rows=", n, "cols=", list(df.columns))
    print(df.head(4).to_string())
    print()
