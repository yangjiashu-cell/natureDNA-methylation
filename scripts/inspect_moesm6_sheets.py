import pandas as pd

xl = pd.ExcelFile(r"D:\naturedna\evoflux-reproduce\evoflux\41586_2025_9374_MOESM6_ESM.xlsx")
print("所有 sheet 名称：", xl.sheet_names)
for sheet in xl.sheet_names:
    df = pd.read_excel(xl, sheet_name=sheet, nrows=3)
    print(f"Sheet '{sheet}' 列名：", df.columns.tolist())
