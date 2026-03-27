import pandas as pd

path = r"d:\naturedna\evoflux-reproduce\evoflux\41586_2025_9374_MOESM5_ESM.xlsx"

# 1g head
dfg = pd.read_excel(path, "Figure 1g", header=None, nrows=25)
print("Figure 1g raw head:")
print(dfg.to_string())

print("\n--- Figure 1f full ---")
dff = pd.read_excel(path, "Figure 1f")
print(dff.to_string())

print("\n--- Figure 1e full ---")
dfe = pd.read_excel(path, "Figure 1e")
print(dfe.to_string())

print("\n--- Figure 1d head ---")
dfd = pd.read_excel(path, "Figure 1d", nrows=5)
print(dfd.head())
print(dfd.columns)

print("\n--- Figure 1c corner ---")
dfc = pd.read_excel(path, "Figure 1c", nrows=3)
print(dfc.iloc[:, :5])
