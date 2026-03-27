# Fig1（c–g）复现整理

## 基本信息
- 主要数据：`D:\naturedna\evoflux-reproduce\evoflux\41586_2025_9374_MOESM5_ESM.xlsx`
- 代码目录：`D:\naturedna\evoflux-reproduce\scripts`
- 输出目录：`D:\naturedna\evoflux-reproduce\figures\reproduced`

## 结果-代码-数据对应

| 结果文件 | 使用代码 | 使用数据 | 说明 |
|---|---|---|---|
| `Fig1_c.pdf` | `scripts/reproduce_fig1c_seaborn_annotated.py` | `MOESM5` sheet `Figure 1c` | 聚类热图（`average` + `euclidean` + `RdBu_r`），含顶部注释条、右侧图例、Tumour/Healthy 放大框。 |
| `Fig1_d.pdf` | `scripts/reproduce_fig1_c_to_g.R` | `MOESM5` sheet `Figure 1d` | 肿瘤/健康 fCpG β 值对比直方图，`binwidth=0.02`。 |
| `Fig1_e.pdf` | `scripts/reproduce_fig1e.py` | `MOESM5` sheet `Figure 1e` | log2 fold-change 热图，`RdBu_r`，左侧 colorbar。 |
| `Fig1_f.pdf` | `scripts/reproduce_fig1_c_to_g.R` | `MOESM5` sheet `Figure 1f` | 染色质状态热图（Healthy/MCL/CLL/MM）。 |
| `Fig1_g.pdf` | `scripts/reproduce_fig1_c_to_g.R` | `MOESM5` sheet `Figure 1g` | 按 `Data_source_Fig.1G` 的 panel 2–3 逻辑复现（log10 轴、Wilcoxon、beeswarm+boxplot）。 |

## 备注
- `Fig1_c.pdf` 已切换为 Python 注释增强版本（不是早期 pheatmap 版本）。
- `Fig1_g.pdf` 不包含 `pp_plot` 子面板（该部分所需独立数据不在 MOESM5 中）。
