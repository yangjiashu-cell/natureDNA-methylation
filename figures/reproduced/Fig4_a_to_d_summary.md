# Fig4（a–d）复现整理

## 基本信息

- 代码文件：`scripts/reproduce_fig4ab_timeline_fish.py`、`scripts/reproduce_fig4c_methylation_hist.py`、`scripts/reproduce_fig4_bottom_methylation_heatmap.R`
- 数据文件：`evoflux/41586_2025_9374_MOESM8_ESM.xlsx`
- 输出目录：`figures/reproduced`

## 结果-代码-数据对应


| 结果文件                             | 使用代码                                                  | 使用数据                       | 说明                                |
| -------------------------------- | ----------------------------------------------------- | -------------------------- | --------------------------------- |
| `Fig4a_top.pdf`                  | `scripts/reproduce_fig4ab_timeline_fish.py`           | `MOESM8` sheet `Figure 4a` | SCLL-012 顶部 timeline + fish plot。 |
| `Fig4b_top.pdf`                  | `scripts/reproduce_fig4ab_timeline_fish.py`           | `MOESM8` sheet `Figure 4b` | SCLL-019 顶部 timeline + fish plot。 |
| `Fig4c.pdf`                      | `scripts/reproduce_fig4c_methylation_hist.py`         | `MOESM8` sheet `Figure 4c` | 三时间点甲基化条形直方图。                     |
| `Fig4d.pdf`                      | `scripts/reproduce_fig4_bottom_methylation_heatmap.R` | `MOESM8` sheet `Figure 4d` | 两时间点甲基化分布图。                       |
| `Fig4a_bottom_methylation_R.pdf` | `scripts/reproduce_fig4_bottom_methylation_heatmap.R` | `MOESM8`                   | Fig4a 底部甲基化热图。                    |
| `Fig4b_bottom_methylation_R.pdf` | `scripts/reproduce_fig4_bottom_methylation_heatmap.R` | `MOESM8`                   | Fig4b 底部甲基化热图。                    |
| `Fig4a_scatter.pdf`              | `scripts/reproduce_fig4_bottom_methylation_heatmap.R` | `MOESM8`                   | Fig4a 底部散点对比。                     |
| `Fig4b_scatter.pdf`              | `scripts/reproduce_fig4_bottom_methylation_heatmap.R` | `MOESM8`                   | Fig4b 底部散点对比。                     |


## Fig.4 底部系统发育（PISCA / BEAST 1.8.4）

### Zenodo 主矩阵（推荐，对齐 harmonised 数据）

- **代码**：`scripts/reproduce_pisca_zenodo_fig4.py`（可选 `scripts/run_pisca_treeannotator_mcc.py` 单独重跑 MCC）
- **数据**：`data/beta_fcpgs.csv`（978 fCpG × 样本）、`data/QC_2204_metadata.csv`（`AGE_SAMPLING`、样本 ID）
- **输出目录**：`figures/reproduced/pisca_zenodo_fig4/`（若已有 `*.log`/`*_MCC.tree`，可只跑 ggtree 节省时间，见 `EXISTING_RESULTS.txt`）
- **BEAST 1.10（SNV 用）**：`scripts/download_beast110.py` 从 GitHub 官方 release 拉取 **v1.10.4** 到 `tools/`。
- **批量 ggtree PDF**：`scripts/plot_pisca_ggtree_batch.py`（依赖 R + ape/ggtree/ggplot2）。
- **默认样本集**：Case12 = `SCLL-545`–`548`（T1–T4）；Case19 = `SCLL-493,494,531,532,533`（T1–T5，不含与 `531` 近重复的 `018`）
- **说明**：论文 Methods 为 **Stan 三态 Beta 混合** 离散化。作者仓库 **CalumGabbutt/evoflux** 仅有 EVOFLUx 推断代码，无 PISCA/BEAST/Stan 脚本；说明见 `figures/reproduced/pisca_zenodo_fig4/AUTHOR_CODE_NOTE.txt`。推荐一键流程：`scripts/run_fig4_pisca_paper_methods_pipeline.py`（Stan→CSV→`--discretize external`）。亦可 `--discretize stan`（需 CmdStan）或分步 `fit_beta_mixture_stan.py --case … --out-csv`。`--paper-methods` = Stan + paper 链长（1e8）+ path sampling 说明；`--also-write-exponential-xml` 另存指数增长群体模型 XML 供边际似然/path sampling 与常数模型对照。未指定时默认 **GMM**；**beta_em** 为同模型 EM 近似。先验按 Methods。MCC：`10% burn-in` + **median** 节点高。ggtree：`scripts/plot_pisca_ggtree.R`（需 ape、ggtree、ggplot2）。
- **已生成 ggtree PDF（基于现有 MCC，免重跑 BEAST）**：
  - `figures/reproduced/pisca_zenodo_fig4/fig4a_case12_zenodo_MCC_ggtree.pdf`
  - `figures/reproduced/pisca_zenodo_fig4/fig4b_case19_zenodo_MCC_ggtree.pdf`



