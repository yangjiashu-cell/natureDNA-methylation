# Fig2（a–e）复现整理

## 基本信息
- 代码文件：`scripts/reproduce_fig2.py`
- 数据文件：`evoflux/41586_2025_9374_MOESM6_ESM.xlsx`
- 输出目录：`figures/reproduced`

## 结果-代码-数据对应
| 结果文件 | 使用代码 | 使用数据 | 说明 |
|---|---|---|---|
| `Fig2_a.pdf` | `scripts/reproduce_fig2.py` | `MOESM6`（示意面板） | EVOFLUx 方法学示意图。 |
| `Fig2_b.pdf` | `scripts/reproduce_fig2.py` | `MOESM6` sheet `Figure 2b` | recent vs distant MRCA 分布对比。 |
| `Fig2_c.pdf` | `scripts/reproduce_fig2.py` | `MOESM6` sheet `Figure 2c` | high vs low growth 分布对比。 |
| `Fig2_d_fixed.pdf` | `scripts/reproduce_fig2.py` | `MOESM6` sheet `Figure 2d`（`header=1`） | 对称 fish 图 + T1/T2/T3/T4 散点对比。 |
| `Fig2_e.pdf` | `scripts/reproduce_fig2.py` | `MOESM6` sheet `Figure 2e`（`header=1`） | T1/T3 直方图叠加 + intermediate 比例柱 + P 值。 |

## 备注
- d 面板当前定稿文件为 `Fig2_d_fixed.pdf`。
