# Fig5（a–c）复现整理

## 基本信息

- 代码文件：`scripts/reproduce_fig5.R`
- 参考代码：`Duran-FerrerM-evoflux/docs/Data_source_Fig.5.html`
- 数据文件：`evoflux/41586_2025_9374_MOESM9_ESM.xlsx`
- 输出目录：`figures/reproduced`

## 结果-代码-数据对应


| 结果文件        | 使用代码                       | 使用数据                                   | 说明                                       |
| ----------- | -------------------------- | -------------------------------------- | ---------------------------------------- |
| `Fig5a.pdf` | `scripts/reproduce_fig5.R` | `MOESM9`（`Figure 5a` + `Figure 5c` 合并） | 单变量 Cox（TTFT/OS）森林图。                     |
| `Fig5b.pdf` | `scripts/reproduce_fig5.R` | `MOESM9`（主要 `Figure 5a` 临床与 `theta`）   | IGHV 分层 + `theta` cutpoint 的 KM 曲线。      |
| `Fig5c.pdf` | `scripts/reproduce_fig5.R` | `MOESM9`（`Figure 5a` + `Figure 5c`）    | 多变量 Cox（`theta + IGHV + TP53 + Age`）森林图。 |
| `Fig5.pdf`  | `scripts/reproduce_fig5.R` | 同上                                     | a/b/c 组合版总图。                             |


