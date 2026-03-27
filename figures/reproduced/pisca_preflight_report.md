# PISCA 复现预检查报告

- 工作目录：`D:\naturedna\evoflux-reproduce`

## BEAST
- `BEAST_ROOT` 环境变量：`None`
- 发现 BEAST 根目录：`D:\naturedna\evoflux-reproduce\tools\BEAST v1.8.4`

## PISCA 插件
- 源码目录：`D:\naturedna\evoflux-reproduce\PISCA`
- dist jar 存在：`False`
- 已安装插件：`D:\naturedna\evoflux-reproduce\tools\BEAST v1.8.4\lib\plugins\PISCA.PISCALoader.jar`
- java in PATH：`True`
- ant in PATH：`True`

## 数据
- Fig4AB Data_source 文件：`D:\naturedna\evoflux-reproduce\Revision\Results\Longitudinal_meth_data\Data_source_methylation_Fig.4_SCLL12-SCLL19.xlsx`
- 是否存在：`False`
- MOESM8 文件：`D:\naturedna\evoflux-reproduce\evoflux\41586_2025_9374_MOESM8_ESM.xlsx`
- 是否存在：`True`

## 下一步必做
1. 补充作者 Fig4AB 原始输入：Revision/Results/Longitudinal_meth_data/Data_source_methylation_Fig.4_SCLL12-SCLL19.xlsx。
2. 若无法提供上述文件，可先用 MOESM8 构建替代输入（可跑流程，但与论文 Data_source 口径有差异）。
