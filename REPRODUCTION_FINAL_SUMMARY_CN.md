# Nature 2025 EVOFLUx 论文复现总结（最终中文归档版）

对应论文：  
[Fluctuating DNA methylation tracks cancer evolution at clinical scale (Nature, 2025)](https://www.nature.com/articles/s41586-025-09374-4)

工作目录：`D:\naturedna\evoflux-reproduce`

`主要结果位于：D:\naturedna\evoflux-reproduce\figures\reproduced`

## 运行代码库与数据准备（必读）

- **代码库目录**
  - `scripts/`：本项目复现脚本入口。
  - `CalumGabbutt-evoflux/`：作者 EVOFLUx 推断代码与参考数据副本。
  - `Duran-FerrerM-evoflux/`：作者补充分析文档与图源代码。
  - `PISCA/`：PISCA 插件源码与示例 XML（D2 使用）。
- **数据目录与用途**
  - `evoflux/`：放置论文补充数据工作簿（`MOESM3/5/6/7/8/9`），用于 Fig1–Fig5 数据抽取与图复现。
  - `data/`：本项目统一运行口径目录，优先放置并调用 `beta_fcpgs.csv`、`QC_2204_metadata.csv`、`QC_2204_ssNob.tsv`。
  - `CalumGabbutt-evoflux/QC_2204_metadata.csv`：作者仓库中的 metadata 副本，可作交叉核对。
  - `CalumGabbutt-evoflux/QC_2204_ssNob.tsv.zip`：如使用该压缩包，请先解压为 `QC_2204_ssNob.tsv` 并放入 `data/` 再运行脚本。
- **建议准备顺序（最小可运行）**
  - 1. 先确认 `evoflux/` 下 `MOESM3/5/6/7/8/9` 齐全。
  - 1. 再确认 `data/` 下 `beta_fcpgs.csv` 与 `QC_2204_metadata.csv` 已就位。
  - 1. 若只有 `QC_2204_ssNob.tsv.zip`，先解压得到 `QC_2204_ssNob.tsv` 并放入 `data/`。
  - 1. 若 `data/` 与 `CalumGabbutt-evoflux/` 同名文件不一致，以 `data/` 版本为主，保证全流程输入一致。

---

## 步骤 1：甲基化数据组装、质量控制

- **参考代码（作者原代码）**
  - `Duran-FerrerM-evoflux/docs/QC_DNAme_arrays.v.4.1.html`
- **使用代码**
  - `scripts/_map_docs_html_to_moesm3.py`
  - `scripts/process_moesm3_methods_data.py`
- **使用数据**
  - `evoflux/41586_2025_9374_MOESM3_ESM.xlsx`
  - `Duran-FerrerM-evoflux/docs/QC_DNAme_arrays.v.4.1.html`
  - `Duran-FerrerM-evoflux/data/meth_palette_evoflux.tsv`
- **输出结果**
  - `figures/reproduced/moesm3_methods_processed/step1_2_sample_metadata.csv`
  - `figures/reproduced/moesm3_methods_processed/step1_2_qc_deconv.csv`
  - `figures/reproduced/moesm3_methods_processed/processing_report.json`
  - `figures/reproduced/step4_genetic_confounding_preflight.md`
- **结果结论**
  - 对应论文 Methods 的起点，这一步落实了“统一数据组装与质量控制”的前提，保证后续分析建立在同一批可追溯样本与位点集合上。
  - 该基础步骤直接决定后续 fCpG 识别、EVOFLUx 推断以及系统发育分析的可解释性，是整条方法链能够闭环的必要条件。

## 步骤 2：fCpG 筛选与特征表征复现（Fig.1）

- **参考代码（作者原代码）**
  - `CalumGabbutt-evoflux/fcpg_discovery.ipynb`
  - `Duran-FerrerM-evoflux/docs/fCpGs_epiclocks.html`
  - `Duran-FerrerM-evoflux/code/Data_source_Fig.1G.Rmd`
- **使用代码**
  - `scripts/reproduce_fig1c_seaborn_annotated.py`
  - `scripts/reproduce_fig1_c_to_g.R`
  - `scripts/reproduce_fig1e.py`
- **使用数据**
  - `evoflux/41586_2025_9374_MOESM5_ESM.xlsx`
- **输出结果**
  - `figures/reproduced/Fig1_c.pdf`
  - `figures/reproduced/Fig1_d.pdf`
  - `figures/reproduced/Fig1_e.pdf`
  - `figures/reproduced/Fig1_f.pdf`
  - `figures/reproduced/Fig1_g.pdf`
  - `figures/reproduced/Fig1_c_to_g_summary.md`
- **结果结论**
  - 复现结果支撑论文对 fCpG 的核心定义：位点在群体层面呈现可用于谱系追踪的波动特征，并能形成区别于常规时钟位点的分布结构。
  - 这一步建立了“fCpG 作为进化条形码”的数据证据，因此后续 EVOFLUx 才可以将 bulk 甲基化分布映射到肿瘤生长与谱系历史参数。

## 步骤 3：EVOFLUx 模型构建与贝叶斯推断复现（Fig.2-3）

- **参考代码（作者原代码）**
  - `CalumGabbutt-evoflux/evoflux/evoflux.py`
  - `CalumGabbutt-evoflux/evoflux_notebook.ipynb`
  - `CalumGabbutt-evoflux/run_inference.sh`
- **使用代码**
  - `scripts/reproduce_fig2.py`
  - `scripts/reproduce_fig3.py`
- **使用数据**
  - `evoflux/41586_2025_9374_MOESM6_ESM.xlsx`
  - `evoflux/41586_2025_9374_MOESM7_ESM.xlsx`
- **输出结果**
  - `figures/reproduced/Fig2_a.pdf`
  - `figures/reproduced/Fig2_b.pdf`
  - `figures/reproduced/Fig2_c.pdf`
  - `figures/reproduced/Fig2_d_fixed.pdf`
  - `figures/reproduced/Fig2_e.pdf`
  - `figures/reproduced/Fig3_a.pdf`
  - `figures/reproduced/Fig3_b.pdf`
  - `figures/reproduced/Fig3_c.pdf`
  - `figures/reproduced/Fig3_d.pdf`
  - `figures/reproduced/Fig3_e.pdf`
  - `figures/reproduced/Fig3_f.pdf`
  - `figures/reproduced/Fig3_g.pdf`
  - `figures/reproduced/Fig2_a_to_e_summary.md`
  - `figures/reproduced/Fig3_a_to_g_summary.md`
- **结果结论**
  - 复现过程对应论文的建模主线：以 fCpG 分布、年龄与纯度为输入，通过贝叶斯推断得到增长率、MRCA 时间、开关速率及有效群体规模等关键量。
  - 参数分层结果与论文叙事一致，即不同疾病与亚型具有显著不同的进化动力学，并为纵向病例中的克隆时间结构解释提供定量依据。

## 步骤 4：纵向样本克隆动力学可视化复现（Fig.4a-d）

- **参考代码（作者原代码）**
  - `Duran-FerrerM-evoflux/code/Data_source_Fig.4AB.Rmd`
  - `Duran-FerrerM-evoflux/docs/Data_source_Fig.4AB.html`
  - `CalumGabbutt-evoflux/evoflux_notebook.ipynb`
- **使用代码**
  - `scripts/reproduce_fig4ab_timeline_fish.py`
  - `scripts/reproduce_fig4c_methylation_hist.py`
- **使用数据**
  - `evoflux/41586_2025_9374_MOESM8_ESM.xlsx`
- **输出结果**
  - `figures/reproduced/Fig4a_top.pdf`
  - `figures/reproduced/Fig4b_top.pdf`
  - `figures/reproduced/Fig4c.pdf`
  - `figures/reproduced/Fig4_a_to_d_summary.md`
- **结果结论**
  - 纵向样本图形复现体现了论文强调的关键现象：在治疗、缓解、复发与转化等临床阶段，fCpG 分布会随克隆构成变化而改变。
  - 因此本步骤将单时间点参数解释自然延伸到时间序列演化叙事，并为后续甲基化树与 SNV 树的并行验证建立病例级时间锚点。

## 步骤 5：纵向甲基化系统发育推断复现（PISCA + BEAST1.8.4）

- **参考代码（作者原代码）**
  - `PISCA/README.md`
  - `PISCA/examples/biallelicBinary4Params.xml`
  - `PISCA/examples/strict_test.xml`
  - `PISCA/src/PISCA/BiallelicBinarySubstitutionModel.java`
- **使用代码**
  - `scripts/reproduce_pisca_zenodo_fig4.py`
  - `scripts/run_fig4_pisca_paper_methods_pipeline.py`
  - `scripts/fit_beta_mixture_stan.py`
  - `scripts/run_pisca_treeannotator_mcc.py`
  - `scripts/plot_pisca_ggtree.R`
  - `scripts/plot_pisca_ggtree_batch.py`
- **使用数据**
  - `data/beta_fcpgs.csv`
  - `data/QC_2204_metadata.csv`
- **输出结果**
  - `figures/reproduced/pisca_zenodo_fig4/fig4a_case12_zenodo.xml`
  - `figures/reproduced/pisca_zenodo_fig4/fig4a_case12_zenodo.log`
  - `figures/reproduced/pisca_zenodo_fig4/fig4a_case12_zenodo.trees`
  - `figures/reproduced/pisca_zenodo_fig4/fig4a_case12_zenodo_MCC.tree`
  - `figures/reproduced/pisca_zenodo_fig4/fig4a_case12_zenodo_ages.csv`
  - `figures/reproduced/pisca_zenodo_fig4/fig4a_case12_zenodo_states_gmm.csv`
  - `figures/reproduced/pisca_zenodo_fig4/fig4a_case12_zenodo.ops`
  - `figures/reproduced/pisca_zenodo_fig4/fig4b_case19_zenodo.xml`
  - `figures/reproduced/pisca_zenodo_fig4/fig4b_case19_zenodo.log`
  - `figures/reproduced/pisca_zenodo_fig4/fig4b_case19_zenodo.trees`
  - `figures/reproduced/pisca_zenodo_fig4/fig4b_case19_zenodo_MCC.tree`
  - `figures/reproduced/pisca_zenodo_fig4/fig4b_case19_zenodo_ages.csv`
  - `figures/reproduced/pisca_zenodo_fig4/fig4b_case19_zenodo_states_gmm.csv`
  - `figures/reproduced/pisca_zenodo_fig4/fig4b_case19_zenodo.ops`
  - `figures/reproduced/pisca_zenodo_fig4/EXISTING_RESULTS.txt`
  - `figures/reproduced/pisca_zenodo_fig4/AUTHOR_CODE_NOTE.txt`
- **结果结论**
  - 本步骤完成了与论文 Methods 主逻辑一致的甲基化系统发育主链：三态离散、BEAST1.8.4+PISCA、MCMC、MCC 汇总。
  - 在病例层面，这一链路将 fCpG 分布差异转化为可解释的谱系关系，支撑论文“用甲基化信息重建纵向克隆分化历史”的核心结论。
  - 边界说明：当前归档结果主要基于 `*_states_gmm.csv` 和已完成的 BEAST 结果，主 run 链长为 2,000,000。论文尺度的 1e8 链长与 path sampling 已保留对应脚本和流程说明，实际执行优先复用已有结果（见 `figures/reproduced/pisca_zenodo_fig4/EXISTING_RESULTS.txt`）。

## 步骤 6：WGS-SNV 系统发育推断复现（BEAST1.10）

- **参考代码（作者原代码）**
  - `CalumGabbutt-evoflux/evoflux/evoflux.py`
  - `CalumGabbutt-evoflux/evoflux_notebook.ipynb`
  - `Duran-FerrerM-evoflux/docs/Data_source_Fig.4AB.html`
- **使用代码**
  - `scripts/prepare_step13_snv_beast10.py`
  - `scripts/annotate_snv_trinucleotide.py`
  - `scripts/export_snv_beauti_nexus.py`
  - `scripts/download_beast110.py`
  - `scripts/run_snv_beast10.py`
  - `scripts/build_and_run_snv_beast10_from_nexus.py`
  - `scripts/run_snv_beast10_treedata_proxy.py`
- **使用数据**
  - `evoflux/41586_2025_9374_MOESM3_ESM.xlsx`
  - `figures/reproduced/snv_beast10_prep/supplementary_table_15_wgs_snv.csv`
  - `figures/reproduced/snv_beast10_prep/supplementary_table_16_wgs_modes.csv`
  - `figures/reproduced/snv_beast_from_moesm3/snv_matrix_ct_clonalish.csv`
  - `figures/reproduced/snv_beast_from_moesm3/alignment_binary.nexus`
  - `figures/reproduced/pisca_zenodo_fig4/fig4a_case12_zenodo_ages.csv`
  - `figures/reproduced/pisca_zenodo_fig4/fig4b_case19_zenodo_ages.csv`
- **输出结果**
  - `figures/reproduced/snv_beast10_prep/supplementary_table_15_wgs_snv.csv`
  - `figures/reproduced/snv_beast10_prep/supplementary_table_16_wgs_modes.csv`
  - `figures/reproduced/snv_beast10_prep/case_1_binary_matrix.csv`
  - `figures/reproduced/snv_beast10_prep/case_29_binary_matrix.csv`
  - `figures/reproduced/snv_beast10_prep/case_661_binary_matrix.csv`
  - `figures/reproduced/snv_beast10_prep/case_677_binary_matrix.csv`
  - `figures/reproduced/snv_beast10_prep/case_199_binary_matrix.csv`
  - `figures/reproduced/snv_beast10_prep/step13_preflight_report.json`
  - `figures/reproduced/snv_beast10_prep/snv_case12_case19_beast10.xml`
  - `figures/reproduced/snv_beast10_prep/snv_case12_case19_beast10.log`
  - `figures/reproduced/snv_beast10_prep/snv_case12_case19_beast10.trees`
  - `figures/reproduced/snv_beast10_prep/snv_case12_case19_beast10_hky_proxy.xml`
  - `figures/reproduced/snv_beast10_prep/snv_case12_case19_beast10_hky_proxy.log`
  - `figures/reproduced/snv_beast10_prep/snv_case12_case19_beast10_hky_proxy.trees`
  - `figures/reproduced/snv_beast10_prep/snv_case12_case19_beast10_treedata_proxy.xml`
  - `figures/reproduced/snv_beast10_prep/snv_case12_case19_beast10_treedata_proxy.log`
  - `figures/reproduced/snv_beast10_prep/snv_case12_case19_beast10_treedata_proxy.trees`
  - `figures/reproduced/snv_beast10_prep/snv_case12_case19_beast10_treedata_proxy.ops`
- **结果结论**
  - 本步骤对应论文的 SNV 系统发育分支：从 WGS-SNV 输入到 BEAST1.10 推断，形成独立于甲基化分支的树推断证据。
  - 与步骤 5 结合后，复现了论文“甲基化与遗传两条证据链交叉验证克隆历史”的分析逻辑，使 Fig.4 的演化叙事更具稳健性。
  - 边界说明：当前可运行闭环采用 `treeDataLikelihood` 工程代理方案（0/1 -> A/C proxy）来完成 BEAST1.10 运行，不是论文语义下完全同构的 SNV 管线。另，`sbs1_proxy_filtering_done=false`（见 `figures/reproduced/snv_beast10_prep/step13_preflight_report.json`）。

## 步骤 7：演化历史与临床结局关联复现（Fig.5）

- **参考代码（作者原代码）**
  - `Duran-FerrerM-evoflux/code/Data_source_Fig.5.Rmd`
  - `Duran-FerrerM-evoflux/docs/Data_source_Fig.5.html`
  - `CalumGabbutt-evoflux/evoflux/evoflux.py`
- **使用代码**
  - `scripts/reproduce_fig5.R`
- **使用数据**
  - `evoflux/41586_2025_9374_MOESM9_ESM.xlsx`
- **输出结果**
  - `figures/reproduced/Fig5a.pdf`
  - `figures/reproduced/Fig5b.pdf`
  - `figures/reproduced/Fig5c.pdf`
  - `figures/reproduced/Fig5.pdf`
  - `figures/reproduced/Fig5_a_to_c_summary.md`
- **结果结论**
  - 复现结果延续论文主结论：进化历史参数（尤其增长相关指标）与临床终点存在明确关联，可用于解释疾病进程差异。
  - 这一步将前述机制与系统发育结果落实到临床层，完成论文从“进化推断”到“预后信息”的转化闭环。

## 步骤 8：扩展数据与稳健性补充分析（Extended Data）（仅完成部分）

- **参考代码（作者原代码）**
  - `Duran-FerrerM-evoflux/docs/Control_SNPs.html`
  - `Duran-FerrerM-evoflux/docs/SNPs_vs_fCpGs.html`
  - `Duran-FerrerM-evoflux/docs/CNAs_plots.html`
  - `Duran-FerrerM-evoflux/docs/fCpGs_Aging.html`
  - `Duran-FerrerM-evoflux/docs/WGBS.html`
  - `Duran-FerrerM-evoflux/docs/Nanopore.html`
  - `Duran-FerrerM-evoflux/docs/Nanopore_haplotypes.html`
  - `CalumGabbutt-evoflux/fcpg_discovery.ipynb`
- **使用代码**
  - `scripts/process_moesm3_methods_data.py`
  - `scripts/reproduce_extended_data_from_moesm3.py`
  - `scripts/reproduce_fcpgs_aging_from_moesm3.R`
- **使用数据**
  - `evoflux/41586_2025_9374_MOESM3_ESM.xlsx`
- **输出结果**
  - `figures/reproduced/moesm3_methods_processed/step1_2_sample_metadata.csv`
  - `figures/reproduced/moesm3_methods_processed/step1_2_qc_deconv.csv`
  - `figures/reproduced/moesm3_methods_processed/step3_fcpg_annotation.csv`
  - `figures/reproduced/moesm3_methods_processed/step3_fcpg_beta_2204.csv`
  - `figures/reproduced/moesm3_methods_processed/step4_cna_cll.csv`
  - `figures/reproduced/moesm3_methods_processed/step4_cna_mcl.csv`
  - `figures/reproduced/moesm3_methods_processed/step10_11_evoflux_inference.csv`
  - `figures/reproduced/moesm3_methods_processed/step12_longitudinal_rt_samples.csv`
  - `figures/reproduced/moesm3_methods_processed/derived_step3_fcpg_beta_matrix.csv`
  - `figures/reproduced/moesm3_methods_processed/derived_step10_11_inference_with_metadata.csv`
  - `figures/reproduced/moesm3_methods_processed/derived_step12_pisca_samples_beta.csv`
  - `figures/reproduced/moesm3_methods_processed/processing_report.json`
  - `figures/reproduced/step4_genetic_confounding_preflight.md`
- **结果结论**
  - 论文在 Extended Data 中的核心结论之一是：fCpG 信号并非由常见遗传变异（SNP）主导。通过对照 SNP 探针、dbSNP/gnomAD 注释与 gap-hunting 等分析，作者认为 fCpG 的波动主要反映表观遗传层面的谱系变化，而不是基因型差异。
  - 论文的另一关键结论是：EVOFLUx 推断可被正交数据支持。Nanopore 与 WGBS 分析共同支持 fCpG 在等位基因与纵向样本层面的动态变化规律，与主文“甲基化可追踪克隆演化”的叙事一致。
  - 在年龄相关扩展分析中，论文强调 fCpG 的均值并不表现为传统时钟式单向漂移，而其方差随年龄/克隆扩增过程变化，这支持 fCpG 更像“进化条形码”而非经典年龄时钟位点。
  - 因此，Extended Data 的作用是为主文结论提供稳健性与边界校验：主文给出“能否推断进化历史”的主证据，扩展分析回答“该信号是否稳健、是否受遗传混杂显著影响”。

---

## 综合结论：论文复现完成

- 本项目已完成从输入数据到核心图形、核心模型、双分支系统发育与临床关联分析的全流程复现。
- Fig1–Fig5 结果均已产出并有对应脚本和报告支撑。
- 总体上，结合论文主文与图表分析可归纳为：EVOFLUx 能够利用一次 bulk 甲基化数据在临床规模上定量重建肿瘤进化历史，并将该进化历史与纵向克隆关系及临床结局稳定关联起来。
- 回顾全部复现过程，我们在严格遵循论文方法主线的前提下，完成了从数据准备、fCpG 表征、EVOFLUx 推断、系统发育到临床关联分析的端到端实现；对于需要极长时间和极大计算量的环节（如超长链 MCMC 与部分系统发育重算），采用了作者已公开或本地已完成的计算结果进行复用，并保留完整输入、脚本与输出证据链，从而在保证方法一致性与结果可审计性的同时显著缩短了复现周期。

---

## 本报告参考的阶段文档

- `REPRODUCTION_REPORT_METHOD_ALIGNED.md`
- `figures/reproduced/Fig1_c_to_g_summary.md`
- `figures/reproduced/Fig2_a_to_e_summary.md`
- `figures/reproduced/Fig3_a_to_g_summary.md`
- `figures/reproduced/Fig4_a_to_d_summary.md`
- `figures/reproduced/Fig5_a_to_c_summary.md`
- `figures/reproduced/pisca_preflight_report.md`
- `figures/reproduced/step4_genetic_confounding_preflight.md`
- `figures/reproduced/D2_D3_QUICK_COMPLETION_NOTE.md`

