# Step4 Genetic Confounding Preflight

## Scope

- HTML workflows checked: Control_SNPs, SNPs_vs_fCpGs, CNAs_plots, Data_source_Fig.4AB, Nanopore, Nanopore_haplotypes.
- Data search roots: `CalumGabbutt-evoflux/data`, `data`, `evoflux`, `Duran-FerrerM-evoflux/data`.

## Summary

- Total referenced external data paths: **20**
- Directly resolvable from HTML relative paths: **0**
- Missing direct paths but filename exists somewhere in allowed roots: **4**
- Unresolved missing dependencies: **16**

## Dependency Matrix


| HTML                       | Referenced path                                                                                          | Status           | Resolved candidate                                                                   |
| -------------------------- | -------------------------------------------------------------------------------------------------------- | ---------------- | ------------------------------------------------------------------------------------ |
| `Control_SNPs.html`        | `../Data/FMC/Final_data/beta_filtered_samples_fcpgs_laplacian.tsv`                                       | `missing`        | `-`                                                                                  |
| `Control_SNPs.html`        | `../Revision/Data/Supplementary_Tables_revision_curated.3.xlsx`                                          | `missing`        | `-`                                                                                  |
| `Control_SNPs.html`        | `../Revision/Data/meth_palette_evoflux.tsv`                                                              | `alt_name_match` | `D:\naturedna\evoflux-reproduce\Duran-FerrerM-evoflux\data\meth_palette_evoflux.tsv` |
| `SNPs_vs_fCpGs.html`       | `../Data/FMC/Final_data/beta_filtered_samples_fcpgs_laplacian.tsv`                                       | `missing`        | `-`                                                                                  |
| `SNPs_vs_fCpGs.html`       | `../Data/FMC/Final_data/candidate_cpgs_24_07_23.txt`                                                     | `missing`        | `-`                                                                                  |
| `SNPs_vs_fCpGs.html`       | `../Results/fCpGs_SNPs/fCpGs.complete.info.xlsx`                                                         | `missing`        | `-`                                                                                  |
| `SNPs_vs_fCpGs.html`       | `../Revision/Data/Supplementary_Tables_revision_curated.3.xlsx`                                          | `missing`        | `-`                                                                                  |
| `SNPs_vs_fCpGs.html`       | `../Revision/Data/meth_palette_evoflux.tsv`                                                              | `alt_name_match` | `D:\naturedna\evoflux-reproduce\Duran-FerrerM-evoflux\data\meth_palette_evoflux.tsv` |
| `CNAs_plots.html`          | `../Manuscript/Submission/SUBMITTED_NATURE_171023/submission/SupplementaryTables.xlsx`                   | `missing`        | `-`                                                                                  |
| `CNAs_plots.html`          | `../Revision/Data/CNAs/CLLs_ICGC_CNAs_fCpGs_03_02_24.xlsx`                                               | `missing`        | `-`                                                                                  |
| `CNAs_plots.html`          | `../Revision/Data/CNAs/MCL_CNAs_merged_fCpGs_04_04_24.xlsx`                                              | `missing`        | `-`                                                                                  |
| `Data_source_Fig.4AB.html` | `../../Revision/Results/Longitudinal_meth_data/Data_source_methylation_Fig.4_SCLL12-SCLL19.xlsx`         | `missing`        | `-`                                                                                  |
| `Nanopore.html`            | `../Data/FMC/Final_data/beta_filtered_samples_fcpgs_laplacian.tsv`                                       | `missing`        | `-`                                                                                  |
| `Nanopore.html`            | `../Revision/Data/Nanopore/Nanopore_CLL_Bcell_curated.xlsx`                                              | `missing`        | `-`                                                                                  |
| `Nanopore.html`            | `../Revision/Data/Nanopore/beds/BCLLatlas_5mC_CpG_meth_cov.txt.gz`                                       | `missing`        | `-`                                                                                  |
| `Nanopore.html`            | `../Revision/Data/Supplementary_Tables_revision_curated.2.xlsx`                                          | `missing`        | `-`                                                                                  |
| `Nanopore.html`            | `../Revision/Data/meth_palette_evoflux.tsv`                                                              | `alt_name_match` | `D:\naturedna\evoflux-reproduce\Duran-FerrerM-evoflux\data\meth_palette_evoflux.tsv` |
| `Nanopore_haplotypes.html` | `../Revision//Data/meth_palette_evoflux.tsv`                                                             | `alt_name_match` | `D:\naturedna\evoflux-reproduce\Duran-FerrerM-evoflux\data\meth_palette_evoflux.tsv` |
| `Nanopore_haplotypes.html` | `../Revision/Results/fcpgs_Nanopore/Fig1_haplotypes/combined_reads_haplotype_chr7_25854120-25856220.csv` | `missing`        | `-`                                                                                  |
| `Nanopore_haplotypes.html` | `../Revision/Results/fcpgs_Nanopore/Fig1_haplotypes/readinfo_chr7_25854120-25856220.csv`                 | `missing`        | `-`                                                                                  |


## Execution Readiness

- End-to-end execution from provided HTML code is **not** currently possible due to unresolved inputs.
- Also requires an R runtime with packages (`data.table`, `openxlsx`, `minfi`, annotation packages).

