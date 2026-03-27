#!/usr/bin/env Rscript
## Nature 2025 Fig. 1 panels c–g from MOESM5 source data only.

suppressPackageStartupMessages({
  if (!requireNamespace("readxl", quietly = TRUE))
    install.packages("readxl", repos = "https://cloud.r-project.org", type = "win.binary")
  if (!requireNamespace("ggplot2", quietly = TRUE))
    install.packages("ggplot2", repos = "https://cloud.r-project.org", type = "win.binary")
  if (!requireNamespace("pheatmap", quietly = TRUE))
    install.packages("pheatmap", repos = "https://cloud.r-project.org", type = "win.binary")
  if (!requireNamespace("RColorBrewer", quietly = TRUE))
    install.packages("RColorBrewer", repos = "https://cloud.r-project.org", type = "win.binary")
  if (!requireNamespace("reshape2", quietly = TRUE))
    install.packages("reshape2", repos = "https://cloud.r-project.org", type = "win.binary")
  if (!requireNamespace("patchwork", quietly = TRUE))
    install.packages("patchwork", repos = "https://cloud.r-project.org", type = "win.binary")
  for (pkg in c("ggpubr", "ggbeeswarm", "cowplot", "scales", "dplyr")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      install.packages(pkg, repos = "https://cloud.r-project.org", type = "win.binary")
  }
  library(readxl)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(reshape2)
  library(patchwork)
  library(ggpubr)
  library(ggbeeswarm)
  library(cowplot)
  library(scales)
  library(dplyr)
  library(grid)
})

xlsx <- normalizePath("d:/naturedna/evoflux-reproduce/evoflux/41586_2025_9374_MOESM5_ESM.xlsx", winslash = "/", mustWork = TRUE)
outdir <- normalizePath("d:/naturedna/evoflux-reproduce/figures/reproduced", winslash = "/", mustWork = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------------------- #
# Fig. 1c — hierarchical clustering heatmap (978 fCpGs × n samples)
# Data: sheet "Figure 1c" — row 1 = sample IDs, row 2 = "Sample group", rows 7+ = β per probe
# --------------------------------------------------------------------------- #
message("Reading Figure 1c …")
raw_c <- as.data.frame(readxl::read_excel(xlsx, sheet = "Figure 1c", col_names = FALSE, .name_repair = "minimal"))
sids <- as.character(unlist(raw_c[1, -1, drop = TRUE]))
grp <- as.character(unlist(raw_c[2, -1, drop = TRUE]))
mat_c <- raw_c[7:nrow(raw_c), -1, drop = FALSE]
rownames(mat_c) <- as.character(unlist(raw_c[7:nrow(raw_c), 1]))
mat_c <- as.matrix(apply(mat_c, 2, function(x) as.numeric(as.character(x))))
colnames(mat_c) <- sids
storage.mode(mat_c) <- "double"

ann_col <- data.frame(Sample_group = grp, row.names = sids)
# Use the paper palette from docs/SNPs_vs_fCpGs.html (cols.paper) when possible,
# to keep consistent colouring across figures.
cols_paper <- c(
  "NA_values" = "#666666",
  "HPC" = "#E6E6E6",
  "pre.Bcell" = "#D5D5D5",
  "NBC" = "#C3C3C3",
  "GCB" = "#AEAEAE",
  "tPC" = "#969696",
  "MBC" = "#787878",
  "bmPC" = "#4D4D4D",
  "Tcell" = "#F0F8FF",
  "PBMCs" = "#FFF8DC",
  "Whole_blood" = "#CDC8B1",
  "Normal_lymphoid_cell" = "#D5D5D5",
  "Normal_myeloid_cell" = "#E6E6FA",
  "T-ALL" = "#CD3333",
  "T-ALL-remission" = "#C7A5A5",
  "B-ALL" = "#E69F00",
  "B-ALL-remission" = "#D8CAAB",
  "MCL" = "#F0E442",
  "C1.or.cMCL" = "#F5D51B",
  "C2.or.nnMCL" = "#0072B2",
  "DLBCL-NOS" = "#56B4E9",
  "MBL" = "#5A9A89",
  "CLL" = "#009E73",
  "RT" = "#02C590",
  "CLL.unmutated" = "#e65100",
  "CLL.mutated" = "#361379",
  "nCLL" = "#006E93",
  "iCLL" = "#FDC010",
  "mCLL" = "#963736",
  "Disease.Unclassified" = "#CCCCCC",
  "MGUS" = "#BF99AE",
  "MM" = "#CC79A7"
)

uniq_grp <- sort(unique(grp))
fallback_cols <- setNames(
  colorRampPalette(brewer.pal(8, "Set3"))(length(uniq_grp)),
  uniq_grp
)
group_colors <- setNames(fallback_cols[uniq_grp], uniq_grp)
for (g in uniq_grp) {
  if (!is.na(cols_paper[g])) group_colors[g] <- cols_paper[g]
}
ann_colors <- list(Sample_group = group_colors)

# row-wise z-score (heatmap “scale”)
mat_z <- t(scale(t(mat_c)))
mat_z[is.nan(mat_z)] <- 0

pdf(file.path(outdir, "Fig1_c.pdf"), width = 28, height = 10)
pheatmap(
  mat_z,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  fontsize = 8,
  border_color = NA,
  main = ""
)
invisible(dev.off())
message("Fig. 1c PDF: ", file.path(outdir, "Fig1_c.pdf"))

# --------------------------------------------------------------------------- #
# Fig. 1d — histograms healthy vs tumour
# Data: sheet "Figure 1d", columns: probe, SW-BCP-ALL-738 (tumour), NBC-01 (healthy)
# --------------------------------------------------------------------------- #
message("Reading Figure 1d …")
d <- as.data.frame(readxl::read_excel(xlsx, sheet = "Figure 1d"))
names(d)[1] <- "probe"
tumour_id <- names(d)[2]
healthy_id <- names(d)[3]
long_d <- data.frame(
  beta = c(d[[2]], d[[3]]),
  sample = rep(c(tumour_id, healthy_id), each = nrow(d))
)
long_d$sample <- factor(long_d$sample, levels = c(healthy_id, tumour_id))

p_d <- ggplot(long_d, aes(x = beta, fill = sample)) +
  geom_histogram(binwidth = 0.02, position = "identity", alpha = 0.5, boundary = 0) +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Fig. 1d — fCpG beta distributions (MOESM5)",
    subtitle = paste0("Healthy: ", healthy_id, " | Tumour: ", tumour_id),
    x = "Beta", y = "Count"
  ) +
  theme_bw()
ggsave(file.path(outdir, "Fig1_d.pdf"), p_d, width = 8, height = 5)
message("Fig. 1d PDF: ", file.path(outdir, "Fig1_d.pdf"))

# --------------------------------------------------------------------------- #
# Fig. 1e — seaborn vertical heatmap (RdBu_r, left colorbar, in-cell stars)
# Implemented in scripts/reproduce_fig1e.py (matches Nature panel layout).
# --------------------------------------------------------------------------- #
message("Fig. 1e via Python reproduce_fig1e.py …")
py_script <- normalizePath("d:/naturedna/evoflux-reproduce/scripts/reproduce_fig1e.py", winslash = "/", mustWork = TRUE)
py_exe <- Sys.which("python")
cp <- Sys.getenv("CONDA_PREFIX", "")
if (nzchar(cp)) {
  pe <- file.path(cp, "python.exe")
  if (file.exists(pe)) py_exe <- pe
}
if (!nzchar(py_exe) || !file.exists(py_exe)) {
  stop("Python not found; conda activate evoflux_nature2025 (or set PATH to python.exe)")
}
out_e <- normalizePath(file.path(outdir, "Fig1_e.pdf"), winslash = "/", mustWork = FALSE)
args <- c(py_script, "--xlsx", xlsx, "--out", out_e)
res <- system2(py_exe, args = args, stdout = TRUE, stderr = TRUE)
message(paste(res, collapse = "\n"))
if (!file.exists(out_e)) stop("Fig1_e.pdf was not created by reproduce_fig1e.py")

# --------------------------------------------------------------------------- #
# Fig. 1f — chromatin state log2 FC heatmap
# Data: sheet "Figure 1f"
# Column 1 = state label; Healthy cell, MCL, CLL, MM
# --------------------------------------------------------------------------- #
message("Reading Figure 1f …")
f <- as.data.frame(readxl::read_excel(xlsx, sheet = "Figure 1f"))
cn <- names(f)
state_col <- cn[1]
cell_cols <- cn[-1]
fl <- data.frame(state = f[[1]], f[, cell_cols, drop = FALSE], check.names = FALSE, stringsAsFactors = FALSE)
names(fl)[1] <- "state"
fl_long <- melt(fl, id.vars = "state", variable.name = "celltype", value.name = "log2fc")

p_f <- ggplot(fl_long, aes(x = celltype, y = state, fill = log2fc)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 1.5) +
  geom_text(aes(label = sprintf("%.2f", log2fc)), size = 2.8) +
  labs(title = "Fig. 1f — Chromatin state log2 FC, MOESM5", x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(outdir, "Fig1_f.pdf"), p_f, width = 7, height = 6)

# --------------------------------------------------------------------------- #
# Fig. 1g — match Duran-FerrerM-evoflux/code/Data_source_Fig.1G.Rmd (panels 2+3).
# Data: MOESM5 sheet "Figure 1g" (same structure as authors' gene_expression xlsx).
# Note: Rmd panel 1 (pp_plot) needs separate workbook; not in MOESM5 — omitted here.
# --------------------------------------------------------------------------- #
message("Reading Figure 1g …")
g_raw <- as.data.frame(readxl::read_excel(xlsx, sheet = "Figure 1g", col_names = FALSE, .name_repair = "minimal"))
g_left <- g_raw[3:nrow(g_raw), 1:4, drop = FALSE]
names(g_left) <- c("ensembl", "sample", "tpm", "fCpG")
g_left$tpm <- as.numeric(g_left$tpm)
g_left$fCpG <- as.logical(tolower(as.character(g_left$fCpG)) %in% c("true", "1", "yes", "t"))
g_left <- g_left[!is.na(g_left$tpm), ]

w <- wilcox.test(tpm ~ fCpG, data = g_left, alternative = "two.sided", exact = FALSE)
message("Fig1g left Wilcoxon W=", w$statistic, " p=", w$p.value)

random_sample_anonymous <- as.character(unique(g_left$sample)[1])
ymax <- max(g_left$tpm, na.rm = TRUE) + 1
# log10 axis cannot include 0; y is tpm+1 so lower limit 1 is valid (tpm >= 0)
ymin_log <- 1

# Left: same as Rmd — y = tpm+1, log10 axis, geom_pwc(wilcox_test), fill colours
plot_genes_sample <- g_left %>%
  ggplot(aes(x = fCpG, y = tpm + 1)) +
  geom_boxplot(aes(fill = fCpG), show.legend = FALSE, linewidth = 0.3) +
  geom_pwc(
    method = "wilcox_test",
    size = 0.3,
    p.adjust.method = "none",
    label = "P={p}",
    y.position = log10(8500),
    label.size = 1
  ) +
  scale_fill_manual(values = c("orange2", "grey70")) +
  coord_cartesian(ylim = c(ymin_log, ymax)) +
  scale_y_log10() +
  annotation_logticks(sides = "l", size = 0.3) +
  ggtitle(paste0("Anonymous=", random_sample_anonymous)) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "pt"),
    text = element_text(size = 6),
    line = element_line(linewidth = 0.3),
    panel.border = element_rect(fill = NA, linewidth = 0.3),
    panel.background = element_blank()
  )

# Right: beeswarm + semi-transparent boxplot + pairwise Wilcoxon (Rmd)
g_right <- g_raw[3:nrow(g_raw), 6:12, drop = FALSE]
names(g_right) <- c("cpg", "sample", "meth", "hgnc", "ensembl", "cpg2genes", "tpm")
g_right <- g_right[rowSums(is.na(g_right)) < ncol(g_right), ]
g_right$tpm <- suppressWarnings(as.numeric(as.character(g_right$tpm)))
g_right <- g_right[!is.na(g_right$tpm), ]
g_right$meth <- factor(as.character(as.integer(g_right$meth)), levels = c("0", "1", "2"))

plot_meth_alleles_sample <- g_right %>%
  ggplot(aes(x = meth, y = tpm + 1)) +
  geom_beeswarm(
    aes(fill = meth),
    color = "grey10",
    corral = "wrap",
    pch = 21,
    stroke = 0.1,
    size = 1
  ) +
  geom_boxplot(
    aes(fill = meth),
    alpha = 0.2,
    outlier.shape = NA,
    linewidth = 0.3
  ) +
  scale_fill_manual(values = c("0" = "#2F7FC2", "2" = "#BC5250", "1" = "grey80")) +
  geom_pwc(
    method = "wilcox_test",
    label = "P={p}",
    size = 0.3,
    label.size = 1,
    step.increase = 0.05
  ) +
  coord_cartesian(ylim = c(ymin_log, ymax)) +
  scale_y_log10() +
  annotation_logticks(sides = "l", size = 0.3) +
  ylab(NULL) +
  xlab("# of fCpG methylated alleles") +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "pt"),
    text = element_text(size = 6),
    line = element_line(linewidth = 0.3),
    panel.border = element_rect(fill = NA, linewidth = 0.3),
    panel.background = element_blank()
  )

# Rmd uses rel_widths c(1, 0.5, 0.75) for three panels; two-panel analogue:
p_g <- plot_grid(
  plot_genes_sample,
  plot_meth_alleles_sample,
  nrow = 1,
  align = "h",
  rel_widths = c(0.5, 0.75)
)

ggsave(file.path(outdir, "Fig1_g.pdf"), p_g, width = 7, height = 2.8)
message("Fig. 1g PDF (Rmd-style): ", file.path(outdir, "Fig1_g.pdf"))

message("All panels c–g written to ", outdir)
