library(readxl)
x <- "d:/naturedna/evoflux-reproduce/evoflux/41586_2025_9374_MOESM5_ESM.xlsx"
g <- read_excel(x, "Figure 1g", col_names = FALSE)
r <- g[3:nrow(g), 6:12]
names(r) <- c("cpg", "sample", "meth", "hgnc", "ensembl", "cpg2genes", "tpm")
r$tpm <- as.numeric(r$tpm)
r$meth <- factor(as.integer(r$meth), levels = c(0, 1, 2))
print(kruskal.test(tpm ~ meth, data = r))
print(pairwise.wilcox.test(r$tpm, r$meth, p.adjust.method = "none", exact = FALSE))
