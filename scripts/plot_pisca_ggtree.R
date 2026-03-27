#!/usr/bin/env Rscript
# Plot PISCA / TreeAnnotator MCC tree with ggtree (Fig.4-style).
# Usage:
#   Rscript scripts/plot_pisca_ggtree.R figures/reproduced/pisca_zenodo_fig4/fig4a_case12_zenodo_MCC.tree \
#     figures/reproduced/pisca_zenodo_fig4/fig4a_case12_zenodo_ggtree.pdf
#
# Requires: ape, ggtree (Bioconductor), ggplot2

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: plot_pisca_ggtree.R <MCC_newick_or_nexus.tree> <out.pdf>")
}

tree_file <- args[[1]]
out_pdf <- args[[2]]

if (!requireNamespace("ape", quietly = TRUE)) stop("Install R package: ape")
if (!requireNamespace("ggtree", quietly = TRUE)) {
  stop("Install Bioconductor package ggtree (BiocManager::install('ggtree'))")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2")

tr <- tryCatch(
  ape::read.nexus(tree_file),
  error = function(e1) {
    tryCatch(ape::read.tree(tree_file), error = function(e2) stop(e1))
  }
)

p <- ggtree::ggtree(tr, ladderize = TRUE) +
  ggtree::geom_tiplab(size = 3) +
  ggtree::theme_tree2()

ggplot2::ggsave(out_pdf, p, width = 8, height = 6)
message("Wrote ", out_pdf)
