suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(tibble)
})

options(stringsAsFactors = FALSE, "openxlsx.dateFormat" = "dd/mm/yyyy")

MOESM3 <- "d:/naturedna/evoflux-reproduce/evoflux/41586_2025_9374_MOESM3_ESM.xlsx"
OUT_DIR <- "d:/naturedna/evoflux-reproduce/figures/Extended_Data_Fig"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# 1) Load data from MOESM3
# ----------------------------
# Equivalent to fCpGs.betas in the HTML: fCpG beta values in 2,204 samples
fCpGs.betas <- openxlsx::read.xlsx(MOESM3, sheet = "supplementary_table_5", startRow = 3)
stopifnot(ncol(fCpGs.betas) >= 2)
fCpGs.betas <- fCpGs.betas %>% column_to_rownames(var = "fCpGs")

# Metadata from the cohort
meta <- openxlsx::read.xlsx(MOESM3, sheet = "supplementary_table_2", startRow = 3)
qc <- openxlsx::read.xlsx(MOESM3, sheet = "supplementary_table_3", startRow = 3)
clocks <- openxlsx::read.xlsx(MOESM3, sheet = "supplementary_table_9", startRow = 3)

meta <- meta %>% rename(PARTICIPANT_ID_ANONYMOUS = participant_id_anonymous)
qc <- qc %>% rename(PARTICIPANT_ID_ANONYMOUS = participant_id_anonymous)
clocks <- clocks %>% rename(PARTICIPANT_ID_ANONYMOUS = PARTICIPANT_ID_ANONYMOUS)

# Merge: a single per-sample metadata table similar to fCpGs.metadata.all in HTML
df <- meta %>%
  left_join(clocks %>% select(PARTICIPANT_ID_ANONYMOUS, Clocks_Horvath), by = "PARTICIPANT_ID_ANONYMOUS") %>%
  mutate(
    CELL_TYPE = cell_type,
    CELL_TYPE_ANNOT.1 = cell_type_annot_1,
    CELL_TYPE_ANNOT.2 = cell_type_annot_2,
    CELL_TYPE_ANNOT.3 = cell_type_annot_3,
    AGE_SAMPLING = age_sampling,
    Horvath = Clocks_Horvath
  )

# Keep only samples that exist in fCpGs.betas
present <- intersect(colnames(fCpGs.betas), df$PARTICIPANT_ID_ANONYMOUS)
df <- df %>% filter(PARTICIPANT_ID_ANONYMOUS %in% present)
fCpGs.betas <- fCpGs.betas[, df$PARTICIPANT_ID_ANONYMOUS, drop = FALSE]

# ----------------------------
# 2) Compute mean/SD of fCpGs
# ----------------------------
df$SD.fCpGs <- apply(fCpGs.betas, 2, sd, na.rm = TRUE)
df$mean.fCpGs <- apply(fCpGs.betas, 2, mean, na.rm = TRUE)

# ----------------------------
# 3) Subset similar to HTML normal cell-type plots
# ----------------------------
# HTML logic:
# - keep grepl("^Normal", CELL_TYPE_ANNOT.1)
# - exclude some in-vitro/precursor categories
df.plot <- df %>%
  filter(
    !is.na(CELL_TYPE_ANNOT.1),
    grepl("^Normal", CELL_TYPE_ANNOT.1),
    !is.na(Horvath),
    !is.na(AGE_SAMPLING)
  ) %>%
  mutate(
    CELL_TYPE_ANNOT.2 = case_when(
      CELL_TYPE == "CD19_positive" ~ "Bcell (CD19+)",
      CELL_TYPE == "Granulocyte" ~ "Granulocyte",
      CELL_TYPE == "NK_cell" ~ "NK",
      CELL_TYPE_ANNOT.2 == "Tcell" ~ "Tcell(CD3+)",
      CELL_TYPE == "Monocyte" ~ "Monocyte",
      TRUE ~ CELL_TYPE_ANNOT.2
    )
  )

# ----------------------------
# 4) Reproduce key plots (from HTML sections around lines ~2086 & ~2140)
# ----------------------------
base_theme <- theme_classic() +
  theme(
    text = element_text(color = "grey0", size = 6),
    line = element_line(linewidth = 0.2),
    axis.text = element_text(color = "grey0"),
    legend.position = "right",
    panel.border = element_rect(fill = NA, color = NA, linewidth = 0.3)
  )

p_age_vs_horvath <-
  ggplot(df.plot, aes(x = Horvath, y = AGE_SAMPLING)) +
  geom_point(size = 0.6, alpha = 0.8) +
  geom_smooth(method = "lm", color = "grey0", linewidth = 0.3, se = FALSE) +
  stat_cor(size = 1.8) +
  xlab("Horvath") +
  ylab("Age at sampling") +
  ggtitle("Normal lymphoid/myeloid cells") +
  base_theme

p_horvath_hist <-
  ggplot(df.plot, aes(x = Horvath)) +
  geom_histogram(bins = 30, color = "grey10", fill = "grey80", linewidth = 0.2) +
  ggtitle(paste0("Normal lymphoid and myeloid cells, N=", nrow(df.plot))) +
  xlab("Horvath") +
  ylab("Count") +
  base_theme +
  theme(legend.position = "none")

p_mean_vs_horvath <-
  ggplot(df.plot, aes(x = Horvath, y = mean.fCpGs)) +
  geom_point(aes(fill = CELL_TYPE_ANNOT.2), size = 0.75, pch = 21, stroke = 0, alpha = 0.85) +
  xlab("Horvath") +
  ylab("mean fCpGs") +
  ggtitle("Mean fCpGs vs Horvath") +
  base_theme

p_sd_vs_horvath <-
  ggplot(df.plot, aes(x = Horvath, y = SD.fCpGs)) +
  geom_point(aes(fill = CELL_TYPE_ANNOT.2), size = 0.75, pch = 21, stroke = 0, alpha = 0.85) +
  xlab("Horvath") +
  ylab("SD fCpGs") +
  ggtitle("SD fCpGs vs Horvath") +
  base_theme

# Layouts: mimic the HTML's patchwork usage with "|" (side-by-side)
fig1 <- p_age_vs_horvath | p_horvath_hist
fig2 <- p_mean_vs_horvath | p_sd_vs_horvath

ggsave(
  filename = file.path(OUT_DIR, "fCpGs_Aging_moesm3_Horvath_vs_Age_and_hist.pdf"),
  plot = fig1,
  width = 7.2,
  height = 3.0,
  units = "in"
)

ggsave(
  filename = file.path(OUT_DIR, "fCpGs_Aging_moesm3_Horvath_vs_mean_and_SD.pdf"),
  plot = fig2,
  width = 7.2,
  height = 3.0,
  units = "in"
)

writeLines("WROTE:")
writeLines(file.path(OUT_DIR, "fCpGs_Aging_moesm3_Horvath_vs_Age_and_hist.pdf"))
writeLines(file.path(OUT_DIR, "fCpGs_Aging_moesm3_Horvath_vs_mean_and_SD.pdf"))

