# =============================================================================
# Fig. 5 — The evolutionary history of a tumour is prognostic of clinical outcome
#
# Reimplements logic from:
#   Duran-FerrerM-evoflux/docs/Data_source_Fig.5.html
# Data:
#   evoflux/41586_2025_9374_MOESM9_ESM.xlsx  (sheets Figure 5a + Figure 5c merged)
#
# Outputs (PDF, figures/reproduced/):
#   Fig5a_univariate_TTFT_OS.pdf
#   Fig5b_KM_TTFT_theta_IGHV.pdf
#   Fig5c_multivariate_Cox_TTFT.pdf
#   Fig5_Nature_layout.pdf   — composite (a/c left, b right) + caption
#
# Required R packages (install missing with install.packages / BiocManager if needed):
#   data.table, openxlsx, dplyr, tidyr, ggplot2, survival, survminer, ggsurvfit,
#   maxstat, ggpubr, ggpp, ggstats, cowplot, patchwork, glue, reshape2, janitor,
#   pals, scales
# =============================================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(survival)
  library(survminer)
  library(ggsurvfit)
  library(maxstat)
  library(ggpubr)
  library(ggpp)
  library(ggstats)
  library(glue)
  library(reshape2)
  library(cowplot)
  library(patchwork)
  library(pals)
  library(scales)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg)) {
  this_file <- sub("^--file=", "", file_arg[1])
  ROOT <- dirname(dirname(normalizePath(this_file, winslash = "/", mustWork = TRUE)))
} else {
  ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
XLSX <- file.path(ROOT, "evoflux/41586_2025_9374_MOESM9_ESM.xlsx")
OUT <- file.path(ROOT, "figures/reproduced")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

pdf_out <- function(path, width, height, plot_obj) {
  dev <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggsave(path, plot_obj, width = width, height = height, device = dev)
}

## ---- Load & merge sheets (MOESM9 has genomics on Figure 5c only) ----
fig5a <- openxlsx::read.xlsx(XLSX, sheet = "Figure 5a")
fig5c <- openxlsx::read.xlsx(XLSX, sheet = "Figure 5c")
join_cols <- c(
  "PARTICIPANT_ID_ANONYMOUS", "DISEASE_SUBTYPE", "AGE_SAMPLING",
  "Genomics.Mutation_TP53", "Genomics.loss_17p13.1", "Genomics.loss_17p"
)
fig5c <- fig5c[, intersect(join_cols, names(fig5c)), drop = FALSE]
cohort.discovery <- dplyr::inner_join(fig5a, fig5c, by = "PARTICIPANT_ID_ANONYMOUS")

cohort <- "Discovery"

## ---- Evolutionary variables: match Nature Fig. 5a (5 rows, not 8 mu/gamma/...) ----
evo.variables <- c(
  theta = "theta",
  Scancer = "Scancer",
  tau = "tau",
  cancerAge = "cancerAge",
  epiRate = "epiRate"
)
evo.variables.desciption <- structure(
  names(evo.variables),
  names = c(
    "Growth rate",
    "Effective population size",
    "Patient's age at MRCA",
    "Cancer age",
    "Mean epigenetic switching rate"
  )
)
y_labels_univar <- stats::setNames(names(evo.variables.desciption), unname(evo.variables.desciption))

## ---- (a) Univariate Cox: TTFT + OS (same formula pattern as Data_source_Fig.5.html) ----
TTFT.uni <- lapply(evo.variables, function(evo) {
  dat <- cohort.discovery %>%
    dplyr::filter(SAMPLING_TIME_TREATMENT_STATUS == "Untreated")
  if (evo == "Scancer") {
    dat[, evo] <- dat[, evo] / 1e6
  }
  if (evo == "tau") {
    dat[, evo] <- dat[, evo] / 10
  }
  if (evo == "epiRate") {
    dat[, evo] <- dat[, evo] * 100
  }
  cox <- coxph(
    Surv(time = as.numeric(Clinics.TTFT_DAYS_SAMPLING / 365.25), event = as.numeric(Clinics.TTFT)) ~ evo,
    data = dat %>%
      dplyr::select(Clinics.TTFT, Clinics.TTFT_DAYS_SAMPLING, evo = dplyr::all_of(evo)) %>%
      tidyr::drop_na()
  )
  ci <- summary(cox)$conf.int[1, c("lower .95", "upper .95")]
  data.frame(
    evo.variable = evo,
    end.point = "TTFT",
    Hazard.ratio = summary(cox)$coefficients[1, "exp(coef)"],
    lower..95 = ci[[1]],
    upper..95 = ci[[2]],
    pval = summary(cox)$coefficients[1, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
})
TTFT.uni <- do.call(rbind, TTFT.uni)

OS.uni <- lapply(evo.variables, function(evo) {
  dat <- cohort.discovery %>%
    dplyr::filter(SAMPLING_TIME_TREATMENT_STATUS == "Untreated")
  if (evo == "Scancer") {
    dat[, evo] <- dat[, evo] / 1e6
  }
  if (evo == "tau") {
    dat[, evo] <- dat[, evo] / 10
  }
  if (evo == "epiRate") {
    dat[, evo] <- dat[, evo] * 100
  }
  cox <- coxph(
    Surv(time = as.numeric(Clinics.OS_DAYS_SAMPLING / 365.25), event = as.numeric(Clinics.OS)) ~ evo,
    data = dat %>%
      dplyr::select(Clinics.OS, Clinics.OS_DAYS_SAMPLING, evo = dplyr::all_of(evo)) %>%
      tidyr::drop_na()
  )
  ci <- summary(cox)$conf.int[1, c("lower .95", "upper .95")]
  data.frame(
    evo.variable = evo,
    end.point = "OS",
    Hazard.ratio = summary(cox)$coefficients[1, "exp(coef)"],
    lower..95 = ci[[1]],
    upper..95 = ci[[2]],
    pval = summary(cox)$coefficients[1, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
})
OS.uni <- do.call(rbind, OS.uni)

cox.data <- rbind(TTFT.uni, OS.uni)

n_un <- cohort.discovery %>%
  dplyr::filter(SAMPLING_TIME_TREATMENT_STATUS == "Untreated") %>%
  nrow()
ev_ttft <- cohort.discovery %>%
  dplyr::filter(SAMPLING_TIME_TREATMENT_STATUS == "Untreated", Clinics.TTFT == 1) %>%
  nrow()
ev_os <- cohort.discovery %>%
  dplyr::filter(SAMPLING_TIME_TREATMENT_STATUS == "Untreated", Clinics.OS == 1) %>%
  nrow()

p_univar <- cox.data %>%
  dplyr::arrange(factor(evo.variable, levels = rev(names(evo.variables)))) %>%
  dplyr::mutate(
    evo.variable = factor(evo.variable, levels = unique(evo.variable)),
    pval = glue::glue("~italic(P) == {signif(pval, 2)}")
  ) %>%
  ggplot(aes(
    x = Hazard.ratio,
    y = evo.variable,
    xmin = lower..95,
    xmax = upper..95,
    color = end.point,
    shape = end.point,
    label = pval
  )) +
  geom_stripped_rows(col = NA, nudge_y = 0.1) +
  geom_vline(xintercept = 1, color = "grey70", lty = "dashed", linewidth = 0.2) +
  geom_errorbar(width = 0.25, position = position_dodge(width = 0.6), linewidth = 0.2) +
  geom_point(size = 1.5, position = position_dodge(width = 0.6)) +
  geom_text(
    alpha = 1,
    position = position_dodgenudge(width = 0.6, y = 0.22),
    size = 2,
    parse = TRUE,
    show.legend = FALSE,
    color = "grey0"
  ) +
  scale_shape_manual(
    "Univariate Cox analyses:",
    values = c("TTFT" = 15, "OS" = 19),
    limits = c("TTFT", "OS")
  ) +
  scale_color_manual(
    "P-value",
    values = structure(pals::coolwarm(2), names = c("TTFT", "OS"))
  ) +
  annotate(
    "text",
    label = c(
      paste0("TTFT, N=", n_un, ", events=", ev_ttft),
      paste0("OS, N=", n_un, ", events=", ev_os)
    ),
    x = 1.75,
    y = c(1.25, 1),
    hjust = 0,
    size = 1.75
  ) +
  guides(
    color = "none",
    shape = guide_legend(override.aes = list(color = structure(pals::coolwarm(2), names = c("TTFT", "OS"))))
  ) +
  scale_x_log10() +
  ylab(NULL) +
  xlab("Hazard ratio (95% CI)") +
  scale_y_discrete(labels = y_labels_univar) +
  theme_bw() +
  theme(
    text = element_text(color = "grey0", size = 6),
    line = element_line(linewidth = 0.2),
    axis.text = element_text(color = "grey0"),
    legend.key.size = unit(0, "pt"),
    legend.position = "top",
    legend.margin = margin(0, 0, 0, 0),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.box.spacing = unit(3, "pt"),
    legend.justification = "left",
    legend.title = element_text(vjust = 1),
    panel.border = element_rect(fill = NA, color = NA, linewidth = 0.3)
  )

pdf_out(file.path(OUT, "Fig5a_univariate_TTFT_OS.pdf"), 4.2, 4.5, p_univar)

## ---- (b) KM: maxstat split on theta within IGHV; TTFT outcome (HTML chunk 4) ----
variables <- c("theta", "Scancer")

df.KM <- lapply(split(seq_len(nrow(cohort.discovery)), cohort.discovery$DISEASE_SUBTYPE), function(indx) {
  dat <- cohort.discovery %>%
    dplyr::slice(indx) %>%
    dplyr::filter(SAMPLING_TIME_TREATMENT_STATUS == "Untreated") %>%
    dplyr::mutate(
      Clinics.TTFT = as.numeric(Clinics.TTFT),
      Clinics.OS = as.numeric(Clinics.OS)
    ) %>%
    dplyr::select(
      Clinics.TTFT, Clinics.TTFT_DAYS_SAMPLING,
      Clinics.OS, Clinics.OS_DAYS_SAMPLING,
      DISEASE_SUBTYPE,
      dplyr::all_of(variables)
    ) %>%
    tidyr::drop_na()

  cp <- surv_cutpoint(
    data = dat,
    time = "Clinics.TTFT_DAYS_SAMPLING",
    event = "Clinics.TTFT",
    variables = variables
  ) %>%
    surv_categorize() %>%
    data.frame()
  variables_cat <- cp %>% dplyr::select(dplyr::all_of(variables))
  colnames(variables_cat) <- paste0(colnames(variables_cat), ".IGHV.maxstat")
  variables_cat <- apply(variables_cat, 2, function(x) paste0(dat$DISEASE_SUBTYPE, "-", x))
  dat <- cbind(dat, variables_cat)
  dat
})
df.KM <- do.call(rbind, df.KM)

fit <- survfit2(
  Surv(time = Clinics.TTFT_DAYS_SAMPLING / 365.25, event = as.numeric(Clinics.TTFT)) ~ IGHV.theta,
  data = df.KM %>%
    dplyr::mutate(
      IGHV.theta = factor(
        theta.IGHV.maxstat,
        levels = c("mutated-low", "mutated-high", "unmutated-low", "unmutated-high"),
        labels = c(
          "M-CLL, low growth rate", "M-CLL, high growth rate",
          "U-CLL, low growth rate", "U-CLL, high growth rate"
        )
      )
    )
)

pairwise.pvals <- pairwise_survdiff(
  Surv(time = Clinics.TTFT_DAYS_SAMPLING / 365.25, event = as.numeric(Clinics.TTFT)) ~ IGHV.theta,
  data = df.KM %>%
    dplyr::mutate(
      IGHV.theta = factor(
        theta.IGHV.maxstat,
        levels = c("mutated-low", "mutated-high", "unmutated-low", "unmutated-high"),
        labels = c(
          "M-CLL, low growth rate", "M-CLL, high growth rate",
          "U-CLL, low growth rate", "U-CLL, high growth rate"
        )
      )
    ),
  p.adjust.method = "none"
)

p_km <- fit %>%
  ggsurvfit(type = "risk", linewidth = 0.3) +
  add_confidence_interval(alpha = 0.1) +
  add_risktable(
    risktable_stats = "n.risk",
    theme = theme_risktable_default(plot.title.size = 5),
    size = 2
  ) +
  add_risktable_strata_symbol() +
  add_censor_mark(size = 3, shape = "|") +
  coord_cartesian(xlim = c(0, 12)) +
  scale_x_continuous(breaks = 0:12) +
  ylab("Probability of treatment") +
  xlab("Time from sampling (years)") +
  scale_color_manual(
    values = c("#757099", "#2B2A65", "#E8A78F", "#C93A17"),
    labels = sapply(seq_along(fit$strata), function(i) {
      paste0(gsub("IGHV.theta=", "", names(fit$strata[i])), " (N=", fit$n[[i]], ")")
    }),
    aesthetics = c("color", "fill")
  ) +
  guides(color = guide_legend(ncol = 2)) +
  annotate(
    geom = "text",
    label = c(
      glue::glue("italic('P') == {signif(surv_pvalue(fit)$pval, 3)}"),
      glue::glue(
        "italic('P') == {signif(pairwise.pvals$p.value %>% reshape2::melt() %>% dplyr::filter(Var1=='U-CLL, high growth rate', Var2=='U-CLL, low growth rate') %>% dplyr::pull(value), 3)}"
      ),
      glue::glue(
        "italic('P') == {signif(pairwise.pvals$p.value %>% reshape2::melt() %>% dplyr::filter(Var1=='M-CLL, high growth rate', Var2=='M-CLL, low growth rate') %>% dplyr::pull(value), 3)}"
      )
    ),
    parse = TRUE,
    x = c(0.5, 10, 10),
    y = c(1, 0.95, 0.35),
    size = 2,
    hjust = 0,
    vjust = 0
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 5, colour = "grey0"),
    axis.text = element_text(size = 5, colour = "grey0"),
    axis.title = element_text(size = 6, colour = "grey0"),
    line = element_line(linewidth = 0.2),
    legend.key.size = unit(0, "pt"),
    legend.position = "top",
    legend.text = element_text(size = 5),
    legend.margin = margin(0, 0, 0, 0),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.spacing = unit(3, "pt"),
    legend.box.spacing = unit(3, "pt")
  )

pdf_out(file.path(OUT, "Fig5b_KM_TTFT_theta_IGHV.pdf"), 5.5, 5.2, p_km)

## ---- (c) Multivariate Cox TTFT (HTML chunk 6) ----
cox.data <- cohort.discovery %>%
  dplyr::filter(SAMPLING_TIME_TREATMENT_STATUS == "Untreated") %>%
  dplyr::mutate(
    Clinics.TTFT = as.numeric(Clinics.TTFT),
    Clinics.OS = as.numeric(Clinics.OS),
    Age = stats::I(AGE_SAMPLING / 10),
    TP53 = factor(
      ifelse(
        (Genomics.Mutation_TP53 != "WT") |
          (Genomics.loss_17p13.1 != "WT") |
          (!is.na(Genomics.loss_17p) & Genomics.loss_17p != "WT"),
        "Altered", "WT"
      ),
      levels = c("WT", "Altered")
    ),
    IGHV = factor(
      DISEASE_SUBTYPE,
      levels = c("mutated", "unmutated"),
      labels = c("M-CLL", "U-CLL")
    ),
    Scancer = Scancer / 1e6
  ) %>%
  dplyr::select(
    Clinics.TTFT,
    Clinics.TTFT_DAYS_SAMPLING,
    Clinics.OS,
    Clinics.OS_DAYS_SAMPLING,
    IGHV,
    TP53,
    theta,
    Scancer,
    Age
  ) %>%
  tidyr::drop_na()

cox.ttft.multi <- coxph(
  Surv(time = Clinics.TTFT_DAYS_SAMPLING / 365.25, event = Clinics.TTFT) ~ theta + IGHV + TP53 + Age,
  data = cox.data
)

df.cox <- cox.ttft.multi
coef_nm <- names(df.cox$coefficients)
coef_var <- vapply(coef_nm, function(nm) {
  if (nm == "theta") {
    return("theta")
  }
  if (grepl("^IGHV", nm)) {
    return("IGHV")
  }
  if (grepl("^TP53", nm)) {
    return("TP53")
  }
  if (nm == "Age") {
    return("Age")
  }
  nm
}, character(1), USE.NAMES = FALSE)

df.cox.multi <- data.frame(
  variable = coef_var,
  coefs = coef_nm,
  Hazard.ratio = summary(df.cox)$coefficients[, "exp(coef)"],
  summary(df.cox)$conf.int[, c("lower .95", "upper .95")],
  pval = summary(df.cox)$coefficients[, "Pr(>|z|)"],
  stringsAsFactors = FALSE
)
colnames(df.cox.multi)[4:5] <- c("lower..95", "upper..95")

df.cox.multi$variable <- sapply(seq_len(nrow(df.cox.multi)), function(i) {
  if (df.cox.multi$variable[i] == df.cox.multi$coefs[i]) {
    df.cox.multi$variable[i]
  } else {
    var <- df.cox.multi$variable[i]
    var.stata <- gsub(var, "", df.cox.multi$coefs[i], fixed = TRUE)
    paste0(var, "\n(N=", length(which(cox.data[[var]] == var.stata)), ")")
  }
})

p_multi <- df.cox.multi %>%
  dplyr::arrange(Hazard.ratio) %>%
  dplyr::mutate(
    variable = factor(variable, levels = unique(variable)),
    pval = ifelse(!is.na(pval), glue::glue("~italic(P) == {signif(pval, 2)}"), NA_character_)
  ) %>%
  ggplot(aes(
    x = Hazard.ratio,
    y = variable,
    xmin = lower..95,
    xmax = upper..95,
    label = pval
  )) +
  geom_stripped_rows(col = NA, nudge_y = 0.1) +
  geom_vline(xintercept = 1, color = "grey70", lty = "dashed", linewidth = 0.2) +
  geom_errorbar(width = 0.25, position = position_dodge(width = 0.6), linewidth = 0.2) +
  geom_point(size = 1.5, shape = 15, position = position_dodge(width = 0.6)) +
  geom_text(
    alpha = 1,
    position = position_dodgenudge(width = 0.6, y = 0.4),
    size = 2,
    parse = TRUE,
    show.legend = FALSE,
    color = "grey0"
  ) +
  annotate(
    geom = "text",
    label = c(
      paste0("N=", df.cox$n),
      paste0("Events=", df.cox$nevent, ", CI=", round(summary(df.cox)$concordance["C"], 3)),
      paste0("Global p-value=", signif(summary(df.cox)$logtest[["pvalue"]], 3))
    ),
    x = 2.75,
    y = c(1.3, 1, 0.7),
    hjust = 0,
    size = 1.5
  ) +
  scale_x_log10() +
  ylab(NULL) +
  xlab("Hazard ratio (95% CI)") +
  theme_bw() +
  theme(
    text = element_text(color = "grey0", size = 6),
    line = element_line(linewidth = 0.2),
    axis.text = element_text(color = "grey0"),
    legend.key.size = unit(0, "pt"),
    legend.position = "top",
    legend.margin = margin(0, 0, 0, 0),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.box.spacing = unit(3, "pt"),
    legend.justification = "left",
    panel.border = element_rect(fill = NA, color = NA, linewidth = 0.3)
  )

pdf_out(file.path(OUT, "Fig5c_multivariate_Cox_TTFT.pdf"), 4.5, 3.2, p_multi)

## ---- Composite layout (paper: b right full height; a & c stacked left) ----
caption <- paste0(
  "Fig. 5 | The evolutionary history of a tumour is prognostic of clinical outcome. ",
  "a, Univariate survival analysis of the TTFT (blue) and overall survival (OS; red) ",
  "in the discovery CLL cohort for evolutionary variables inferred via EVOFLUx. ",
  "b, Kaplan–Meier curves comparing the TTFT between patients with high versus low ",
  "inferred cancer growth rates, separated by IGHV mutational status. ",
  "c, Multivariate Cox regression model of the TTFT shows that the cancer growth rate ",
  "is significant when controlling for IGHV status, TP53 alterations and age at sampling. ",
  "The error bars represent 95% confidence intervals. ",
  "A log-rank test was used in the Kaplan–Meier curves and Wald tests for Cox models. ",
  "A Schoenfeld residuals test was used to test proportional hazard assumptions. ",
  "No multiple comparison adjustments were done."
)

layout_fig <- (p_univar / p_multi) | p_km +
  patchwork::plot_layout(widths = c(1, 1.05), guides = "collect") +
  patchwork::plot_annotation(
    caption = caption,
    theme = theme(
      plot.caption = element_text(size = 5.5, hjust = 0, lineheight = 1.15),
      plot.caption.position = "plot"
    )
  )

pdf_out(file.path(OUT, "Fig5_Nature_layout.pdf"), 10.5, 7.2, layout_fig)

message("Saved PDFs to: ", OUT)
