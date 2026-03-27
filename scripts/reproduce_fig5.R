# =============================================================================
# Figure 5 — Nature 2025 (EVOFLUx / CLL discovery cohort)
# Based on: Duran-FerrerM-evoflux/docs/Data_source_Fig.5.html (original R logic)
# Data:     MOESM9 41586_2025_9374_MOESM9_ESM.xlsx (sheets Figure 5a + Figure 5c)
#
# Outputs (PDF, figures/reproduced/):
#   Fig5a.pdf — Panel (a): univariate Cox, TTFT (blue) & OS (red), 5 EVOFLUx variables
#               (aligned with main figure: growth rate, N_e, MRCA age, cancer age, mean epi rate)
#   Fig5b.pdf — Panel (b): KM cumulative incidence of TTFT by IGHV × maxstat split on theta
#   Fig5c.pdf — Panel (c): multivariate Cox for TTFT (theta + IGHV + TP53 + age)
#
# Caption reference:
# "Fig. 5 | The evolutionary history of a tumour is prognostic of clinical outcome..."
# =============================================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(survival)
  library(survminer)
  library(ggsurvfit)
  library(ggplot2)
  library(ggpp)
  library(patchwork)
  library(glue)
  library(pals)
  library(reshape2)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg)) {
  ROOT <- dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  ROOT <- normalizePath(getwd())
}
# scripts/reproduce_fig5.R -> project root evoflux-reproduce
ROOT <- normalizePath(file.path(ROOT, ".."), winslash = "/", mustWork = FALSE)

XLSX <- file.path(ROOT, "evoflux", "41586_2025_9374_MOESM9_ESM.xlsx")
OUT_DIR <- file.path(ROOT, "figures", "reproduced")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Load & merge sheets (5a: EVOFLUx + clinics; 5c: subtype, age, genomics) ----------
fig5a <- read.xlsx(XLSX, sheet = "Figure 5a")
fig5c <- read.xlsx(XLSX, sheet = "Figure 5c")

join_cols <- c(
  "PARTICIPANT_ID_ANONYMOUS",
  "DISEASE_SUBTYPE",
  "AGE_SAMPLING",
  "Genomics.Mutation_TP53",
  "Genomics.loss_17p13.1",
  "Genomics.loss_17p"
)
fig5c_sub <- fig5c[, intersect(join_cols, names(fig5c))]

cohort.discovery <- fig5a %>%
  left_join(fig5c_sub, by = "PARTICIPANT_ID_ANONYMOUS", suffix = c("", ".from5c"))

if (!("DISEASE_SUBTYPE" %in% names(cohort.discovery)) &&
    ("DISEASE_SUBTYPE.from5c" %in% names(cohort.discovery))) {
  cohort.discovery$DISEASE_SUBTYPE <- cohort.discovery$DISEASE_SUBTYPE.from5c
}

if (("DISEASE_SUBTYPE" %in% names(cohort.discovery)) &&
    ("DISEASE_SUBTYPE.from5c" %in% names(cohort.discovery))) {
  cohort.discovery$DISEASE_SUBTYPE <- ifelse(
    is.na(cohort.discovery$DISEASE_SUBTYPE),
    cohort.discovery$DISEASE_SUBTYPE.from5c,
    cohort.discovery$DISEASE_SUBTYPE
  )
}

if ("DISEASE_SUBTYPE.from5c" %in% names(cohort.discovery)) {
  cohort.discovery <- cohort.discovery %>% select(-DISEASE_SUBTYPE.from5c)
}

stopifnot(nrow(cohort.discovery) > 0)

# Paper panel (a): five evolutionary quantities (matches figure labelling)
evo.variables <- c(
  theta = "theta",
  Scancer = "Scancer",
  tau = "tau",
  cancerAge = "cancerAge",
  epiRate = "epiRate"
)

evo.labels <- c(
  theta = "Growth rate",
  Scancer = "Effective population size",
  tau = "Patient's age at MRCA",
  cancerAge = "Cancer age",
  epiRate = "Mean epigenetic switching rate"
)

# --- Helpers: scale predictors like HTML notebook -----------------------------------
scale_evo_col <- function(dat, evo) {
  if (evo %in% c("Scancer")) {
    dat[[evo]] <- dat[[evo]] / 1e6
  }
  if (evo == "tau") {
    dat[[evo]] <- dat[[evo]] / 10
  }
  if (evo %in% c("mu", "gamma", "nu")) {
    dat[[evo]] <- dat[[evo]] * 100
  }
  if (evo == "zeta") {
    dat[[evo]] <- dat[[evo]] * 1000
  }
  if (evo == "epiRate") {
    dat[[evo]] <- dat[[evo]] * 100
  }
  dat
}

# --- (a) Univariate Cox TTFT + OS ---------------------------------------------------
run_uni_cox <- function(endpoint) {
  rows <- lapply(names(evo.variables), function(evo) {
    evo <- evo.variables[[evo]]
    dat <- cohort.discovery %>%
      filter(SAMPLING_TIME_TREATMENT_STATUS == "Untreated") %>%
      scale_evo_col(evo)

    if (endpoint == "TTFT") {
      fml <- as.formula(paste0(
        "Surv(time = as.numeric(Clinics.TTFT_DAYS_SAMPLING/365.25), ",
        "event = as.numeric(Clinics.TTFT)) ~ ", evo
      ))
      dat <- dat %>%
        select(Clinics.TTFT, Clinics.TTFT_DAYS_SAMPLING, all_of(evo)) %>%
        drop_na()
    } else {
      fml <- as.formula(paste0(
        "Surv(time = as.numeric(Clinics.OS_DAYS_SAMPLING/365.25), ",
        "event = as.numeric(Clinics.OS)) ~ ", evo
      ))
      dat <- dat %>%
        select(Clinics.OS, Clinics.OS_DAYS_SAMPLING, all_of(evo)) %>%
        drop_na()
    }

    cox <- coxph(fml, data = dat)
    s <- summary(cox)
    ci <- s$conf.int
    data.frame(
      evo.variable = evo,
      end.point = endpoint,
      Hazard.ratio = ci[, "exp(coef)"],
      lower..95 = ci[, "lower .95"],
      upper..95 = ci[, "upper .95"],
      pval = s$coefficients[, "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

TTFT.uni <- run_uni_cox("TTFT")
OS.uni <- run_uni_cox("OS")
cox.data <- rbind(TTFT.uni, OS.uni)

p_a <- cox.data %>%
  mutate(
    evo.variable = factor(evo.variable, levels = rev(names(evo.variables))),
    pval = glue("~italic(P) == {signif(pval, 2)}")
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
  NULL +
  geom_vline(xintercept = 1, color = "grey70", lty = "dashed", lwd = 0.2) +
  geom_errorbar(width = 0.25, position = position_dodge(width = 0.6), lwd = 0.2) +
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
      paste0(
        "TTFT, N=",
        nrow(filter(cohort.discovery, SAMPLING_TIME_TREATMENT_STATUS == "Untreated")),
        ", events=",
        nrow(filter(cohort.discovery, SAMPLING_TIME_TREATMENT_STATUS == "Untreated", Clinics.TTFT == 1))
      ),
      paste0(
        "OS, N=",
        nrow(filter(cohort.discovery, SAMPLING_TIME_TREATMENT_STATUS == "Untreated")),
        ", events=",
        nrow(filter(cohort.discovery, SAMPLING_TIME_TREATMENT_STATUS == "Untreated", Clinics.OS == 1))
      )
    ),
    x = 1.75,
    y = c(1.25, 1),
    hjust = 0,
    size = 1.75
  ) +
  guides(
    color = "none",
    alpha = "none",
    shape = guide_legend(override.aes = list(color = structure(pals::coolwarm(2), names = c("TTFT", "OS"))))
  ) +
  scale_x_log10() +
  ylab(NULL) +
  xlab("Hazard ratio (95% CI)") +
  scale_y_discrete(labels = function(x) unname(evo.labels[x])) +
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

ggsave(file.path(OUT_DIR, "Fig5a.pdf"), p_a, width = 5.2, height = 4.2)

# --- (b) KM: maxstat on theta within IGHV, TTFT outcome -----------------------------
variables <- c("theta", "Scancer")

cohort.km <- cohort.discovery %>% filter(!is.na(DISEASE_SUBTYPE))

df.KM <- lapply(split(seq_len(nrow(cohort.km)), cohort.km$DISEASE_SUBTYPE), function(indx) {
  dat <- cohort.km %>%
    slice(indx) %>%
    filter(SAMPLING_TIME_TREATMENT_STATUS == "Untreated") %>%
    mutate(
      Clinics.TTFT = as.numeric(Clinics.TTFT),
      Clinics.OS = as.numeric(Clinics.OS)
    ) %>%
    select(
      Clinics.TTFT, Clinics.TTFT_DAYS_SAMPLING,
      Clinics.OS, Clinics.OS_DAYS_SAMPLING,
      DISEASE_SUBTYPE,
      all_of(variables)
    ) %>%
    drop_na()

  if (nrow(dat) < 10) {
    return(NULL)
  }

  cp <- tryCatch(
    surv_cutpoint(
      data = dat,
      time = "Clinics.TTFT_DAYS_SAMPLING",
      event = "Clinics.TTFT",
      variables = variables,
      minprop = 0.1
    ),
    error = function(e) NULL
  )
  if (is.null(cp)) {
    return(NULL)
  }

  sc <- surv_categorize(cp) %>% data.frame()
  vars_df <- sc[, variables, drop = FALSE]
  colnames(vars_df) <- paste0(colnames(vars_df), ".IGHV.maxstat")
  vars_df <- apply(vars_df, 2, function(x) paste0(dat$DISEASE_SUBTYPE, "-", x))
  cbind(dat, vars_df)
})

df.KM <- do.call(rbind, df.KM[!vapply(df.KM, is.null, logical(1))])
rownames(df.KM) <- NULL

ighv_theta_factor <- function(dat) {
  dat %>%
    mutate(
      IGHV.theta = factor(
        theta.IGHV.maxstat,
        levels = c("mutated-low", "mutated-high", "unmutated-low", "unmutated-high"),
        labels = c(
          "M-CLL, low growth rate", "M-CLL, high growth rate",
          "U-CLL, low growth rate", "U-CLL, high growth rate"
        )
      )
    )
}

df_km_f <- ighv_theta_factor(df.KM)

fit <- survfit2(
  Surv(time = Clinics.TTFT_DAYS_SAMPLING / 365.25, event = as.numeric(Clinics.TTFT)) ~ IGHV.theta,
  data = df_km_f
)

pairwise.pvals <- pairwise_survdiff(
  Surv(time = Clinics.TTFT_DAYS_SAMPLING / 365.25, event = as.numeric(Clinics.TTFT)) ~ IGHV.theta,
  data = df_km_f,
  p.adjust.method = "none"
)

p_u <- pairwise.pvals$p.value
pv_u_high_low <- tryCatch(
  reshape2::melt(p_u) %>%
    filter(as.character(Var1) == "U-CLL, high growth rate", as.character(Var2) == "U-CLL, low growth rate") %>%
    dplyr::slice(1) %>%
    pull(value),
  error = function(e) NA_real_
)
pv_m_high_low <- tryCatch(
  reshape2::melt(p_u) %>%
    filter(as.character(Var1) == "M-CLL, high growth rate", as.character(Var2) == "M-CLL, low growth rate") %>%
    dplyr::slice(1) %>%
    pull(value),
  error = function(e) NA_real_
)

pv_global <- tryCatch(
  survminer::surv_pvalue(fit)$pval,
  error = function(e) {
    sd <- survival::survdiff(
      Surv(time = Clinics.TTFT_DAYS_SAMPLING / 365.25, event = as.numeric(Clinics.TTFT)) ~ IGHV.theta,
      data = df_km_f
    )
    stats::pchisq(sd$chisq, df = length(sd$n) - 1, lower.tail = FALSE)
  }
)

p_b <- fit %>%
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
    "text",
    label = c(
      glue("italic('P') == {signif(pv_global, 3)}"),
      glue("italic('P') == {signif(pv_u_high_low, 3)}"),
      glue("italic('P') == {signif(pv_m_high_low, 3)}")
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

ggsave(file.path(OUT_DIR, "Fig5b.pdf"), p_b, width = 6.5, height = 5.5)

# --- (c) Multivariate Cox TTFT (HTML + paper) --------------------------------------
cox.data.multi <- cohort.discovery %>%
  filter(SAMPLING_TIME_TREATMENT_STATUS == "Untreated") %>%
  mutate(
    Clinics.TTFT = as.numeric(Clinics.TTFT),
    Clinics.OS = as.numeric(Clinics.OS),
    Age = AGE_SAMPLING / 10,
    TP53 = factor(
      ifelse(
        (!is.na(Genomics.Mutation_TP53) & Genomics.Mutation_TP53 != "WT") |
          (!is.na(Genomics.loss_17p13.1) & Genomics.loss_17p13.1 != "WT") |
          (!is.na(Genomics.loss_17p) & Genomics.loss_17p != "WT"),
        "Altered", "WT"
      ),
      levels = c("WT", "Altered")
    ),
    IGHV = factor(DISEASE_SUBTYPE, levels = c("mutated", "unmutated"), labels = c("M-CLL", "U-CLL")),
    Scancer = Scancer / 1e6
  ) %>%
  filter(!is.na(DISEASE_SUBTYPE)) %>%
  select(
    Clinics.TTFT, Clinics.TTFT_DAYS_SAMPLING,
    Clinics.OS, Clinics.OS_DAYS_SAMPLING,
    IGHV, TP53, theta, Scancer, Age
  ) %>%
  drop_na()

cox.ttft.multi <- coxph(
  Surv(time = Clinics.TTFT_DAYS_SAMPLING / 365.25, event = Clinics.TTFT) ~ theta + IGHV + TP53 + Age,
  data = cox.data.multi
)

s_multi <- summary(cox.ttft.multi)
df.cox.multi <- data.frame(
  variable = names(coef(cox.ttft.multi)),
  coefs = names(coef(cox.ttft.multi)),
  Hazard.ratio = s_multi$conf.int[, "exp(coef)"],
  lower..95 = s_multi$conf.int[, "lower .95"],
  upper..95 = s_multi$conf.int[, "upper .95"],
  pval = s_multi$coefficients[, "Pr(>|z|)"],
  stringsAsFactors = FALSE
)

df.cox.multi$ylab <- vapply(seq_len(nrow(df.cox.multi)), function(i) {
  v <- df.cox.multi$variable[i]
  if (grepl("^IGHV", v)) {
    return("IGHV U-CLL")
  }
  if (grepl("^TP53", v)) {
    return("TP53")
  }
  if (v == "theta") {
    return("Growth rate")
  }
  if (v == "Age") {
    return("Age")
  }
  v
}, character(1))

df.cox.multi <- df.cox.multi %>%
  mutate(
    ylab = factor(ylab, levels = c("Age", "TP53", "Growth rate", "IGHV U-CLL")),
    pval = ifelse(!is.na(pval), glue("~italic(P) == {signif(pval, 2)}"), NA)
  ) %>%
  arrange(ylab)

p_c <- df.cox.multi %>%
  ggplot(aes(
    x = Hazard.ratio,
    y = ylab,
    xmin = lower..95,
    xmax = upper..95,
    label = pval
  )) +
  NULL +
  geom_vline(xintercept = 1, color = "grey70", lty = "dashed", lwd = 0.2) +
  geom_errorbar(width = 0.25, position = position_dodge(width = 0.6), lwd = 0.2) +
  geom_point(size = 1.5, shape = 15, position = position_dodge(width = 0.6), color = "#3B5998") +
  geom_text(
    alpha = 1,
    position = position_dodgenudge(width = 0.6, y = 0.4),
    size = 2,
    parse = TRUE,
    show.legend = FALSE,
    color = "grey0"
  ) +
  annotate(
    "text",
    label = c(
      paste0("N=", cox.ttft.multi$n),
      paste0("Events=", cox.ttft.multi$nevent, ", C-index=", round(unname(s_multi$concordance["C"]), 3)),
      paste0("Global p-value=", signif(s_multi$logtest[["pvalue"]], 3))
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
    legend.position = "none",
    panel.border = element_rect(fill = NA, color = NA, linewidth = 0.3)
  )

ggsave(file.path(OUT_DIR, "Fig5c.pdf"), p_c, width = 5.2, height = 3.8)

## --- Combined Nature-style layout (a/c left, b right) --------------------------------
caption_txt <- paste0(
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

layout_fig <- ((p_a + ggtitle("a")) / (p_c + ggtitle("c"))) | (p_b + ggtitle("b"))
layout_fig <- layout_fig +
  plot_layout(widths = c(1, 1.35)) +
  plot_annotation(
    caption = caption_txt,
    theme = theme(
      plot.caption = element_text(size = 5.5, hjust = 0, lineheight = 1.15),
      plot.caption.position = "plot"
    )
  )

ggsave(file.path(OUT_DIR, "Fig5.pdf"), layout_fig, width = 10.6, height = 7.2)

message("Wrote: ", file.path(OUT_DIR, "Fig5a.pdf"))
message("Wrote: ", file.path(OUT_DIR, "Fig5b.pdf"))
message("Wrote: ", file.path(OUT_DIR, "Fig5c.pdf"))
message("Wrote: ", file.path(OUT_DIR, "Fig5.pdf"))

message("\n--- Schoenfeld residuals (TTFT multivariate Cox; proportional hazards) ---")
print(summary(survival::cox.zph(cox.ttft.multi)))
message("\n--- sessionInfo() ---")
print(sessionInfo())
