# Supplementary binary calibration.
# Reproduces:
#   SuppFig/sex_manhattan.pdf
#   SuppFig/compare_binary_logp_sex.png
#   SuppFig/binary_trait_qq_plot.png

library(data.table)
library(ggplot2)
library(cowplot)

sex_plink_file <- "data/binary/sex_plink.table"
sex_ukc_file <- "data/binary/sex_ukc.table"
snoring_files <- list(plink = "data/binary/X1210_plink.table", ukc = "data/binary/X1210_ukc.table")
broken_files <- list(plink = "data/binary/X2040_plink.table", ukc = "data/binary/X2040_ukc.table")
risk_files <- list(plink = "data/binary/X2463_plink.table", ukc = "data/binary/X2463_ukc.table")
supp_dir <- "SuppFig"

load_binary <- function(path) {
  dt <- fread(path)
  if ("TEST" %in% names(dt)) dt <- dt[TEST %in% c("ADD", "ADD_BOLT", "ADDitive", "ADD_ALT")]
  if (!("-log10_P" %in% names(dt)) && "P" %in% names(dt)) {
    dt[, `-log10_P` := -log10(P)]
  }
  dt
}

manhattan_basic <- function(dt, title) {
  if (!all(c("#CHROM", "POS", "-log10_P") %in% names(dt))) stop("Need #CHROM, POS, -log10_P for Manhattan plot.")
  ggplot(dt, aes(x = POS, y = `-log10_P`, colour = factor(`#CHROM` %% 2))) +
    geom_point(size = 0.4, alpha = 0.6) +
    labs(x = "Position", y = expression(-log[10](p)), title = title) +
    theme_bw() +
    theme(legend.position = "none", panel.grid = element_blank())
}

compare_logp <- function(plink_dt, ukc_dt, title) {
  id_col <- if ("ID" %in% names(plink_dt) && "ID" %in% names(ukc_dt)) "ID" else "SNP"
  merged <- merge(plink_dt, ukc_dt, by = id_col, suffixes = c("_plink", "_ukc"))
  ggplot(merged, aes(x = `-log10_P_plink`, y = `-log10_P_ukc`)) +
    geom_point(alpha = 0.4, size = 0.6, colour = "#1663A9") +
    geom_abline(intercept = 0, slope = 1, colour = "grey60") +
    labs(x = "PLINK -log10(p)", y = "UKC -log10(p)", title = title) +
    theme_bw() +
    theme(panel.grid = element_blank())
}

compute_lambda_gc <- function(p) {
  p <- p[!is.na(p) & p > 0 & p < 1]
  chisq <- qchisq(1 - p, df = 1)
  median(chisq) / qchisq(0.5, df = 1)
}

make_qq_data <- function(p) {
  p <- p[!is.na(p) & p > 0 & p < 1]
  n <- length(p)
  p <- sort(p)
  data.frame(
    exp = -log10((1:n - 0.5) / n),
    obs = -log10(p)
  )
}

qq_panel <- function(p, label, colour) {
  d <- make_qq_data(p)
  lambda <- round(compute_lambda_gc(p), 3)
  ggplot(d, aes(x = exp, y = obs)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60") +
    geom_point(size = 0.4, alpha = 0.5, colour = colour) +
    labs(x = "Expected -log10(p)", y = "Observed -log10(p)",
         title = glue::glue("{label} (lambda={lambda})")) +
    theme_bw() +
    theme(panel.grid = element_blank())
}

# Sex Manhattan + comparison.
sex_plink <- load_binary(sex_plink_file)
sex_ukc <- load_binary(sex_ukc_file)

sex_manhattan <- manhattan_basic(sex_ukc, "Sex GWAS (UKC/WBC)")
dir.create(supp_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(supp_dir, "sex_manhattan.pdf"), sex_manhattan, width = 8, height = 4)

sex_compare <- compare_logp(sex_plink, sex_ukc, "Sex GWAS: PLINK vs UKC")
ggsave(file.path(supp_dir, "compare_binary_logp_sex.png"), sex_compare, width = 5, height = 4, dpi = 300)

# QQ plots for Snoring, Broken bones, Risk taking.
qq_plots <- list()
colour_pair <- c("#1663A9", "#B9181A")
trait_list <- list(snoring = snoring_files, broken_bones = broken_files, risk_taking = risk_files)
for (trait_name in names(trait_list)) {
  plink_dt <- load_binary(trait_list[[trait_name]]$plink)
  ukc_dt <- load_binary(trait_list[[trait_name]]$ukc)
  p_plink <- if ("P" %in% names(plink_dt)) plink_dt$P else 10^(-plink_dt$`-log10_P`)
  p_ukc <- if ("P" %in% names(ukc_dt)) ukc_dt$P else 10^(-ukc_dt$`-log10_P`)
  qq_plots[[length(qq_plots) + 1]] <- qq_panel(p_plink, glue::glue("{trait_name} (PLINK)"), colour_pair[1])
  qq_plots[[length(qq_plots) + 1]] <- qq_panel(p_ukc, glue::glue("{trait_name} (UKC)"), colour_pair[2])
}

qq_grid <- cowplot::plot_grid(plotlist = qq_plots, ncol = 2)
ggsave(file.path(supp_dir, "binary_trait_qq_plot.png"), qq_grid, width = 10, height = 8, dpi = 300)

message("Binary calibration supplementary figures saved to ", supp_dir)
