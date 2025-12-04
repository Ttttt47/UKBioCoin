# Supplementary WBBC validation.
# Reproduces:
#   SuppFig/FigS_WBC.pdf


library(data.table)
library(ggplot2)
library(cowplot)

plink_file <- "data/wbbc/height_plink.table"
ukc_file <- "data/wbbc/height_ukc.table"
supp_dir <- "SuppFig"

load_gwas <- function(path) {
  dt <- fread(path)
  if ("TEST" %in% names(dt)) dt <- dt[TEST %in% c("ADD", "ADD_BOLT", "ADDitive", "ADD_ALT")]
  if (!("-log10_P" %in% names(dt)) && "P" %in% names(dt)) {
    dt[, `-log10_P` := -log10(P)]
  }
  dt
}

plink_dt <- load_gwas(plink_file)
ukc_dt <- load_gwas(ukc_file)

id_col <- if ("ID" %in% names(plink_dt) && "ID" %in% names(ukc_dt)) "ID" else "SNP"
merged <- merge(plink_dt, ukc_dt, by = id_col, suffixes = c("_plink", "_ukc"))

beta_plot <- ggplot(merged, aes(x = BETA_plink, y = BETA_ukc)) +
  geom_point(alpha = 0.4, size = 0.6, colour = "#1663A9") +
  geom_abline(intercept = 0, slope = 1, colour = "grey60") +
  labs(x = "PLINK beta", y = "UKC beta",
       title = "WBBC estimation of beta") +
  theme_bw() +
  theme(panel.grid = element_blank())

logp_plot <- ggplot(merged, aes(x = `-log10_P_plink`, y = `-log10_P_ukc`)) +
  geom_point(alpha = 0.4, size = 0.6, colour = "#1663A9") +
  geom_abline(intercept = 0, slope = 1, colour = "grey60") +
  labs(x = "PLINK -log10(p)", y = "UKC -log10(p)",
       title = "WBBC estimation of -log10(p)") +
  theme_bw() +
  theme(panel.grid = element_blank())

fig <- cowplot::plot_grid(beta_plot, logp_plot, ncol = 2)
dir.create(supp_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(supp_dir, "WBC.pdf")
cowplot::save_plot(out_file, fig, base_width = 10, base_height = 4)
message("Saved WBBC validation figure to ", out_file)
