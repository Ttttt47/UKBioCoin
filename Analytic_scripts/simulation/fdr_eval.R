#!/usr/bin/env Rscript
# --------------------------------------------------------------
# plot_eval.R      11 May 2025
# --------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(glue)
})
# Rscript plot_eval.R n5000_m10000 10000 0.05/10000
# ---------- command‑line arg: prefix ---------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript plot_eval.R <prefix> <m> <nominal_fdr>")
prefix <- args[1]
m <- as.numeric(args[2])
nominal_fdr <- as.numeric(args[3])

eval_fn <- file.path("evaluation", sprintf("eval_results_%s.rds", prefix))
eval_results <- readRDS(eval_fn)

# ---------- helper to extract FDR + Power ----------------------
conf2metrics <- function(mat) {
  # rows: pred FALSE/TRUE; cols: true 0/1
  TP <- mat["TRUE",  "1"]
  FP <- mat["TRUE",  "0"]
  FN <- mat["FALSE", "1"]
  denom_fdr   <- TP + FP
  denom_power <- TP + FN
  tibble(
    FDR   = ifelse(denom_fdr   > 0, FP / denom_fdr,   0),
    Power = ifelse(denom_power > 0, TP / denom_power, 0)
  )
}

# ---------- build tidy data frame (fixed) -----------------------
records <- lapply(names(eval_results), function(key) {
  pieces <- unlist(strsplit(key, "_", fixed = TRUE))  # <- unlist -> character
  tibble(
    trait =  pieces[1],
    maf   =  as.numeric(pieces[2]),
    h2    =  as.numeric(pieces[3]),
    mcar  =  as.numeric(pieces[4]),
    method = c("UKC", "PLINK_correct", "PLINK_wrong"),
    FDR   = NA_real_,   # placeholders
    Power = NA_real_
  ) |>                                      # fill metrics row‑wise
    group_by(method) |>
    mutate(across(c(FDR, Power), ~{
      mat <- switch(method,
        "UKC"            = eval_results[[key]]$ukc,
        "PLINK_correct"  = eval_results[[key]]$plink_correct,
        "PLINK_wrong"    = eval_results[[key]]$plink_wrong)
      if (cur_column() == "FDR") {
        TP <- mat["TRUE","1"]; FP <- mat["TRUE","0"]
        ifelse((TP+FP)>0, FP/m, 0)
      } else {
        TP <- mat["TRUE","1"]; FN <- mat["FALSE","1"]
        ifelse((TP+FN)>0, TP/(TP+FN), 0)
      }
    })) |> ungroup()
})

df <- bind_rows(records)  # all scalar columns => no list column

# tidy factors for facets
df <- df |>
  mutate(
    trait = factor(trait, levels = c("continuous","binary","ordinal")),
    maf   = factor(maf,  levels = c(0.01,0.10,0.30),
                   labels = c("MAF 0.01","MAF 0.10","MAF 0.30")),
    mcar  = factor(mcar, levels = c(0.01,0.10,0.30),
                   labels = c("MCAR 0.01","MCAR 0.10","MCAR 0.30")),
    method = factor(method,
                    levels = c("PLINK_wrong","PLINK_correct","UKC"),
                    labels = c("PLINK wrong","PLINK correct","UKC"))
  )

# ---------- plotting helper (linewidth instead of size) ---------
plot_metric <- function(trait_sel, metric_sel, ylab) {
  ggplot(df %>% filter(trait == trait_sel),
         aes(x = h2, y = .data[[metric_sel]], colour = method)) +
    geom_line(linewidth = 0.5) +              # <- linewidth
    geom_point(size = 1) +
    facet_grid(rows = vars(maf), cols = vars(mcar)) +
    scale_x_continuous(breaks = c(0.1,0.3,0.5)) +
    ylab(ylab) + xlab(expression(h^2)) +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 8)) +
    scale_colour_brewer(palette = "Set1") +
    { if (metric_sel == "FDR")
        geom_hline(yintercept = nominal_fdr, linetype = "dashed")
      else
        geom_hline(yintercept = 1.0,     linetype = "dotted")
    }
}

# ---------- generate four figures ------------------------------
dir.create(glue("plots/{prefix}"), showWarnings = FALSE)

fig_cont_fdr   <- plot_metric("continuous", "FDR",   "False Discovery Rate")
fig_cont_power <- plot_metric("continuous", "Power", "Power")
fig_bin_fdr    <- plot_metric("binary",     "FDR",   "False Discovery Rate")
fig_bin_power  <- plot_metric("binary",     "Power", "Power")
fig_ord_fdr    <- plot_metric("ordinal",     "FDR",   "False Discovery Rate")
fig_ord_power  <- plot_metric("ordinal",     "Power", "Power")

ggsave(glue("plots/{prefix}/cont_FDR.pdf"),   fig_cont_fdr,   width = 8, height = 6)
ggsave(glue("plots/{prefix}/cont_power.pdf"), fig_cont_power, width = 8, height = 6)
ggsave(glue("plots/{prefix}/bin_FDR.pdf"),    fig_bin_fdr,    width = 8, height = 6)
ggsave(glue("plots/{prefix}/bin_power.pdf"),  fig_bin_power,  width = 8, height = 6)
ggsave(glue("plots/{prefix}/ord_FDR.pdf"),    fig_ord_fdr,    width = 8, height = 6)
ggsave(glue("plots/{prefix}/ord_power.pdf"),  fig_ord_power,  width = 8, height = 6)


# extract legend once
legend_shared <- get_legend(
  fig_cont_fdr + theme(legend.position = "bottom", legend.title = element_text(size=12), 
                        legend.text = element_text(size=12))
)

# drop legends from individual plots
figs_noleg <- lapply(
  list(fig_cont_fdr, fig_cont_power, fig_bin_fdr, fig_bin_power, fig_ord_fdr, fig_ord_power),
  function(p) p + theme(legend.position = "none")
)

# create row labels
row_label <- function(text) {
  ggdraw() +
    draw_label(text, angle = 90, fontface = "bold", x = 0.5, y = 0.5) +
    theme_void()
}
label_cont <- row_label("Continuous trait")
label_bin  <- row_label("Binary trait")
label_ord  <- row_label("Ordinal trait")

# assemble rows with labels in first column
row1 <- plot_grid(
  label_cont, 
  figs_noleg[[1]], figs_noleg[[2]],
  ncol = 3,
  rel_widths = c(0.1, 1, 1)
)
row2 <- plot_grid(
  label_bin,  
  figs_noleg[[3]], figs_noleg[[4]],
  ncol = 3,
  rel_widths = c(0.1, 1, 1)
)
row3 <- plot_grid(
  label_ord,  
  figs_noleg[[5]], figs_noleg[[6]],
  ncol = 3,
  rel_widths = c(0.1, 1, 1)
)

# combine rows
panel_g <- plot_grid(
  row1, row2, row3,
  ncol = 1,
  rel_heights = c(1,1,1)
)

# title
title_g <- ggdraw() +
  draw_label(
    "UKC vs PLINK: FDR and Power across Heritability, MAF and MCAR",
    fontface = "bold", size = 14
  )

# final layout: title, panels, legend
final_g <- plot_grid(
  title_g,
  panel_g,
  legend_shared,
  ncol = 1,
  rel_heights = c(0.06, 1, 0.04)
)

# save

pdf(glue("plots/{prefix}/all_metrics.pdf"),
    width = 9, height = 10)
final_g
dev.off()


cat("Saved figures to 'plots/' directory\n")
