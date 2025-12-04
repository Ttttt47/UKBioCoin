# Figure 6 (Mendelian randomization).


library(data.table)
library(ggplot2)

# ---- Edit these paths ----
merged_dir <- "data/mr/merged"  

if (file.exists(file.path(merged_dir, "snp_data.csv")) &&
    file.exists(file.path(merged_dir, "mr_res.csv"))) {
  point_data <- fread(file.path(merged_dir, "snp_data.csv"))
  res_data   <- fread(file.path(merged_dir, "mr_res.csv"))

  p_scatter <- ggplot(point_data, aes(x = beta.exposure, y = beta.outcome)) +
    geom_point(color = "black", size = 2, alpha = 0.6) +
    geom_hline(yintercept = 0, color = "red", alpha = 0.2) +
    geom_vline(xintercept = 0, color = "red", alpha = 0.2) +
    geom_abline(
      data = res_data,
      aes(intercept = Intercept, slope = Estimate, color = Method),
      linetype = "solid",
      show.legend = TRUE,
      linewidth = 1
    ) +
    xlab("SNP effect on exposure") +
    ylab("SNP effect on outcome") +
    theme(
      plot.title = element_text(size = rel(1.2), face = "bold"),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(colour = "gray80", linetype = "dotted"),
      panel.grid.minor = element_line(colour = "gray90", linetype = "dotted"),
      legend.key = element_rect(fill = "white"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
    )
  dir.create(dirname(out_scatter), showWarnings = FALSE, recursive = TRUE)
  ggsave(out_scatter, p_scatter, width = 6, height = 5)
  message("Saved MR scatter: ", out_scatter)
}
