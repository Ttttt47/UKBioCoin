# Figure 1 (panels A–C)

library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(grid)

# ---- Edit these paths to your files ----
allele_freq_file <- "data/0.basic/0.freq.afreq"
cov_yy_file <- "data/2.matrix/10M_cov_yy.table" # full matrix; subset handled below
xld_file <- "data/XLD/UKB_10M.X-LD-Plus.txt"
ukb_pcs_file <- "data/pcs/22009.tab" # IID, PC1, PC2, ...
ukc_pcs_file <- "data/pcs/10M_PCs.eigenvec" # IID, PC1, PC2, ...
cov_slice <- 1:520 # rows/cols to keep for the 505x505 heatmap
out_file <- file.path("fig", "fig1_basic.pdf")

# ---- Panel A: MAF distribution ----
allele_freq <- fread(allele_freq_file)
allele_freq[, MAF := pmin(ALT_FREQS, 1 - ALT_FREQS)]

p_maf <- ggplot(allele_freq, aes(x = MAF)) +
  geom_histogram(binwidth = 0.005, fill = "grey") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    legend.position = "none",
    axis.text = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 12, colour = "black")
  ) +
  xlab("MAF") +
  ylab("Count")

ggsave("0.freq.maf.pdf", p_maf, width = 9, height = 8, dpi = 300, units = "cm")

# ---- Panel B: chromosome-wise LD vs inverse length ----
xld <- fread(xld_file, sep = "\t")
xld <- xld[chri == chrj]
xld[, inv_chr_len := 1 / marker.number]
plot_dt <- xld[, .(x = inv_chr_len, LD, chri)]

lm_fit <- lm(LD ~ inv_chr_len, data = xld)
r_label <- cor(xld$LD, xld$inv_chr_len, use = "complete.obs")

p_ld <- ggplot(plot_dt, aes(x = x, y = LD)) +
  geom_point(color = "#357EBD99", size = 6) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "#357EBD99",
    fill = "#ADB6B699"
  ) +
  geom_text(aes(label = chri), family = "serif") +
  scale_y_continuous(expand = c(0, 0), limits = c(1e-4, 8e-4)) +
  scale_x_continuous(expand = c(0, 0), limits = c(1e-6, 9.8e-06)) +
  labs(x = "1/Number of SNPs", y = "LD") +
  theme_minimal() +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(8, 8, 4, 4),
    panel.background = element_rect(fill = "transparent", colour = NA),
    panel.grid = element_blank(),
    text = element_text(colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    axis.text = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 12, colour = "black")
  )

ggsave("0.UKB_10M.XLD.pdf", p_ld, width = 9, height = 8, units = "cm", dpi = 300)

# ---- Panel C: 505x505 phenotype correlation heatmap ----
cov_mat <- as.matrix(read.table(cov_yy_file))
cov_sub <- cov_mat[cov_slice, cov_slice, drop = FALSE]

png("0.cov_yy.png", width = 11, height = 8, units = "cm", res = 300)
Heatmap(
  cov_sub,
  name = "Correlation",
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_title = "Phenotypes",
  column_title = "Phenotypes",
  width = unit(7, "cm"),
  height = unit(7, "cm"),
  row_names_side = "left",
  column_names_side = "top",
  column_names_gp = gpar(fontsize = 12),
  row_names_gp = gpar(fontsize = 12),
  heatmap_legend_param = list(
    title = "Correlation",
    at = seq(-1, 1, 0.2),
    labels = seq(-1, 1, 0.2),
    title_gp = gpar(fontsize = 12),
    labels_gp = gpar(fontsize = 12)
  )
)
dev.off()
