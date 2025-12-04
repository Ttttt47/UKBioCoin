# Figure 2 reproduction script (Missingness QC).

library(data.table)
library(ggplot2)
library(cowplot)
library(glue)
library(egg)
library(grid)

# ---- Edit this working directory to your local path ----
setwd("path/to/3.analysis/reproducing_GWAS/near_com/")

phe <- "X20015.0.0"
covariates <- c("None", "5 PCs", "10 PCs")
color_plates <- c("#0071C2", "#614099", "#369F2D")

for (PC_num in c(0, 5, 10)) {
  eval(parse(text = glue("plink_res_{PC_num}pc <- fread('plink-{phe}-{PC_num}PC.{phe}.glm.linear')")))
  eval(parse(text = glue("plink_res_{PC_num}pc <- plink_res_{PC_num}pc[plink_res_{PC_num}pc$TEST == 'ADD', ]")))
  eval(parse(text = glue("ukc_res_{PC_num}pc   <- fread('UKC-{phe}-{PC_num}PC_results.table')")))
  print(PC_num)
}

# ---- Standing height ----
options(digits = 5)
plot_list <- list()

for (i in 1:3) {
  PC_num <- c(0, 5, 10)[i]

  eval(parse(text = glue("VIF_filtered_flag <- ukc_res_{PC_num}pc$VIF < 50")))
  ids <- (1:dim(ukc_res_0pc)[1])[VIF_filtered_flag]

  plink_points <- eval(parse(text = glue("plink_res_{PC_num}pc$BETA[ids]")))
  ukc_points <- eval(parse(text = glue("ukc_res_{PC_num}pc$BETA[ids]")))

  lr <- lm(ukc_points ~ plink_points)
  cor_val <- cor(plink_points, ukc_points, use = "complete.obs")
  cor_val <- round(cor_val, digits = 3)

  if (i == 3) {
    VIF_filtered_flag <- which(
      ukc_res_10pc$VIF > 50 &
        abs(ukc_res_10pc$BETA - plink_res_10pc$BETA) > 0.2
    )
    highlight_outliers <- setdiff(
      which(abs(ukc_res_10pc$BETA - plink_res_10pc$BETA) > 0.01),
      VIF_filtered_flag
    )
    ids <- setdiff(ids, union(VIF_filtered_flag, highlight_outliers))
    p <- ggplot() +
      geom_point(
        data = data.frame(
          x = plink_res_10pc$BETA[highlight_outliers],
          y = ukc_res_10pc$BETA[highlight_outliers]
        ),
        aes(x = x, y = y),
        color = "red",
        size = 0.3,
        alpha = 0.2
      )

    # highlight VIF > 50 & bias > 0.2 points with a different shape/size
    VIF_filtered_flag <- which(
      ukc_res_10pc$VIF > 50 &
        abs(ukc_res_10pc$BETA - plink_res_10pc$BETA) > 0.2
    )
    p <- p +
      geom_point(
        data = data.frame(
          x = plink_res_10pc$BETA[VIF_filtered_flag],
          y = ukc_res_10pc$BETA[VIF_filtered_flag]
        ),
        aes(x = x, y = y),
        color = "red",
        size = 1.5,
        alpha = 1,
        shape = 2
      )
  } else {
    p <- ggplot()
  }

  # plot only a subset of points with p-value < 0.05
  normal_ids <- intersect(
    which(eval(parse(text = glue("plink_res_{PC_num}pc$P"))) > 0.05),
    which(eval(parse(text = glue("ukc_res_{PC_num}pc$`-log10_P`"))) < -log10(0.05))
  )
  ids <- setdiff(ids, sample(normal_ids, length(normal_ids) * 0.95))
  plink_points <- eval(parse(text = glue("plink_res_{PC_num}pc$BETA[ids]")))
  ukc_points <- eval(parse(text = glue("ukc_res_{PC_num}pc$BETA[ids]")))

  p <- p +
    scale_x_continuous(
      name = expression("UKB estimates of " * beta),
      limits = c(-0.5, 0.5)
    ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linewidth = 10,
      color = "grey",
      alpha = 0.3
    ) +
    geom_abline(
      intercept = lr$coeff[[1]],
      slope = lr$coeff[[2]],
      linewidth = 0.5,
      color = "black",
      alpha = 1,
      linetype = "dashed"
    ) +
    geom_point(
      data = data.frame(x = plink_points, y = ukc_points),
      aes(x = x, y = y),
      color = color_plates[i],
      size = 0.3,
      alpha = 0.4
    ) +
    theme(
      panel.background = element_rect(fill = "transparent", colour = NA),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.5),
      axis.ticks = element_line(colour = "black", linewidth = 0.5)
    ) +
    ggtitle(glue("Covariates: {covariates[i]}"))

  p <- p + scale_y_continuous(
    name = expression("UKC estimates of " * beta),
    limits = c(-0.5, 0.5)
  )

  print(glue("PC{i}, y = {lr$coeff[[2]]} x + {lr$coeff[[1]]}, cor = {cor_val}"))
  plot_list[[length(plot_list) + 1]] <- p
  ggsave(glue("nearcom_beta_compare_{i}.png"), height = 900, width = 900, units = "px")
}

# Comparing PLINK and UKC results: -log10(p)
for (i in 1:3) {
  PC_num <- c(0, 5, 10)[i]

  eval(parse(text = glue("VIF_filtered_flag <- ukc_res_{PC_num}pc$VIF < 50")))
  ids <- (1:dim(ukc_res_0pc)[1])[VIF_filtered_flag]
  plink_points <- eval(parse(text = glue("-log(plink_res_{PC_num}pc$P[ids], 10)")))
  ukc_points <- eval(parse(text = glue("ukc_res_{PC_num}pc$`-log10_P`[ids]")))
  cor_val <- cor(plink_points, ukc_points, use = "complete.obs")
  cor_val <- round(cor_val, digits = 3)
  lr <- lm(ukc_points ~ plink_points)

  if (i == 3) {
    VIF_filtered_flag <- which(
      ukc_res_10pc$VIF > 50 &
        abs(ukc_res_10pc$BETA - plink_res_10pc$BETA) > 0.2
    )
    highlight_outliers <- setdiff(
      which(abs(ukc_res_10pc$BETA - plink_res_10pc$BETA) > 0.01),
      VIF_filtered_flag
    )
    ids <- setdiff(ids, union(VIF_filtered_flag, highlight_outliers))
    p <- ggplot() +
      geom_point(
        data = data.frame(
          x = -log(plink_res_10pc$P[highlight_outliers], 10),
          y = ukc_res_10pc$`-log10_P`[highlight_outliers]
        ),
        aes(x = x, y = y),
        color = "red",
        size = 0.3,
        alpha = 0.2
      ) +
      geom_point(
        data = data.frame(
          x = -log(plink_res_10pc$P[VIF_filtered_flag], 10),
          y = ukc_res_10pc$`-log10_P`[VIF_filtered_flag]
        ),
        aes(x = x, y = y),
        color = "red",
        size = 1.5,
        alpha = 1,
        shape = 2
      )
  } else {
    p <- ggplot()
  }

  normal_ids <- intersect(
    which(eval(parse(text = glue("plink_res_{PC_num}pc$P"))) > 0.05),
    which(eval(parse(text = glue("ukc_res_{PC_num}pc$`-log10_P`"))) < -log10(0.05))
  )
  ids <- setdiff(ids, sample(normal_ids, length(normal_ids) * 0.95))
  plink_points <- eval(parse(text = glue("-log(plink_res_{PC_num}pc$P[ids], 10)")))
  ukc_points <- eval(parse(text = glue("ukc_res_{PC_num}pc$`-log10_P`[ids]")))

  p <- p +
    geom_abline(
      intercept = 0,
      slope = 1,
      linewidth = 10,
      color = "grey",
      alpha = 0.3
    ) +
    geom_abline(
      intercept = lr$coeff[[1]],
      slope = lr$coeff[[2]],
      linewidth = 0.5,
      color = "black",
      alpha = 1,
      linetype = "dashed"
    ) +
    geom_point(
      data = data.frame(x = plink_points, y = ukc_points),
      aes(x = x, y = y),
      color = color_plates[i],
      size = 0.3,
      alpha = 0.2
    ) +
    scale_x_continuous(
      name = expression("UKB estimates of " * -log[10] * "(p)"),
      limits = c(0, 100)
    ) +
    theme(
      panel.background = element_rect(fill = "transparent", colour = NA),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.5),
      axis.ticks = element_line(colour = "black", linewidth = 0.5)
    ) +
    ggtitle(glue("Covariates: {covariates[i]}"))

  p <- p + scale_y_continuous(
    name = expression("UKC estimates of " * -log[10] * "(p)"),
    limits = c(0, 100)
  )

  print(glue("PC{i}, y = {lr$coeff[[2]]} x + {lr$coeff[[1]]}, cor = {cor_val}"))
  ggsave(glue("nearcom_logp_compare_{i}.png"), height = 900, width = 900, units = "px")
  plot_list[[length(plot_list) + 1]] <- p
}

jpeg(
  filename = "trueBeta_and_pvalues_vs_UKC_height.png",
  height = 1800,
  width = 900 * 3,
  quality = 100,
  res = 300
)
plot_grid(plotlist = plot_list, nrow = 2, align = "hv", axis = "lr")
dev.off()

# plot VIF vs beta_bias, and the regression line, denote the correlation and the regression equation
ids <- which(VIF_logp$beta_bias > 0.001)
highlight <- which(VIF_logp$beta_bias > 0.2 & VIF_logp$VIF > 50)

p <- ggplot(VIF_logp[setdiff(ids, highlight), ], aes(x = VIF, y = beta_bias)) +
  geom_point() +
  geom_point(
    data = VIF_logp[highlight, ],
    color = "red",
    shape = 2,
    size = 3
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "grey",
    alpha = 0.5
  ) +
  labs(x = "VIF", y = expression("Bias of " * beta)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    legend.position = "none",
    axis.text = element_text(size = 9, colour = "black"),
    axis.title = element_text(size = 11, colour = "black")
  )

ggsave(
  filename = "VIF_vs_beta_bias.png",
  plot = p,
  width = 1200,
  height = 1800,
  units = "px"
)

# ---- Manhattan helper ----
my_manhattan <- function(
    x,
    chr = "CHR",
    bp = "BP",
    p = "P",
    snp = "SNP",
    col = c("gray10", "gray60"),
    chrlabs = NULL,
    suggestiveline = -log10(1e-05),
    genomewideline = -log10(5e-08),
    highlight = NULL,
    logp = TRUE,
    highlight.color = "green3",
    highlight.pch = 20,
    point.cex = 1,
    point.pch = 20,
    annotatePval = NULL,
    annotateTop = TRUE,
    ...
  ) {
  CHR <- BP <- P <- index <- NULL
  if (!(chr %in% names(x))) {
    stop(paste("Column", chr, "not found!"))
  }
  if (!(bp %in% names(x))) {
    stop(paste("Column", bp, "not found!"))
  }
  if (!(p %in% names(x))) {
    stop(paste("Column", p, "not found!"))
  }
  if (!(snp %in% names(x))) {
    warning(paste("No SNP column found. OK unless you're trying to highlight."))
  }
  if (!is.numeric(x[[chr]])) {
    stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  }
  if (!is.numeric(x[[bp]])) {
    stop(paste(bp, "column should be numeric."))
  }
  if (!is.numeric(x[[p]])) {
    stop(paste(p, "column should be numeric."))
  }

  if (!is.null(x[[snp]])) {
    d <- data.frame(
      CHR = x[[chr]],
      BP = x[[bp]],
      P = x[[p]],
      pos = NA,
      index = NA,
      SNP = x[[snp]],
      stringsAsFactors = FALSE
    )
  } else {
    d <- data.frame(
      CHR = x[[chr]],
      BP = x[[bp]],
      P = x[[p]],
      pos = NA,
      index = NA
    )
  }

  d <- d[order(d$CHR, d$BP), ]
  if (logp) {
    d$logp <- -log10(d$P)
  } else {
    d$logp <- d$P
  }

  d$index <- rep.int(seq_along(unique(d$CHR)), times = tapply(d$SNP, d$CHR, length))
  nchr <- length(unique(d$CHR))

  if (nchr == 1) {
    d$pos <- d$BP
    xlabel <- paste("Chromosome", unique(d$CHR), "position")
  } else {
    lastbase <- 0
    ticks <- NULL
    for (i in unique(d$index)) {
      if (i == 1) {
        d[d$index == i, ]$pos <- d[d$index == i, ]$BP
      } else {
        lastbase <- lastbase + max(d[d$index == (i - 1), "BP"])
        d[d$index == i, "BP"] <- d[d$index == i, "BP"] - min(d[d$index == i, "BP"]) + 1
        d[d$index == i, "pos"] <- d[d$index == i, "BP"] + lastbase
      }
    }
    ticks <- tapply(d$pos, d$index, quantile, probs = 0.5)
    xlabel <- "Chromosome"
    labs <- unique(d$CHR)
  }

  xmax <- ceiling(max(d$pos) * 1.03)
  xmin <- floor(max(d$pos) * -0.03)
  def_args <- list(
    xaxt = "n",
    bty = "n",
    xaxs = "i",
    yaxs = "i",
    las = 1,
    pch = 20,
    xlim = c(xmin, xmax),
    ylim = c(0, ceiling(max(d$logp))),
    xlab = xlabel,
    ylab = expression(-log[10](italic(p)))
  )
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))

  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(labs)) {
        labs <- chrlabs
      } else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    } else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }

  if (nchr == 1) {
    axis(1, ...)
  } else {
    axis(1, at = ticks, labels = labs, ...)
  }

  col <- rep_len(col, max(d$index))

  if (nchr == 1) {
    with(d, points(pos, logp, pch = point.pch, col = col[1], cex = point.cex, ...))
  } else {
    icol <- 1
    for (i in unique(d$index)) {
      points(
        d[d$index == i, "pos"],
        d[d$index == i, "logp"],
        col = col[icol],
        pch = point.pch,
        cex = point.cex,
        ...
      )
      icol <- icol + 1
    }
  }

  if (suggestiveline) {
    abline(h = suggestiveline, col = "blue")
  }
  if (genomewideline) {
    abline(h = genomewideline, col = "red")
  }

  if (!is.null(highlight)) {
    if (any(!(highlight %in% d$SNP))) {
      warning("You're trying to highlight SNPs that don't exist in your results.")
    }
    d.highlight <- d[which(d$SNP %in% highlight), ]
    with(d.highlight, points(pos, logp, col = highlight.color, pch = highlight.pch, ...))
  }

  if (!is.null(annotatePval)) {
    if (logp) {
      topHits <- subset(d, P <= annotatePval)
    } else {
      topHits <- subset(d, P >= annotatePval)
    }

    par(xpd = TRUE)
    if (annotateTop == FALSE) {
      if (logp) {
        with(
          subset(d, P <= annotatePval),
          textxy(pos, -log10(P), offset = 0.625, labs = topHits$SNP, cex = 0.45),
          ...
        )
      } else {
        with(
          subset(d, P >= annotatePval),
          textxy(pos, P, offset = 0.625, labs = topHits$SNP, cex = 0.45),
          ...
        )
      }
    } else {
      topHits <- topHits[order(topHits$P), ]
      topSNPs <- NULL
      for (i in unique(topHits$CHR)) {
        chrSNPs <- topHits[topHits$CHR == i, ]
        topSNPs <- rbind(topSNPs, chrSNPs[1, ])
      }
      if (logp) {
        textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625, labs = topSNPs$SNP, cex = 0.5, ...)
      } else {
        textxy(topSNPs$pos, topSNPs$P, offset = 0.625, labs = topSNPs$SNP, cex = 0.5, ...)
      }
    }
  }
  par(xpd = FALSE)
}

# ---- VIF vs beta_bias miami plot ----
VIF_bias <- data.frame(
  SNP = ukc_res_10pc$ID,
  POS = ukc_res_10pc$POS,
  CHR = ukc_res_10pc$`#CHROM`,
  beta_bias = abs(ukc_res_10pc$BETA - plink_res_10pc$BETA),
  VIF = ukc_res_10pc$VIF
)

VIF_bias_selected <- VIF_bias[which(abs(ukc_res_10pc$BETA - plink_res_10pc$BETA) > 0.001 | VIF_bias$VIF > 1.01), ]

jpeg("./miami_VIF_vs_bias.jpg", height = 1200, width = 3600, quality = 100, res = 400)
par(mfrow = c(2, 1))
par(mar = c(2.5, 5, 1, 1))
my_manhattan(
  data.frame(
    CHR = VIF_bias_selected$CHR,
    BP = VIF_bias_selected$POS,
    SNP = VIF_bias_selected$SNP,
    P = VIF_bias_selected$beta_bias
  ),
  ylim = c(0, 0.35),
  highlight = VIF_bias_selected$SNP[which(VIF_bias_selected$beta_bias > 0.2 & VIF_bias_selected$VIF > 50)],
  highlight.color = "red",
  highlight.pch = 2,
  logp = FALSE,
  ylab = expression("Bias of " * beta)
)

par(mar = c(1, 5, 0, 1))
my_manhattan(
  data.frame(
    CHR = VIF_bias_selected$CHR,
    BP = VIF_bias_selected$POS,
    SNP = VIF_bias_selected$SNP,
    P = VIF_bias_selected$VIF
  ),
  ylim = c(260, 0),
  highlight = VIF_bias_selected$SNP[which(VIF_bias_selected$beta_bias > 0.2 & VIF_bias_selected$VIF > 50)],
  highlight.color = "red",
  highlight.pch = 2,
  xlab = "",
  logp = FALSE,
  ylab = "VIF",
  genomewideline = 50,
  suggestiveline = FALSE,
  xaxt = "n"
)
dev.off()

# ---- Neuroticism score ----
setwd("/public3/project_users/chengb/hjc/projects/UKBioCoin/0729/3.analysis/reproducing_GWAS/high_miss/")

phe <- "X20127.0.0"
color_plates <- c("#0071C2", "#369F2D", "#EDB11A", "#D75615")
covariates <- c("5 PCs", "5 PCs, C1", "5 PCs, C1, C2", "5 PCs, C1, C2, C3")

for (PC_num in c(5)) {
  for (i in 0:3) {
    eval(parse(text = glue("plink_res_{PC_num}pc_{i}=fread('./plink-{phe}-{PC_num}PC-{i}.{phe}.glm.linear')")))
    eval(parse(text = glue("plink_res_{PC_num}pc_{i}=plink_res_{PC_num}pc_{i}[plink_res_{PC_num}pc_{i}$TEST=='ADD',]")))
    eval(parse(text = glue("ukc_res_{PC_num}pc_{i}=fread('./ukc-{phe}-{PC_num}PC-{i}_results.table')")))
    print(i)
  }
}

plot_list <- list()
for (i in c(0:3)) {
  PC_num <- 5
  eval(parse(text = glue("VIF_filtered_flag = ukc_res_{PC_num}pc_{i}$VIF<50")))
  ids <- (1:dim(ukc_res_5pc_1)[1])[VIF_filtered_flag]
  plink_points <- eval(parse(text = glue("plink_res_{PC_num}pc_{i}$BETA[ids]")))
  ukc_points <- eval(parse(text = glue("ukc_res_{PC_num}pc_{i}$BETA[ids]")))
  cor <- cor(plink_points, ukc_points, use = "complete.obs")
  cor <- round(cor, digits = 3)
  lr <- lm(ukc_points ~ plink_points)

  p <- ggplot()

  # plot only the points with p-value < 0.05 to save time.
  normal_ids <- intersect(
    which(eval(parse(text = glue("plink_res_{PC_num}pc_{i}$P"))) > 0.05),
    which(eval(parse(text = glue("ukc_res_{PC_num}pc_{i}$`-log10_P`"))) < -log10(0.05))
  )
  ids <- setdiff(ids, sample(normal_ids, length(normal_ids) * 0.8))
  plink_points <- eval(parse(text = glue("plink_res_{PC_num}pc_{i}$BETA[ids]")))
  ukc_points <- eval(parse(text = glue("ukc_res_{PC_num}pc_{i}$BETA[ids]")))

  p <- p +
    scale_x_continuous(name = expression("UKB estimates of " * beta), limits = c(-0.5, 0.5)) +
    geom_abline(intercept = 0, slope = 1, linewidth = 10, color = "grey", alpha = 0.3) +
    geom_abline(intercept = lr$coeff[[1]], slope = lr$coeff[[2]], linewidth = 0.5, color = "black", alpha = 1, linetype = "dashed") +
    geom_point(
      data = data.frame(x = plink_points, y = ukc_points),
      aes(x = x, y = y),
      color = color_plates[i + 1],
      size = 0.3,
      alpha = 0.4
    ) +
    theme(
      panel.background = element_rect(fill = "transparent", colour = NA),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.5),
      axis.ticks = element_line(colour = "black", linewidth = 0.5)
    ) +
    ggtitle(glue("Covariates: {covariates[i + 1]}")) +
    theme(plot.title = element_text(size = 12))

  p <- p + scale_y_continuous(name = expression("UKC estimates of " * beta), limits = c(-0.5, 0.5))

  print(glue("PC{i},y={lr$coeff[[2]]}x+{lr$coeff[[1]]},cor={cor}"))
  plot_list[[length(plot_list) + 1]] <- p
  ggsave(glue("highmiss_beta_compare_{i}.png"), height = 900, width = 900, units = "px")
}

for (i in c(0:3)) {
  PC_num <- 5
  eval(parse(text = glue("VIF_filtered_flag = ukc_res_{PC_num}pc_{i}$VIF<50")))
  ids <- (1:dim(ukc_res_5pc_1)[1])[VIF_filtered_flag]
  plink_points <- eval(parse(text = glue("-log(plink_res_{PC_num}pc_{i}$P[ids],10)")))
  ukc_points <- eval(parse(text = glue("ukc_res_{PC_num}pc_{i}$`-log10_P`[ids]")))
  cor <- cor(plink_points, ukc_points, use = "complete.obs")
  cor <- round(cor, digits = 3)
  lr <- lm(ukc_points ~ plink_points)

  p <- ggplot()

  # plot only the points with p-value < 0.05 to save time.
  normal_ids <- intersect(
    which(eval(parse(text = glue("plink_res_{PC_num}pc_{i}$P"))) > 0.05),
    which(eval(parse(text = glue("ukc_res_{PC_num}pc_{i}$`-log10_P`"))) < -log10(0.05))
  )
  ids <- setdiff(ids, sample(normal_ids, length(normal_ids) * 0.8))
  plink_points <- eval(parse(text = glue("-log(plink_res_{PC_num}pc_{i}$P[ids],10)")))
  ukc_points <- eval(parse(text = glue("ukc_res_{PC_num}pc_{i}$`-log10_P`[ids]")))

  p <- p +
    geom_abline(intercept = 0, slope = 1, linewidth = 10, color = "grey", alpha = 0.3) +
    geom_abline(intercept = lr$coeff[[1]], slope = lr$coeff[[2]], linewidth = 0.5, color = "black", alpha = 1, linetype = "dashed") +
    geom_point(
      data = data.frame(x = plink_points, y = ukc_points),
      aes(x = x, y = y),
      color = color_plates[i + 1],
      size = 0.3,
      alpha = 0.4
    ) +
    theme(
      panel.background = element_rect(fill = "transparent", colour = NA),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.5),
      axis.ticks = element_line(colour = "black", linewidth = 0.5)
    ) +
    scale_x_continuous(name = expression("UKB estimates of " * -log[10] * "(p)"), limits = c(0, 25)) +
    ggtitle(glue("Covariates: {covariates[i + 1]}")) +
    theme(plot.title = element_text(size = 12))

  p <- p + scale_y_continuous(name = expression("UKC estimates of " * -log[10] * "(p)"), limits = c(0, 25))

  print(glue("PC{i},y={lr$coeff[[2]]}x+{lr$coeff[[1]]},cor={cor}"))
  ggsave(glue("highmiss_logp_compare_{i}.png"), height = 900, width = 900, units = "px")
  plot_list[[length(plot_list) + 1]] <- p
}

jpeg("./trueBeta&pvalues.vs.ukc_neuro.png", height = 1800, width = 900 * 4, quality = 100, res = 300)
plot_grid(plotlist = plot_list, nrow = 2, align = "hv", axis = "lr")
dev.off()
