#!/usr/bin/env Rscript
# run_simulation_pipeline.R
# End-to-end: simulate G/C/Y, build PLINK pfile, run GWAS.

# ---- 1. Dependencies ----
if (!requireNamespace("data.table", quietly=TRUE)) {
  install.packages("data.table")
}
library(data.table)
library(genio)

# ---- 2. Simulation function (generate_GxCyY) ----
generate_GxCyY <- function(
    n        = 5000,
    m        = 20000,
    s_y      = 50,
    s_c2     = 12,
    s_c3     = 200,
    maf_sig  = 0.30,
    maf_bg   = c(0.05,0.50),
    h2       = 0.3,
    R2_c2    = 0.25,
    R2_c3    = 0.05,
    gamma    = c(0.4,0.4,0.4),
    trait    = c("continuous","binary", "ordinal"),
    K_ord = 5,
    prevalence = 0.10,
    mcar     = 0.20              # overall MCAR rate
) {
  trait <- match.arg(trait)

  # index sets
  stopifnot(s_y + s_c2 + s_c3 < m)
  S_Y  <- 1:s_y
  S_C2 <- (s_y+1):(s_y+s_c2)
  S_C3 <- (s_y+s_c2+1):(s_y+s_c2+s_c3)

  # MAFs and genotypes
  maf <- numeric(m)
  maf[S_Y]  <- maf_sig
  maf[-S_Y] <- runif(m-s_y, maf_bg[1], maf_bg[2])
  G <- sapply(maf, function(p) rbinom(n, 2, p))
  storage.mode(G) <- "integer"

  # Covariates
  C1 <- rnorm(n)
  # C2 few-strong
  z2 <- scale(G[,S_C2,drop=FALSE]) %*% rnorm(s_c2)
  sd2 <- sqrt((1-R2_c2)/R2_c2 * var(z2))
  C2 <- scale(z2 + rnorm(n, 0, sd2))
  # C3 many-weak
  z3 <- scale(G[,S_C3,drop=FALSE]) %*% rnorm(s_c3, sd=1/5)
  sd3 <- as.numeric(sqrt((1-R2_c3)/R2_c3 * var(z3)))
  C3 <- scale(z3 + rnorm(n, 0, sd3))

  # Phenotype
  beta  <- rnorm(s_y)
  zy    <- scale(G[,S_Y,drop=FALSE]) %*% beta
  beta  <- beta * as.numeric(sqrt(h2 / var(zy)))
  zy    <- scale(G[,S_Y,drop=FALSE]) %*% beta

  if (trait == "continuous") {
    var_cov <- sum(gamma^2 * c(var(C1), var(C2), var(C3)))
    sd_e    <- sqrt(1 - h2 - var_cov)
    Y       <- as.numeric(zy + gamma[1]*C1 + gamma[2]*C2 + gamma[3]*C3 +
                         rnorm(n,0,sd_e))
    Y <- scale(Y)
  } else if (trait == "binary") {
    # 1) Compute total variance from covariates
    var_cov <- sum(gamma^2 * c(var(C1), var(C2), var(C3)))
    # 2) At this point we've already scaled beta so that
    #    Var(zy) == h2, i.e. genetic variance = h2
    # 3) Residual SD so that total latent variance = 1
    sd_eps <- sqrt(1 - h2 - var_cov)

    # 4) Simulate latent liability (mean 0, var 1)
    L <- as.vector(
      zy +
      gamma[1]*C1 +
      gamma[2]*C2 +
      gamma[3]*C3 +
      rnorm(n, 0, sd_eps)
    )
    # 5) Threshold at the (1 – prevalence) quantile
    thr <- qnorm(1 - prevalence, mean = 0, sd = 1)
    Y   <- as.integer(L > thr)  # 1 = case, 0 = control
  } else {
      var_cov <- sum(gamma^2 * c(var(C1), var(C2), var(C3)))
      sd_eps  <- sqrt(1 - h2 - var_cov)
      L <- as.vector(zy + gamma[1]*C1 + gamma[2]*C2 + gamma[3]*C3 +
                     rnorm(n, 0, sd_eps))
      L <- scale(L)                         # mean 0, var 1

      probs <- seq(0, 1, length.out = K_ord + 1)
      cuts  <- quantile(L, probs, names = FALSE)
      Y     <- as.integer(cut(L, breaks = cuts, include.lowest = TRUE))
  }

  # MCAR missingness for each covariate
  miss   <- matrix(runif(n*3) < mcar, nrow=n, ncol=3)
  Cmat   <- cbind(C1, C2, C3)
  Cmat[miss] <- NA
  colnames(Cmat) <- c("C1","C2","C3")

  list(
    G              = G,
    maf            = maf,
    sets           = list(Y=S_Y, C2=S_C2, C3=S_C3),
    beta           = beta,
    covariates     = Cmat,
    Y              = as.numeric(Y)
  )
}

write_plink_pfile_via_bigsnpr <- function(G, prefix = "sim_bigsnpr") {
  # 1) install/load bigsnpr if needed
  if (!requireNamespace("bigsnpr", quietly = TRUE)) {
    install.packages("bigsnpr")
  }
  library(bigsnpr)

  n <- nrow(G); m <- ncol(G)

  # 2) prepare sample & SNP IDs
  sample_id <- sprintf("I%04d", seq_len(n))
  snp_id    <- paste0("rs", seq_len(m))
  CHR       <- rep(1, m)
  POS       <- seq_len(m)

  # 3) write a backing FBM (file‐backed big matrix)
  #    this only writes a small .bk + .rds file
  fbm <- FBM(n, m, backingfile = paste0(prefix, "_FBM"))
  fbm[] <- G

  # 4) use bigsnpr’s fast bed writer
  snp_writeBed(
    bedfile    = paste0(prefix, ".bed"),

  )
  # This creates prefix.bed, prefix.bim, prefix.fam in seconds.

  message("✔ BED files ready: ", prefix, ".{bed,bim,fam}")

  # 5) call PLINK 2 to convert to its pfile
  #    assumes plink2 is on your PATH
  system2("plink2", c(
    "--bfile", prefix,
    "--make-pgen",
    "--out", prefix,
    "--silent"
  ))
  message("✔ PGEN files ready: ", prefix, ".{pgen,pvar,psam}")
}
