

## === User settings ====================================
data_file   <- "/path/to/genotype_data_prefix"   # PLINK2 pfile prefix
freq_file   <- "/path/to/allele_frequency.afreq" # PLINK2 .afreq file
pheno_file  <- "/path/to/phenotype_table.tsv"    # phenotype file
phe_names   <- c('X50.0.0')

PC_nums     <- c(0, 5, 10)   # CHANGE PC SETS HERE

ukc_matrix  <- "/path/to/ukc_matrix_prefix"      # UKC matrix prefix

out_dir     <- getwd()
plink2_bin  <- "plink2"                          # or absolute path to plink2
ukc_bin     <- "/path/to/UKBioCoin_binary"       # UKC executable

## =======================================================================

## If no frequency file is provided, generate one from the pfile
if (is.null(freq_file) || freq_file == "") {
  auto_freq_prefix <- file.path(out_dir, "auto_freq")
  message("No freq_file provided; generating allele frequency file with PLINK2...")
  cmd <- sprintf(
    '%s --pfile %s --freq --out %s',
    plink2_bin,
    data_file,
    auto_freq_prefix
  )
  system(cmd)
  freq_file <- paste0(auto_freq_prefix, ".afreq")
}

for (phe in phe_names) {

  for (PC in PC_nums) {

    ## -------------------------- Prepare covariates ----------------------
    if (PC == 0) {
      covars <- ""
      covar_opt_plink <- ""
      covar_opt_ukc   <- ""
      covar_tag       <- "0PC"
    } else {
      PCs <- paste0("PC", 1:PC, collapse = ",")
      covars <- paste0(PCs)
      covar_opt_plink <- paste0(" --covar-name ", PCs)
      covar_opt_ukc   <- paste0(" --covar ", PCs)
      covar_tag       <- paste0(PC, "PC")
    }

    ## ========= 1. PLINK2 job (16 threads) ==============================
  taskname <- paste0("bench-plink-", phe, "-", covar_tag)
  qsubfile <- paste0(taskname, ".pbs")
  logfile  <- paste0(taskname, ".log")
  errfile  <- paste0(taskname, ".err")
  timefile <- paste0(taskname, ".time.log")
  
  sink(qsubfile)
  cat("#!/bin/bash\n")
  cat(sprintf("#PBS -N %s\n", taskname))
  cat(sprintf("#PBS -o %s\n", logfile))
  cat(sprintf("#PBS -e %s\n", errfile))
  cat("#PBS -q QUEUE_NAME\n")                    # replace with your queue
  cat("#PBS -l nodes=1:ppn=16\n\n")
  cat("cd $PBS_O_WORKDIR\n")
  cat("source ~/.bashrc\n\n")
  
  cat(sprintf("/usr/bin/time -v %s \\\n", plink2_bin))
  cat("  --glm hide-covar allow-no-covars \\\n")   # <-- FIXED
  cat(sprintf("  --pfile %s \\\n", data_file))
  cat(sprintf("  --pheno %s \\\n", pheno_file))
  cat(sprintf("  --pheno-name %s \\\n", phe))
  
  ## include covariates only if PC > 0
  if (PC > 0) cat(sprintf("  %s \\\n", covar_opt_plink))
  
  cat(sprintf("  --read-freq %s \\\n", freq_file))
  cat(sprintf("  --out %s/%s \\\n", out_dir, taskname))
  cat(sprintf("  &> %s\n", timefile))
  sink()
  
  # system(paste("qsub", qsubfile))


    ## ========= 2. UKC job (1 thread) ==================================
    taskname2 <- paste0("bench-ukc-", phe, "-", covar_tag)
    qsubfile2 <- paste0(taskname2, ".pbs")
    logfile2  <- paste0(taskname2, ".log")
    errfile2  <- paste0(taskname2, ".err")
    timefile2 <- paste0(taskname2, ".time.log")

    sink(qsubfile2)
    cat("#!/bin/bash\n")
    cat(sprintf("#PBS -N %s\n", taskname2))
    cat(sprintf("#PBS -o %s\n", logfile2))
    cat(sprintf("#PBS -e %s\n", errfile2))
    cat("#PBS -q QUEUE_NAME\n")                  # replace with your queue
    cat("#PBS -l nodes=1:ppn=1\n\n")
    cat("cd $PBS_O_WORKDIR\n")
    cat("source ~/.bashrc\n\n")

    cat(sprintf("/usr/bin/time -v %s \\\n", ukc_bin))
    cat(sprintf("  --file %s \\\n", ukc_matrix))
    cat(sprintf("  --phe %s \\\n", phe))
    if (PC > 0) cat(sprintf("  %s \\\n", covar_opt_ukc))
    cat(sprintf("  --out %s/%s \\\n", out_dir, taskname2))
    cat("  --use-missing-rate-file --use-missing-rate-estimate \\\n")
    cat(sprintf("  &> %s\n", timefile2))
    sink()

    # system(paste("qsub", qsubfile2))
  }
}
