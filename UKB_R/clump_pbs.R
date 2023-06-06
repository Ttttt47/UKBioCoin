library(ggplot2)
library(reshape2)
library(circlize)
library(data.table)
library(glue)
library(cowplot)
library(ieugwasr)

setwd('/public3/project_users/chengb/hjc/projects/UKBioCoin/0525/analysis0529/adjusting_covar/')
Args <- commandArgs(trailingOnly=TRUE)
i = as.integer(Args[1])
log_sig_level = -log(5*10^{-8},10)

my_ld_clump_local <- function(dat, clump_kb, clump_r2, clump_p, bfile, plink_bin,nThreads=1)
{
  
  # Make textfile
  shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")
  fn <- tempfile()
  write.table(data.frame(SNP=dat[["rsid"]], P=dat[["pval"]]), file=fn, row.names=F, col.names=T, quote=F)
  
  fun2 <- paste0(
    shQuote(plink_bin, type=shell),
    " --bfile ", shQuote(bfile, type=shell),
    " --clump ", shQuote(fn, type=shell), 
    " --clump-p1 ", clump_p, 
    " --clump-r2 ", clump_r2, 
    " --clump-kb ", clump_kb, 
    " --threads ", nThreads, 
    " --out ", shQuote(fn, type=shell)
  )
  system(fun2)
  res <- read.table(paste(fn, ".clumped", sep=""), header=T)
  unlink(paste(fn, "*", sep=""))
  y <- subset(dat, !dat[["rsid"]] %in% res[["SNP"]])
  if(nrow(y) > 0)
  {
    message("Removing ", length(y[["rsid"]]), " of ", nrow(dat), " variants due to LD with other variants or absence from LD reference panel")
  }
  return(subset(dat, dat[["rsid"]] %in% res[["SNP"]]))
}


eval(parse(text=glue('oath_WC_res_{i}=fread(\'./oath-WC-{i}_results.table\')')))
eval(parse(text=glue('sig_WC_snps_{i} = oath_WC_res_{i}[oath_WC_res_{i}$`-log10_P`>log_sig_level,]')))
eval(parse(text=glue('da = data.frame(rsid=sig_WC_snps_{i}$ID,pval=10^(-sig_WC_snps_{i}$`-log10_P`))')))
res = my_ld_clump_local(
  dat = da,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p = 0.99,
  bfile = '/public3/project_users/chengb/DATA/ukb_white_impute_14M/UKB_white_Impute_10M_erase_dosage_excluede_palindrome_BF',
  plink_bin = '/public3/project_users/chengb/software/bin/plink',
  nThreads = 16
)

saveRDS(res,file = glue('oath_WC_{i}_clump_res.rds'))
