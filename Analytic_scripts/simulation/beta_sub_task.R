

source("functions.R")

args = commandArgs(trailingOnly = TRUE)
trait = as.character(args[1])
sig_maf = as.numeric(args[2])
h2 = as.numeric(args[3])
mcar = as.numeric(args[4])
B_start = as.numeric(args[5])
B_end = as.numeric(args[6])
n = as.numeric(args[7])
prefix = as.character(args[8])

# test parameters
B = B_end - B_start + 1

m = 10000
if (!dir.exists(prefix)) {
  dir.create(prefix)
}
setwd(prefix)

out_pref = paste0(trait, "_", sig_maf, "_", h2, "_", mcar, "_n", n)
if (!dir.exists(out_pref)) {
  dir.create(out_pref)
}
setwd(out_pref)
time1 = Sys.time()
for (b in B_start:B_end){

    out_pref_b = paste0(trait, "_", sig_maf, "_", h2, "_", mcar, "_", b)
    sim <- generate_GxCyY(n=n,
                        m=m,
                        s_y=10,
                        s_c2=20,
                        s_c3=200,
                        maf_sig=sig_maf,
                        maf_bg=c(0.01,0.50),
                        h2=h2,
                        R2_c2=0.3,
                        R2_c3=0.3,
                        gamma=c(0.3,0.3,0.3),
                        trait=trait,
                        prevalence=0.20,
                        mcar=mcar)

    # write PLINK1 bed files

    write_plink(out_pref_b,t(sim$G))
    # phenotype + covariate table
    ids <- data.table(FID = 1:nrow(sim$G),
                    IID = 1:nrow(sim$G))
    phc <- cbind(ids,
                PHENO = sim$Y + 1,
                sim$covariates)
    scaled_pheno = cbind(phc[,c(1,2)],scale(phc[,-c(1,2)]))
    fwrite(scaled_pheno, file=paste0(out_pref_b,".pheno_cov.txt"),na='NA',quote=F,sep=' ',row.names=F)
    fwrite(phc, file=paste0(out_pref_b,".pheno_cov.raw"),na='NA',quote=F,sep=' ',row.names=F)
    # run GWAS: correct model
    if (trait == "binary") {
        plink_pheno = paste0(out_pref_b,".pheno_cov.raw")
    } else {
        plink_pheno = paste0(out_pref_b,".pheno_cov.txt")
    }
    system2("plink2", args=c("--bfile", out_pref_b,
                            "--pheno", plink_pheno,
                            "--pheno-name","PHENO",
                            "--covar", plink_pheno,
                            "--covar-name","C1,C2,C3",
                            "--glm hide-covar","--out", paste0(out_pref_b,"_correct"),
                            "--silent"))

    # run GWAS: wrong model
    # system2("plink2", args=c("--bfile", out_pref_b,
    #                         "--pheno", plink_pheno,
    #                         "--pheno-name","PHENO",
    #                         "--glm allow-no-covars","--out", paste0(out_pref_b,"_wrong"),
    #                         "--silent"))

    # run UKC procedure
    # generate NSS
    NSS_dir = paste0(out_pref_b,"_NSS")
    if (!dir.exists(NSS_dir)) {
        dir.create(NSS_dir)
    }
    script_dir = "/data/jingcheng/project/UKBioCoin/simulation/script.R"
    system(paste("Rscript",script_dir, "--bfile", out_pref_b,
                            "--pheno", paste0(out_pref_b,".pheno_cov.txt"),
                            "--novisualize","T",
                            "--threads", "1",
                            "--memory", "1000",
                            "--out", NSS_dir,
                            "> /dev/null 2>&1"))
    # run UKC
    system(paste("UKBioCoin", "--file", paste0(NSS_dir,"/2.matrix/UKC"),
                            "--phe", "PHENO",
                            "--covar", "C1,C2,C3",
                            "--use-missing-rate-estimate --use-missing-rate-file",
                            "--totalsize", sprintf("%d", n),
                            "--out", paste0(out_pref_b,"_UKC"),
                            "> /dev/null 2>&1"))
    sim$G = NULL
    saveRDS(sim, file=paste0(out_pref_b,".rds"))
    # delete intermediate files
    system(paste("rm", paste0(out_pref_b,"_NSS -rf")))
    system(paste("rm", paste0(out_pref_b,".bim")))
    system(paste("rm", paste0(out_pref_b,".fam")))
    system(paste("rm", paste0(out_pref_b,".bed")))
    system(paste("rm", paste0(out_pref_b,".pheno_cov.txt")))
    time2 = Sys.time()
    message("Done: ", b, " / ", B, " in ", round(difftime(time2,time1,units="mins"),2), " mins")
    time1 = time2
}