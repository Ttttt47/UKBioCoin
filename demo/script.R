# packages requirement
options(warn = -1)
necessary=c('data.table', 'getopt')
installed=necessary %in% installed.packages()[, 'Package']
if (length(necessary[!installed]) >=1)
  install.packages(necessary[!installed],repos='http://cran.us.r-project.org')
info=lapply(necessary, function(x){library(x,character.only = T,warn.conflicts = F,quietly = T)})
options(warn = 0)

cat_header = function(){
  version="V1.0"
  cat("###################################################\n")
  cat("## Generating NSS for UKBioCoin (UKC)\n")
  cat(paste0("## Version: ", version,"\n"))
  cat("## Written by: Jing-cheng He, Guo-An Qi, Zhejiang University\n")
  cat("## Bug report: jc_he@zju.edu.cn\n")
  cat("###################################################\n\n")
}

command = matrix(c("pfile","pf",2,"character", "PLINK2 binary genotype table, specify at least this option or --bfile as the genotype input",
                   "bfile","bf",2,"character", "PLINK1.9 binary genotype table, specify at least this option or --pfile as the genotype input",
                   "pheno","p",1,"character","ID of donor parental line",
                   "novisualize","v",2,"logical","by default, figures used to evaluate UKC results against the general results based on individual-level data is generated, specify this option to turn off this behaviour",
                   "nonormalize","n",2,"logical","by default, UKC will normalized the input phenotypes by mean 0 and sd 1, specify this option to turn off this behaviour. WARNING: Un-normalized phenotype would lead to bad NSS, this option should not be specified in most cases",
                   "threads","t",2,"numeric","Specify to set threads for PLINK, defaulted 4",
                   "memory","m",2,"numeric","Specify to set memory limits for PLINK, in MB, defaulted 8000",
                   "prefix","pr",2,"character","prefix of the generated NSS, defaulted UKC",
                   "ukc","ue",2,"character","full directory of UKC excutable, otherwise the scrip would try UKBioCoin",
                   "out", 'o',2,"character", "output directory",
                   "help","h",0,"logical", "parameters input instruction"),
                 byrow=T,ncol=5)
args = getopt(spec = command)

if (!is.null(args$help) || (is.null(args$pfile) & is.null(args$bfile)) || is.null(args$pheno)) {
  cat_header()
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}

# grab arguements
pheno = args$pheno
ukc = 'UKBioCoin'
prefix = 'UKC'
novisualize = FALSE
nonormalize = FALSE
threads = 4
memory = 8000
out = getwd()

if(is.null(args$pfile)){
  file = args$bfile
  plink_header = "plink2 --bfile "
} else {
  file = args$pfile
  plink_header = "plink2 --pfile "
}
if (!is.null(args$novisualize)) {
  novisualize = TRUE
}
if (!is.null(args$nonormalize)) {
  nonormalize = TRUE
}
if (!is.null(args$threads)) {
  threads = args$threads
}
if (!is.null(args$memory)) {
  memory = args$memory
}
if (!is.null(args$prefix)) {
  prefix = args$prefix
}
if (!is.null(args$ukc)) {
  ukc = args$ukc
}
if (!is.null(args$out)) {
  out = args$out
}
cat_header()
if(is.null(args$pfile)){
  sample_size = nrow(fread(paste0(file,".fam")))
  cat(paste0("Options: \n  --bfile ",file," \n  --pheno ",pheno," \n  --threads ",threads," \n  --memory ", memory," \n  --novisualize ",novisualize," \n  --nonormalized ", nonormalize, " \n  --prefix ", prefix, " \n  --ukc ", ukc, "\n\n"))
} else {
  sample_size = nrow(fread(paste0(file,".psam"))) - 1
  cat(paste0("Options: \n  --pfile ",file," \n  --pheno ",pheno," \n  --threads ",threads," \n  --memory ", memory," \n  --novisualize ",novisualize," \n  --nonormalized ", nonormalize, " \n  --prefix ", prefix, " \n  --ukc ", ukc, "\n\n"))
}


if(!dir.exists(out)) dir.create(out)
setwd(out)
if(!dir.exists('./0.basic')) dir.create('./0.basic')

# generating basic
TIME = 0
time1 = proc.time()
plink_command = paste0(plink_header, 
                       file,
                       ' --silent --geno-counts ',
                       ' --out ./0.basic/0.count',
                       ' --threads ', threads,
                       ' --memory ', memory)
cat(plink_command)
cat("\n")
system(plink_command)

plink_command = paste0(plink_header, 
                       file,
                       ' --silent --freq ',
                       ' --out ./0.basic/0.freq',
                       ' --threads ', threads,
                       ' --memory ', memory)
cat(plink_command)
cat("\n")
system(plink_command)

time2 = proc.time()
time = (time2-time1)/60

time = time["elapsed"]
TIME = TIME+time
cat(paste0("Generating basic statistics done. Elapsed time: ", time, " minute.\n\n"))

# generating plink temp
phe_names = paste0(colnames(fread(pheno, header = T, nrows = 1, nThread = threads))[-c(1,2)], collapse = ',')
phes = colnames(fread(pheno, header = T, nrows = 1, nThread=threads))[-c(1,2)]
pheDT = fread(pheno, header = T, data.table = F)
if(!nonormalize){
  pheDT[,phes] = apply(pheDT[,phes],2,scale)
  write.table(pheDT,paste0(pheno,".scaled"),quote=F,col.names = T,row.names = F)
  pheno = paste0(pheno,".scaled")
}

# setting up the directory
# setwd(out)
if(!dir.exists('./1.plink_temp')) dir.create('./1.plink_temp')

# run the script in command line
time1 = proc.time()
plink_command = paste0(plink_header, 
                       file,
                       ' --silent --glm allow-no-covars skip-invalid-pheno --out ./1.plink_temp/single_reg --pheno ', 
                       pheno, ' --pheno-name ', 
                       phe_names, ' --threads ', 
                       threads, ' --memory ', 
                       memory, ' --read-freq ./0.basic/0.freq.afreq')
cat(plink_command)
cat("\n")
system(plink_command)

time2 = proc.time()
time = (time2-time1)/60

time = time["elapsed"]
TIME = TIME+time
cat(paste0("Generating PLINK statistics done. Elapsed time: ", time, " minute.\n\n"))


# setwd(out)
if(!dir.exists('./2.matrix')) dir.create('./2.matrix')
time1 = proc.time()
# phenotype file
phes = colnames(fread(pheno, header = T, nrows = 1, nThread=threads))[-c(1,2)]
gcount = fread("0.basic/0.count.gcount",nThread=threads)
Ns = rowSums(gcount[,c('HOM_REF_CT','HET_REF_ALT_CTS','TWO_ALT_GENO_CTS')])
var_x = (4*gcount[,'HOM_REF_CT']+gcount[,'HET_REF_ALT_CTS'])/Ns - 
  (2*gcount[,'HOM_REF_CT']/Ns+gcount[,'HET_REF_ALT_CTS']/Ns)^2
var_x = var_x$HOM_REF_CT
cov_xy = data.table(matrix(NA,length(var_x),length(phes)))

# this is slow
for (i in 1:length(phes)){
  cat(paste0("  Proceeding ",i,"/",length(phes)," phenotypes ... \n"))
  # try to read the glm.linear file, if it is not there, then skip it.
  if (!file.exists(paste0('1.plink_temp/single_reg.', phes[i], '.glm.linear'))) next
  lm = fread(paste0('1.plink_temp/single_reg.', phes[i], '.glm.linear'), header = T, nThread = threads)
  cov_xy_slice = var_x * lm$BETA
  cov_xy[,i] = cov_xy_slice
}

scaled_pheno = fread(pheno)
scaled_pheno = data.frame(scaled_pheno)[,phes]


if(!nonormalize){
  scaled_pheno = apply(scaled_pheno,2,scale)
}

colnames(cov_xy) = phes
cov_yy = cov(scaled_pheno,use='pairwise')
colnames(cov_yy) = colnames(cov_xy)
rownames(cov_yy) = colnames(cov_xy)


# write the first row of cov_xy to a file
write.table(cov_xy[1:2,], file=paste0('2.matrix/', prefix, '_cov_xy.table'), row.names = T, col.names = T, sep = ' ')
# write the rest of cov_xy to a file using fwrite()
fwrite(cov_xy[3:dim(cov_xy)[1],], file=paste0('2.matrix/', prefix, '_cov_xy.table'), sep=' ',na='NA',row.names = T, col.names = F, quote=T, append=T)
write.table(cov_yy, file=paste0('2.matrix/', prefix, '_cov_yy.table'), row.names = T, col.names = T, sep = ' ')
write.table(var_x, file=paste0('2.matrix/', prefix, '_var_x.table'))


## meta, you may use other meta info. other than these.
afreq = fread("0.basic/0.freq.afreq", nThread = threads)
colnames(afreq)[3] = 'REF_Allele'
colnames(afreq)[4] = 'ALT_Allele'
afreq$REF_FREQ = 1-afreq$ALT_FREQS
if(is.null(args$bfile)){
  pvar = fread(paste0(file,".pvar"), nThread = threads)
} else {
  pvar = fread(paste0(file,".bim"), nThread = threads)
  colnames(pvar) = c('#CHROM', 'ID', 'GENETICPOS','POS',"ALT","REF")
}

afreq$POS = pvar$POS
afreq$Effect_Allele = ifelse(afreq$REF_FREQ > 0.5, afreq$ALT_Allele, afreq$REF_Allele)
fwrite(afreq[,c('#CHROM', 'ID', 'POS', 'Effect_Allele', 'REF_Allele',
                'ALT_Allele', 'REF_FREQ')], file=paste0('2.matrix/', prefix, '_meta.table'),sep=' ',na='NA',row.names = F, col.names = T, quote=F)


## missing rates of y and x, these are needed to estimate the sample size
y_missing_rate = sapply((scaled_pheno[phes]),FUN=function(x) mean(is.na(x)))
write.table(y_missing_rate, row.names = T, col.names = 'missing_rate',
            file=paste0('2.matrix/', prefix, '_y_missing.table'))
x_missing = gcount$MISSING_CT /(gcount$MISSING_CT + Ns)
write.table(x_missing, row.names = T, col.names = 'x_missing',
            file=paste0('2.matrix/', prefix, '_x_missing.table'))

time2 = proc.time()
time = (time2-time1)/60

time = time["elapsed"]
TIME = TIME+time
cat(paste0("Generating NSS done. Elapsed time: ", time, " minute.\n\n"))


if(!novisualize){
  # setwd(out)
  if(!dir.exists('./3.analysis')) dir.create('./3.analysis')
  time1 = proc.time()
  # test UKC
  UKC_command = paste0(paste0(ukc,' --file 2.matrix/', prefix, ' --phe ', phes[1], ' --use-missing-rate-file --use-missing-rate-estimate --totalsize ', sample_size, ' --out 3.analysis/test.', phes[1]))
  cat(UKC_command)
  system(UKC_command)
  
  res.file = paste0("3.analysis/test.",phes[1],"_results.table")
  dt1 = fread(res.file,nThread = threads, header = T, data.table = F)
  dt2 = fread(paste0("1.plink_temp/single_reg.",phes[1],".glm.linear"),nThread = threads, header = T, data.table = F)
  
  ambiguous_ids = unique(c(which(dt1$REF_FREQ==0.5), which(dt2$A1_FREQ==0.5)))
  if (length(ambiguous_ids) > 0){
    cat(paste0("Warning: ", length(ambiguous_ids), " ambiguous SNPs found, removing them.\n"))
    dt1 = dt1[-ambiguous_ids,]
    dt2 = dt2[-ambiguous_ids,]
  }

  #dt1$BETA = ifelse(dt1$Effect_Allele==dt2$A1, dt1$BETA, ifelse(dt1$REF_FREQ==0.5, dt1$BETA, -1*dt1$BETA))
  dt1$BETA = ifelse(dt1$Effect_Allele==dt2$A1, dt1$BETA, -1*dt1$BETA)
  dt1 = dt1[,c('ID','BETA','SE','T-STAT','-log10_P')]
  colnames(dt1) = c('SNP','UKC_BETA','UKC_SE','UKC_TSTAT','UKC_P')
  dt1$UKC_P = 10^(-1*dt1$UKC_P)
  
  dt2 = dt2[,c('ID','BETA','SE','T_STAT','P')]
  colnames(dt2) = c('SNP','UKB_BETA','UKB_SE','UKB_TSTAT','UKB_P')
  
  dt = merge(dt1,dt2,'SNP')
  png(paste0("3.analysis/Validation.", phes[1],".png"),width = 16*300, height = 16*300,res = 300)
  par(mfrow=c(2,2))
  plot(dt$UKB_BETA, dt$UKC_BETA, xlab = "UKB_BETA", ylab = "UKC_BETA", 
       xlim = 1.05*c(min(dt$UKB_BETA),max(dt$UKB_BETA)),
       ylim = 1.05*c(min(dt$UKC_BETA),max(dt$UKC_BETA)))
  abline(a = 0, b = 1, col = 'red')
  usr = par("usr")
  xpos = usr[1] + 0.2 * (usr[2] - usr[1])
  ypos = usr[3] + 0.8 * (usr[4] - usr[3])
  text(xpos, ypos, labels = paste0("Cor = ",round(cor(dt$UKB_BETA, dt$UKC_BETA),digits = 2)))
  
  plot(dt$UKB_SE, dt$UKC_SE, xlab = "UKB_SE", ylab = "UKC_SE", 
       xlim = 1.05*c(min(dt$UKB_SE),max(dt$UKB_SE)),
       ylim = 1.05*c(min(dt$UKC_SE),max(dt$UKC_SE)))
  abline(a = 0, b = 1, col = 'red')
  usr = par("usr")
  xpos = usr[1] + 0.2 * (usr[2] - usr[1])
  ypos = usr[3] + 0.8 * (usr[4] - usr[3])
  text(xpos, ypos, labels = paste0("Cor = ",round(cor(dt$UKB_SE, dt$UKC_SE),digits = 2)))
  
  plot(dt$UKB_TSTAT, dt$UKC_TSTAT, xlab = "UKB_TSTAT", ylab = "UKC_TSTAT", 
       xlim = 1.05*c(min(dt$UKB_TSTAT),max(dt$UKB_TSTAT)),
       ylim = 1.05*c(min(dt$UKC_TSTAT),max(dt$UKC_TSTAT)))
  abline(a = 0, b = 1, col = 'red')
  usr = par("usr")
  xpos = usr[1] + 0.2 * (usr[2] - usr[1])
  ypos = usr[3] + 0.8 * (usr[4] - usr[3])
  text(xpos, ypos, labels = paste0("Cor = ",round(cor(dt$UKB_TSTAT, dt$UKC_TSTAT),digits = 2)))
  
  plot(dt$UKB_P, dt$UKC_P, xlab = "UKB_P", ylab = "UKC_P", 
       xlim = 1.05*c(min(dt$UKB_P),max(dt$UKB_P)),
       ylim = 1.05*c(min(dt$UKC_P),max(dt$UKC_P)))
  abline(a = 0, b = 1, col = 'red')
  usr = par("usr")
  xpos = usr[1] + 0.2 * (usr[2] - usr[1])
  ypos = usr[3] + 0.8 * (usr[4] - usr[3])
  text(xpos, ypos, labels = paste0("Cor = ",round(cor(dt$UKB_P, dt$UKC_P),digits = 2)))
  
  dev.off()
  
  time2 = proc.time()
  time = (time2-time1)/60
  
  time = time["elapsed"]
  TIME = TIME+time
  cat(paste0("Plotting validated figures done. Elapsed time: ", time, " minute.\n\n"))
}
cat(paste0("ALL UKC NSS generation done. Elapsed time: ", TIME, " minute.\n\n"))
