plink2  --pfile /data5/UKB_Eighty_Million/Human_UKB_Impute  --freq --out 'matrix/allel.freq'

plink2 --pfile /data5/UKB_Eighty_Million/Human_UKB_Impute --keep /home/hjc/projects/MR/UKB_data/whiteID.id --make-pgen --out /data5/UKB_Eighty_Million/UKB_white_Impute

# test speed on one snps
plink2   --pfile /data5/UKB_Eighty_Million/UKB_white_Impute \
--out 'UKB_data/ukb_pca'   \
--pca 30 approx

# making phenotype for plink2 calculation
scaled_covar = fread('pheno/scaled/all_selected_phe_139_withcolnames_scaled.table')
write.table(scaled_covar, file='pheno/scaled/all_selected_phe_139_withcolnames_scaled.table', row.names = F, col.names = T, sep = ' ', quote = F)
# also need to change the second colname to IID
pheno = read.csv('./hearing_exam/hearing_8_pheno.csv') # first col are ids
pcs = read.table('./ukb_3_pca/ukb_pca.eigenvec')[,2:(2+total_pcs_num)] # first 1 cols are ids

# test speed on one snps
plink2   --pfile /data5/UKB_Eighty_Million/UKB_white_Impute  --glm allow-no-covars skip-invalid-pheno \
--out 'plink_temp/test_single_reg'   \
--pheno pheno/scaled/all_selected_phe_139_withcolnames_scaled.table \
--pheno-name X21001.0.0


########### test seleted SNPs approach ############
#
# require: 
# pheno_table, pheno_names

write.table(IA_clump$rsid,'./selected_snp/IA/IA_snp.table',col.names=F,row.names=F,quote=F)

### 1. plink prepare
if (!dir.exists('selected_snp/IA/1_plink_prepare')) dir.create('selected_snp/IA/1_plink_prepare')
plink2  --pfile /data5/UKB_Eighty_Million/UKB_white_Impute  \
--extract /home/hjc/projects/MR/selected_snp/IA/IA_snp.table \
--glm allow-no-covars skip-invalid-pheno \
--out 'selected_snp/IA/1_plink_prepare/1_test_single_reg'   \
--pheno hearing_exam/8pheno_20pcs_withcolnames_scaled.table 

plink2  --pfile /data5/UKB_Eighty_Million/UKB_white_Impute  \
--extract /home/hjc/projects/MR/selected_snp/IA/IA_snp.table \
--freq \
--out 'selected_snp/IA/1_plink_prepare/allele_freq'   

### 2. calulating covariance
library(glue)
if (!dir.exists('selected_snp/IA/2_covariance_matrix')) dir.create('selected_snp/IA/2_covariance_matrix')
total_pcs_num=20
pheno_names = c('X2247.0.0','X21003.0.0','X31.0.0','X6138.0.0','X2050.0.0','X1200.0.0','X48.0.0','X49.0.0',sapply(c(1:total_pcs_num),function(i){glue('PC{i}')}))
covar_num = length(pheno_names) - total_pcs_num

freq = read.table('selected_snp/IA/1_plink_prepare/allele_freq.afreq',header=T,comment.char='$')
var_x = 2*freq$ALT_FREQS*(1-freq$ALT_FREQS) # for human and other diploid, var=2pq
cov_xy = matrix(0,length(var_x),covar_num+total_pcs_num)
for (i in c(1:(covar_num+total_pcs_num))){
  print(i)
  lm = read.table(glue('./selected_snp/IA/1_plink_prepare/1_test_single_reg.{pheno_names[i]}.glm.linear'),header = T,comment.char='$')
  cov_xy_slice = var_x * lm$BETA
  cov_xy[,i] = cov_xy_slice
}
colnames(cov_xy) = pheno_names
rownames(cov_xy) = c(freq$ID)
cov_yy = cov(scaled_pheno[,-c(1,2)],use='pairwise')
colnames(cov_yy) = colnames(cov_xy)
rownames(cov_yy) = colnames(cov_xy)
# saving covariance matrix
write.table(cov_xy, file='selected_snp/IA/2_covariance_matrix/cov_x_yz.table', row.names = T, col.names = T, sep = ' ')
write.table(cov_yy, file='selected_snp/IA/2_covariance_matrix/cov_yz_yz.table', row.names = T, col.names = T, sep = ' ')
write.table(var_x, file='selected_snp/IA/2_covariance_matrix/var_x.table')

### 3. OATH algorithm
library(glue)
if (!dir.exists('selected_snp/IA/3_OATH_estimated')) dir.create('selected_snp/IA/3_OATH_estimated')
# loading covariance
var_x = read.table(glue('selected_snp/IA/2_covariance_matrix/var_x.table'))[[1]]
cov_xy = as.data.frame(as.matrix(fread(file=glue('selected_snp/IA/2_covariance_matrix/cov_x_yz.table')),rownames=1))
cov_yy = as.data.frame(as.matrix(fread(file=glue('selected_snp/IA/2_covariance_matrix/cov_yz_yz.table')),rownames=1))

p = length(var_x)
n = 270000 # est. effective sample size
total_pcs_num=20
using_pcs_num=5
phe_name = 'X6138.0.0'  # phenotype to reg
all_covar_name = c(c('X2247.0.0','X21003.0.0','X31.0.0','X6138.0.0','X2050.0.0','X1200.0.0','X48.0.0','X49.0.0')[-4],sapply(c(1:total_pcs_num),function(i){glue('PC{i}')}))
covar_num = length(all_covar_name) - total_pcs_num
all_possible_3_combn = combn(all_covar_name[1:covar_num],3)

for (j in c(1:dim(all_possible_3_combn)[2])){
  reg_summary = as.data.frame(matrix(,p,5))
  colnames(reg_summary) = c('ID','BETA','SE','T-STAT','P')
  # select one combination of covar
  covar_name = c(all_possible_3_combn[,j],sapply(c(1:using_pcs_num),function(i){glue('PC{i}')}))
  d = 2+length(covar_name)
  Theta = matrix(0,d,d)
  Theta[c(1,3:d),c(1,3:d)] = as.matrix(cov_yy[c(phe_name,covar_name),c(phe_name,covar_name)])
  cov_xy_part = cov_xy[,c(phe_name,covar_name)]
  
  start = Sys.time()
  for (i in c(1:p)){
    if (i%%10000==0) print(glue('{100*(i/p)}% calculated'))
    Theta[2,-2] = as.matrix(cov_xy_part[i,])
    Theta[-2,2] = as.matrix(t(cov_xy_part[i,]))
    Theta[2,2] = var_x[i]
    
    # get Omega, Lambda, b
    Omega = Theta[2:d,2:d]
    Lambda = diag(diag(Omega))
    b = matrix(Theta[2:d,1] / diag(Theta[2:d,2:d]),d-1,1)
    
    # compute regression parameter beta. and its variance
    rv_Omega = solve(Omega)
    beta = rv_Omega %*% Lambda %*% b
    var_beta = as.numeric((Theta[1,1]-t(beta)%*%Lambda%*%b)/(n-(d-2)-1))*rv_Omega
    
    # beta/SE follow a t-dis with df=n-2
    t_stats = beta/sqrt(diag(var_beta))
    p_val = pt(-abs(t_stats),df=n-2)*2
    
    # adding result to summary.
    reg_summary[i,2:5] = as.numeric(cbind(beta[1],sqrt(diag(var_beta)[1]),t_stats[1],p_val[1]))
  }
  
  lm = read.table(glue('./selected_snp/IA/1_plink_prepare/1_test_single_reg.{covar_name[1]}.glm.linear'),header = T,comment.char='$')
  reg_summary$ID = lm$ID
  reg_summary = data.frame(reg_summary)
  freq = read.table('selected_snp/IA/1_plink_prepare/allele_freq.afreq',header=T,comment.char='$')
  reg_summary$effect_allele = lm$A1
  reg_summary$CHR = lm$X.CHROM
  reg_summary$pos = lm$POS
  reg_summary$N = lm$OBS_CT
  reg_summary$eaf = freq$ALT_FREQS
  reg_summary$non_eff = lm$REF
  reg_summary$non_eff[lm$REF==lm$A1] = lm$ALT[lm$REF==lm$A1]
  reg_summary$non_eff[lm$ALT==lm$A1] = lm$REF[lm$ALT==lm$A1]
  
  saveRDS(list(covar_name,reg_summary),file = glue('selected_snp/IA/3_OATH_estimated/3_possible_covar_{j}.rds'))
  end = Sys.time()
  
  print(glue('Finished calculation for covariates combination {j} out of total {dim(all_possible_3_combn)[2]} combinations, using {difftime(end,start,units=\'mins\')} mins'))
}





########### apply