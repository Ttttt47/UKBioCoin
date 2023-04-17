library(glue)
library(data.table)
phes = read.table('/public3/project_users/chengb/hjc/projects/MR/pheno/')$V1
phes = paste0('X',phes,'.0.0')

PCs = paste0('PC',1:20)
paste0(c(PCs,phes),collapse = ',')
# PLINK 
subfile="/public3/project_users/chengb/hjc/projects/UKBioCoin/single_reg.pbs"

# 查找算完的表型。
library(stringr)
string = "X1190.0.0,X1200.0.0,"
pattern = "X\\d+\\.0\\.0"
str_extract_all(string, pattern)

## calculating covariance
setwd('/public3/project_users/chengb/hjc/projects/UKBioCoin')
library(glue)
library(data.table)
phes = read.table('/public3/project_users/chengb/hjc/projects/MR/pheno/pheno_set_0320_3_162')$V1
phes = paste0('X',phes,'.0.0')

PCs = paste0('PC',1:20)
gcount = fread("/public3/project_users/chengb/hjc/projects/UKBioCoin/UKB_14M_geno_count.gcount")
Ns = rowSums(gcount[,c(5,6,7)])
var_x = (4*gcount[,5]+gcount[,6])/Ns - (2*gcount[,5]/Ns+gcount[,6]/Ns)^2
var_x = var_x$HOM_REF_CT
# freq = read.table("/public3/project_users/chengb/hjc/projects/UKBioCoin/10M.allel.freq.afreq",header=T,comment.char='$')
# var_x = 2*freq$ALT_FREQS*(1-freq$ALT_FREQS) # for human and other diploid, var=2pq
cov_xy = matrix(0,length(var_x),length(PCs)+length(phes))
for (i in 1:length(c(PCs,phes))){
  print(i)
  lm = fread(glue('/public3/project_users/chengb/hjc/projects/UKBioCoin/plink_temp/single_reg.{c(PCs,phes)[i]}.glm.linear'),header = T,nThread=16)
  cov_xy_slice = var_x * lm$BETA
  cov_xy[,i] = cov_xy_slice
}

scaled_pheno = fread(file="/public3/project_users/chengb/hjc/projects/MR/pheno/valid_584pheno_20pcs_withcolnames_scaled.table")
scaled_pheno = data.frame(scaled_pheno)[,c(PCs,phes)]
colnames(cov_xy) = c(PCs,phes)

# all_pheno = cbind(scaled_phe136[,-c(1,2)],scaled_pcs[,-c(1,2)])
cov_yy = cov(scaled_pheno,use='pairwise')
colnames(cov_yy) = colnames(cov_xy)
rownames(cov_yy) = colnames(cov_xy)
# saving covariance matrix
fwrite(cov_xy[1:2,],file='/public3/project_users/chengb/hjc/projects/UKBioCoin/matrix/14M_cov_xy.table',sep=' ',na='NA',row.names = T, col.names = T, nThread=16)
# edit the header manually and run
fwrite(cov_xy,file='/public3/project_users/chengb/hjc/projects/UKBioCoin/matrix/14M_cov_xy.table',sep=' ',na='NA',row.names = T, col.names = F, append=T , nThread=16)

# fwrite(cov_yy,file='/public3/project_users/chengb/hjc/projects/UKBioCoin/matrix/14M_cov_yy.table',sep=' ',na='NA',row.names = T, col.names = T)
# write.table(cov_xy, file='./matrix/14M_cov_xy.table', row.names = T, col.names = T, sep = ' ')
write.table(cov_yy, file='/public3/project_users/chengb/hjc/projects/UKBioCoin/matrix/14M_cov_yy.table', row.names = T, col.names = T, sep = ' ')
write.table(var_x, file='./matrix/14M_var_x.table')

#### meta
afreq = fread("/public3/project_users/chengb/hjc/projects/UKBioCoin/UKB_14M_allel_freq.afreq")

colnames(afreq)[3] = 'REF_Allele'
colnames(afreq)[4] = 'ALT_Allele'
afreq$REF_FREQ = 1-afreq$ALT_FREQS
pvar = fread("/public3/project_users/chengb/DATA/ukb_white_impute_14M/UKB_white_Impute_10M_erase_dosage.pvar")
afreq$POS = pvar$POS

fwrite(afreq[,c(1,2,8,3,4,7)], file='./matrix/14M_meta.table',sep=' ',na='NA',row.names = F, col.names = T, quote=F)

#### docker


# push 
docker login -u cn-east-3@CW7MEZSWZIWWTG4ZK0TB -p \
10347289f6f305b44ba1d10b3291adbe7cf99be1a15b3a6e3afbd7ae4bbef8d4 \
swr.cn-east-3.myhuaweicloud.com

# pull
docker login -u cn-east-3@CFB4B8RZY5ZZ5RYWIHRI -p \
3f70cdbb9813cc57ad3a4b08de08217bcd5575ba09daef2b41b6388f0f336bc1 \
swr.cn-east-3.myhuaweicloud.com

docker pull swr.cn-east-3.myhuaweicloud.com/ukbiocoin/ukbcoin_10m:v1


########### apply