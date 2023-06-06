
### input:
# {sample}_geno_count.gcount
# plink_temp/single_reg.xxx.glm.linear # all single variates regression results
# pheno_set # phenotype names.
# {sample}.phe
# {sample}.afreq

## calculating covariance
library(glue)
library(data.table)
phes = read.table('pheno_set')$V1
PCs = paste0('PC',1:20)

gcount = fread("0.basic/{sample}_geno_count.gcount")
N = 292216 # modified it with your sample size
var_x = (4*gcount[,5] + gcount[,6])/(N-gcount[,10]) - ((2*gcount[,5] + gcount[,6])/(N-gcount[,10]))^2
var_x = var_x[[1]]
cov_xy = matrix(0,length(var_x),length(PCs)+length(phes))
for (i in 1:length(c(PCs,phes))){
  print(i)
  lm = fread(glue('1.plink_temp/single_reg.{c(PCs,phes)[i]}.glm.linear'),header = T,nThread=10)
  cov_xy_slice = var_x * lm$BETA
  cov_xy[,i] = cov_xy_slice
}

scaled_pheno = fread(file="{sample}.phe")
scaled_pheno = data.frame(scaled_pheno)[,c(PCs,phes)]
colnames(cov_xy) = c(PCs,phes)

cov_yy = cov(scaled_pheno,use='pairwise')
colnames(cov_yy) = colnames(cov_xy)
rownames(cov_yy) = colnames(cov_xy)

# this could be slow and take hours.
write.table(cov_xy, file=glue('./matrix/{sample}_cov_xy.table'), row.names = T, col.names = T, sep = ' ')
write.table(cov_yy, file=glue('matrix/{sample}_cov_yy.table'), row.names = T, col.names = T, sep = ' ')
write.table(var_x, file=glue('./matrix/{sample}_var_x.table'))

# the './matrix/{sample}_meta.table' file could in tailored with your interests, 
# for example, you can produce it from a plink generated frequency file

afreq = fread("{sample}.afreq")
colnames(afreq)[3] = 'REF_Allele'
colnames(afreq)[4] = 'ALT_Allele'
afreq$REF_FREQ = 1-afreq$ALT_FREQS
pvar = fread("{sample}.pvar")
afreq$POS = pvar$POS
fwrite(afreq[,c(1,2,8,3,4,7)], file='./matrix/{sample}_meta.table',sep=' ',na='NA',row.names = F, col.names = T, quote=F)

write.table(meta$missing_rate, row.names = c(paste0('X',meta$FieldID,'.0.0')), file=glue('/public3/project_users/chengb/hjc/projects/UKBioCoin/0422/2.matrix/10M_y_missing.table'))
write.table(x_missing, file=glue('/public3/project_users/chengb/hjc/projects/UKBioCoin/0422/2.matrix/10M_x_missing.table'))


### output:
# {sample}_cov_yy.table
# {sample}_cov_xy.table
# {sample}_var_x.table
# {sample}_meta.table
# 
# all the needed NSS for BioCoin algorithm is now ready.
