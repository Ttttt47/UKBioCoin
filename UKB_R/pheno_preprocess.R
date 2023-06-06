library(data.table)
pheno = fread("/public3/project_users/chengb/hjc/projects/MR/pheno/584pheno_20pcs_withcolnames_nonscaled.table",data.table=F)
meta = fread("/public3/project_users/chengb/hjc/projects/MR/pheno/description of phenotype0421.csv")
meta = meta[,c(3,8,20)]
meta = meta[ValueType!='Categorical multiple'&V20!='multiple']


meta_c = meta[ValueType=='Continuous']
meta_i = meta[ValueType=='Integer']
meta_s = meta[ValueType=='Categorical single']

newpheno = pheno

newpheno[,paste0('X',meta_i[V20=='Data-Coding 100291']$FieldID,'.0.0')][
  pheno[,paste0('X',meta_i[V20=='Data-Coding 100291']$FieldID,'.0.0')]<0] = NA

# less than one set as NA.
newpheno[,paste0('X',meta_i[V20=='Data-Coding 100373']$FieldID,'.0.0')][
  pheno[,paste0('X',meta_i[V20=='Data-Coding 100373']$FieldID,'.0.0')]<0] = NA

# unabole to walk set to NA.
newpheno[,paste0('X',meta_i[V20=='Data-Coding 100307']$FieldID,'.0.0')][
  pheno[,paste0('X',meta_i[V20=='Data-Coding 100307']$FieldID,'.0.0')]<0] = NA

# Less than an hour a day set to 0.5.
newpheno[,paste0('X',meta_i[V20=='Data-Coding 100329']$FieldID,'.0.0')][
  pheno[,paste0('X',meta_i[V20=='Data-Coding 100329']$FieldID,'.0.0')]==-10] = 0.5

newpheno[,paste0('X',meta_i[V20=='Data-Coding 100329']$FieldID,'.0.0')][
  pheno[,paste0('X',meta_i[V20=='Data-Coding 100329']$FieldID,'.0.0')] %in% c(-1,-3)] = NA


newpheno[,paste0('X',meta_i[V20=='Data-Coding 100504']$FieldID,'.0.0')][
  pheno[,paste0('X',meta_i[V20=='Data-Coding 100504']$FieldID,'.0.0')]<0] = NA

# Less than once a year set to NA.
newpheno[,paste0('X',meta_i[V20=='Data-Coding 100537']$FieldID,'.0.0')][
  pheno[,paste0('X',meta_i[V20=='Data-Coding 100537']$FieldID,'.0.0')]<0] = NA

# Less than a year set to NA.
newpheno[,paste0('X',meta_i[V20=='Data-Coding 100290']$FieldID,'.0.0')][
  pheno[,paste0('X',meta_i[V20=='Data-Coding 100290']$FieldID,'.0.0')]<0] = NA

newpheno[,paste0('X',meta_s[V20=='ordinal']$FieldID,'.0.0')][
  pheno[,paste0('X',meta_s[V20=='ordinal']$FieldID,'.0.0')] < 0] = NA

# Prefer not to answer set to NA
newpheno[,paste0('X',meta_s[V20=='Data-Coding 100349']$FieldID,'.0.0')][
  pheno[,paste0('X',meta_s[V20=='Data-Coding 100349']$FieldID,'.0.0')]<0] = NA

# Prefer not to answer set to NA
newpheno[,paste0('X',meta_s[V20=='Data-Coding 100349.']$FieldID,'.0.0')][
  pheno[,paste0('X',meta_s[V20=='Data-Coding 100349.']$FieldID,'.0.0')]<0] = NA

# Is not ordinal but and change it to, but we prefer just delete it.
meta_s = meta_s[V20!='Data-Coding 100428']
meta_s = meta_s[V20!='Data-Coding 100429']

# Prefer not to answer set to NA
newpheno[,paste0('X',meta_s[V20=='Data-Coding 100352']$FieldID,'.0.0')][
  pheno[,paste0('X',meta_s[V20=='Data-Coding 100352']$FieldID,'.0.0')]<0] = NA

# Prefer not to answer set to NA
newpheno[,paste0('X',meta_s[V20=='Data-Coding 100508']$FieldID,'.0.0')][
  pheno[,paste0('X',meta_s[V20=='Data-Coding 100508']$FieldID,'.0.0')]<0] = NA

# 'I am completely deaf' set to NA
newpheno[,paste0('X',meta_s[V20=='Data-Coding 100631']$FieldID,'.0.0')][
  pheno[,paste0('X',meta_s[V20=='Data-Coding 100631']$FieldID,'.0.0')]==99] = NA
newpheno[,paste0('X',meta_s[V20=='Data-Coding 100631']$FieldID,'.0.0')][
  pheno[,paste0('X',meta_s[V20=='Data-Coding 100631']$FieldID,'.0.0')]<0] = NA

# Prefer not to answer set to NA
newpheno[,paste0('X',meta_s[V20=='Data-Coding 100603']$FieldID,'.0.0')][
  pheno[,paste0('X',meta_s[V20=='Data-Coding 100603']$FieldID,'.0.0')]<0] = NA

# Prefer not to answer set to NA
newpheno[,paste0('X',meta_s[V20=='Data-Coding 100417']$FieldID,'.0.0')][
  pheno[,paste0('X',meta_s[V20=='Data-Coding 100417']$FieldID,'.0.0')]<0] = NA

# "It variesr set to" 0.5
newpheno[,paste0('X',meta_s[V20=='Data-Coding 100416']$FieldID,'.0.0')][
  pheno[,paste0('X',meta_s[V20=='Data-Coding 100416']$FieldID,'.0.0')]<0] = NA
newpheno[,paste0('X',meta_s[V20=='Data-Coding 100416']$FieldID,'.0.0')][
  pheno[,paste0('X',meta_s[V20=='Data-Coding 100416']$FieldID,'.0.0')]==-6] = 0.5


# Is not ordinal but and change it to, but we prefer just delete it.
meta_s = meta_s[V20!='Data-Coding 100347']

# Prefer not to answer set to NA
newpheno[,paste0('X',meta_s[V20=='Data-Coding 90']$FieldID,'.0.0')][
  pheno[,paste0('X',meta_s[V20=='Data-Coding 90']$FieldID,'.0.0')]<0] = NA

# Prefer not to answer set to NA
newpheno[,paste0('X',meta_s[V20=='Data-Coding 100370']$FieldID,'.0.0')][
  pheno[,paste0('X',meta_s[V20=='Data-Coding 100370']$FieldID,'.0.0')]<0] = NA

# Prefer not to answer set to NA
newpheno[,paste0('X',meta_s[V20=='Data-Coding 100370']$FieldID,'.0.0')][
  pheno[,paste0('X',meta_s[V20=='Data-Coding 100370']$FieldID,'.0.0')]<0] = NA

# Is not ordinal but and change it to, but we prefer just delete it.
meta_s = meta_s[V20!='Data-Coding 100435']

new_meta = rbind(meta_c,meta_i,meta_s)
new_meta$missing_rate = sapply((newpheno[,paste0('X',new_meta$FieldID,'.0.0')]),FUN=function(x) mean(!is.na(x)))
fwrite(new_meta,"/public3/project_users/chengb/hjc/projects/MR/pheno/pheno_set_0525_129.table",quote=F,sep=' ')
newpheno = cbind(newpheno[,c(1,2)],newpheno[,paste0('X',new_meta$FieldID,'.0.0')],pheno[,paste0('PC',1:20)])
fwrite(newpheno,"/public3/project_users/chengb/hjc/projects/MR/pheno/129_filtered_pheno_20pcs_non_scaled_0525.table",na='NA',quote=F,sep=' ',row.names=F)
scaledpheno = cbind(newpheno[,c(1,2)],scale(newpheno[,paste0('X',new_meta$FieldID,'.0.0')]),scale(pheno[,paste0('PC',1:20)]))
fwrite(scaledpheno,"/public3/project_users/chengb/hjc/projects/MR/pheno/129_filtered_pheno_20pcs_scaled_0525.table",na='NA',quote=F,sep=' ',row.names=F)


meta = fread("/public3/project_users/chengb/hjc/projects/MR/pheno/description of phenotype0421.csv")
meta = meta[FieldID %in% new_meta$FieldID,-c(16:19)]
meta = merge(x=meta,by.x='FieldID',y=new_meta[,c(1,4)],by.y='FieldID')

library(openxlsx)
colnames(meta)[16] = 'type'
write.xlsx(meta,"/public3/project_users/chengb/hjc/projects/MR/pheno/description_of_phenotype_0525.xlsx")
