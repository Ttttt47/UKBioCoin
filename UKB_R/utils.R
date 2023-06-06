
ms = data.frame(missing=meta$missing_rate, col.names = c(paste0('X',meta$FieldID,'.0.0')))
write.table(meta$missing_rate, row.names = c(paste0('X',meta$FieldID,'.0.0')), file=glue('/public3/project_users/chengb/hjc/projects/UKBioCoin/0422/2.matrix/10M_y_missing.table'))
write.table(x_missing, file=glue('/public3/project_users/chengb/hjc/projects/UKBioCoin/0422/2.matrix/10M_x_missing.table'))

# create sam for testing
library(data.table)
library(glue)

dir.create('./sample')
for (m in c('cov_xy','meta','var_x','x_missing')){
  print(m)
  infcon <- file(glue('./10M_{m}.table'), open="rt") 
  outfcon <- file(glue('./sample/sam_10M_{m}.table'), open="wt") 
  batch <- 1000

  lines <- readLines(infcon, n=batch) 
  writeLines(lines, con=outfcon) 

  close(outfcon) 
  close(infcon)
}

for (m in c('cov_yy','y_missing')){
  print(m)
  system(glue('cp 10M_{m}.table ./sample/sam_10M_{m}.table'))
}
