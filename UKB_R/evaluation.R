library(ggplot2)
library(reshape2)
library(circlize)
library(data.table)
library(glue)
library(cowplot)

setwd('/public3/project_users/chengb/hjc/projects/UKBioCoin/0422/')

# Allel freq plot ---------------------------------------------------------
freq = fread("/public3/project_users/chengb/hjc/projects/UKBioCoin/0422/0.basic/0.freq.afreq")
MAF = sapply(freq$ALT_FREQS,function(x) min(x,1-x))

pdf('./analysis0523/basics/allele_freq.pdf',height = 6, width = 7)
hist(MAF,breaks=100)
dev.off()

#  near complete -----------------------------------------------------------
setwd('/public3/project_users/chengb/hjc/projects/UKBioCoin/0525/analysis0529/reproducing_GWAS/nearcom/')
phe = 'X20015.0.0'

color_plates = c("#0071C2","#D75615","#EDB11A")
for (PC_num in c(0,5,10)) {
  eval(parse(text=glue('plink_res_{PC_num}pc=fread(\'./plink-{phe}-{PC_num}PC.{phe}.glm.linear\')')))
  eval(parse(text=glue('plink_res_{PC_num}pc=plink_res_{PC_num}pc[plink_res_{PC_num}pc$TEST==\'ADD\',]')))
  eval(parse(text=glue('oath_res_{PC_num}pc=fread(\'./oath-{phe}-{PC_num}PC_results.table\')')))
  print(PC_num)
}

# Comparing PLINK results and OATH results
options(digits = 5)
plot_list = list()
for (i in c(1:3)) {
  PC_num = c(0,5,10)[i]
  # cor = eval(parse(text = glue('cor(plink_res_{PC_num}pc$BETA,oath_res_{PC_num}pc$BETA)')))
  # cor = round(cor,digits = 3)
  # ids = sample(1:dim(oath_res_0pc)[1],100000)
  ids = 1:dim(oath_res_0pc)[1]
  plink_points = eval(parse(text = glue('plink_res_{PC_num}pc$BETA[ids]')))
  oath_points = eval(parse(text = glue('oath_res_{PC_num}pc$BETA[ids]')))
  p = ggplot(data.frame(x=plink_points,y=oath_points), aes(x = x, y = y)) +
    scale_x_continuous(name = 'PLINK estimates',limits = c(-0.5,0.5)) + 
    scale_y_continuous(name = 'UKC estimates',limits = c(-0.5,0.5)) + 
    geom_point(color = color_plates[i], size=0.3, alpha = 0.2) +
    ggtitle(glue('Beta, {PC_num} PC', )) +
    # annotate("text", x = -1.5, y = 1.5, label = glue("{cor}"), size = 4) +
    geom_abline(intercept=0, slope=1, linewidth= 0.4, color='grey', alpha=0.7)
  # theme_minimal()
  plot_list[[length(plot_list)+1]] = p
}

for (i in c(1:3)) {
  PC_num = c(0,5,10)[i]
  plink_points = eval(parse(text = glue('-log(plink_res_{PC_num}pc$P[ids],10)')))
  oath_points = eval(parse(text = glue('oath_res_{PC_num}pc$`-log10_P`[ids]')))
  p = ggplot(data.frame(x=plink_points,y=oath_points), aes(x = x, y = y)) +
    scale_x_continuous(name = 'PLINK estimates',limits = c(0,100)) + 
    scale_y_continuous(name = 'UKC estimates',limits = c(0,100)) + 
    geom_point(color = color_plates[i], size=0.3, alpha = 0.2) +
    ggtitle(glue('-log10(p), {PC_num} PC')) +
    # annotate("text", x = -1.5, y = 1.5, label = glue("{cor}"), size = 4) +
    geom_abline(intercept=0, slope=1, linewidth= 0.4, color='grey', alpha=0.7)
  # theme_minimal()
  plot_list[[length(plot_list)+1]] = p
}

jpeg('./trueBeta&pvalues.vs.oath.jpg', height = 1800, width = 900*3, quality = 100, res = 300)
plot_grid(plotlist = plot_list,nrow = 2)
dev.off()



#  high missing -----------------------------------------------------------
setwd('/public3/project_users/chengb/hjc/projects/UKBioCoin/0525/analysis0529/reproducing_GWAS/highmissing/')
phe = 'X20127.0.0'
color_plates = c("#0071C2","#D75615","#EDB11A","#369F2D")

for (PC_num in c(5)) {
  for (i in 0:3) {
    eval(parse(text=glue('plink_res_{PC_num}pc_{i}=fread(\'./plink-{phe}-{PC_num}PC-{i}.{phe}.glm.linear\')')))
    eval(parse(text=glue('plink_res_{PC_num}pc_{i}=plink_res_{PC_num}pc_{i}[plink_res_{PC_num}pc_{i}$TEST==\'ADD\',]')))
    eval(parse(text=glue('oath_res_{PC_num}pc_{i}=fread(\'./oath-{phe}-{PC_num}PC-{i}_results.table\')')))
    print(i)
  }
}

# Comparing PLINK results and OATH results

# ids = sample(1:dim(oath_res_5pc_1)[1],10000)
ids = 1:dim(oath_res_5pc_1)[1]


covars = c('None','ID.1269','ID.1269, ID.1210','ID.1269, ID.1210, ID.1618')
options(digits = 5)
plot_list = list()
for (i in c(0:3)) {
  PC_num = 5
  # cor = eval(parse(text = glue('cor(plink_res_{PC_num}pc$BETA,oath_res_{PC_num}pc$BETA)')))
  # cor = round(cor,digits = 3)
  
  plink_points = eval(parse(text = glue('plink_res_{PC_num}pc_{i}$BETA[ids]')))
  oath_points = eval(parse(text = glue('oath_res_{PC_num}pc_{i}$BETA[ids]')))
  p = ggplot(data.frame(x=plink_points,y=oath_points), aes(x = x, y = y)) +
    scale_x_continuous(name = 'PLINK estimates',limits = c(-0.5,0.5)) + 
    scale_y_continuous(name = 'UKC estimates',limits = c(-0.5,0.5)) + 
    geom_point(color = color_plates[i+1], size=0.3, alpha = 0.2) +
    ggtitle(glue('Beta, {PC_num} PC & {i} covarates adjusted')) +
    # annotate("text", x = -1.5, y = 1.5, label = glue("{cor}"), size = 4) +
    geom_abline(intercept=0, slope=1, linewidth= 0.4, color='grey', alpha=0.7)+
    theme(plot.title = element_text(size = 10))
  # theme_minimal()
  plot_list[[length(plot_list)+1]] = p
}

for (i in c(0:3)) {
  PC_num = 5
  plink_points = eval(parse(text = glue('-log(plink_res_{PC_num}pc_{i}$P[ids],10)')))
  oath_points = eval(parse(text = glue('oath_res_{PC_num}pc_{i}$`-log10_P`[ids]')))
  p = ggplot(data.frame(x=plink_points,y=oath_points), aes(x = x, y = y)) +
    scale_x_continuous(name = 'PLINK estimates',limits = c(0,25)) + 
    scale_y_continuous(name = 'UKC estimates',limits = c(0,25)) + 
    geom_point(color = color_plates[i+1], size=0.3, alpha = 0.2) +
    ggtitle(glue('-log10(p), {PC_num} PC & {i} covarates adjusted')) +
    # annotate("text", x = -1.5, y = 1.5, label = glue("{cor}"), size = 4) +
    geom_abline(intercept=0, slope=1, linewidth= 0.4, color='grey', alpha=0.7)+
    theme(plot.title = element_text(size = 10))
  # theme_minimal()
  plot_list[[length(plot_list)+1]] = p
}

jpeg('./X20127_0530_trueBeta&pvalues.vs.oath.jpg', height = 1800, width = 900*4, quality = 100, res = 300)
plot_grid(plotlist = plot_list,nrow = 2)
dev.off()





# table -------------------------------------------------------------------
cut = -log10(5*10^(-8))

df_all = data.frame()
for (i in c(1:3)) {
  PC_num = c(0,5,10)[i]
  eval(parse(text=glue('B=plink_res_{PC_num}pc[-log(plink_res_{PC_num}pc$P,10)>cut,]')))
  eval(parse(text=glue('A=oath_res_{PC_num}pc[oath_res_{PC_num}pc$`-log10_P`>cut,]')))
  both = merge(A,by.x='ID',B,by.y='ID')
  df_all = rbind(df_all,c(i,
                          PC_num,
                          length(setdiff(A$ID,B$ID)),
                          length(intersect(A$ID,B$ID)),
                          length(setdiff(B$ID,A$ID)),
                          length(intersect(A$ID,B$ID))/length(union(A$ID,B$ID)),
                          sqrt(mean((both$BETA.x-both$BETA.y)^2))))
}

for (i in c(0:3)) {
  PC_num = 5
  eval(parse(text=glue('B=plink_res_{PC_num}pc_{i}[-log(plink_res_{PC_num}pc_{i}$P,10)>cut,]')))
  eval(parse(text=glue('A=oath_res_{PC_num}pc_{i}[oath_res_{PC_num}pc_{i}$`-log10_P`>cut,]')))
  both = merge(A,by.x='ID',B,by.y='ID')
  df_all = rbind(df_all,c(i,
                          PC_num,
                          length(setdiff(A$ID,B$ID)),
                          length(intersect(A$ID,B$ID)),
                          length(setdiff(B$ID,A$ID)),
                          length(intersect(A$ID,B$ID))/length(union(A$ID,B$ID)),
                          sqrt(mean((both$BETA.x-both$BETA.y)^2))))
}

colnames(df_all) = c('i','PC','A\\B','AUB','B\\A','jacard','r.m.s.e of BETA')



for (i in c(0:3)) {
  PC_num = 5
  eval(parse(text=glue('print(paste0(\'plink-oath: \', length(setdiff(plink_res_{PC_num}pc_cut$ID,oath_res_{PC_num}pc_cut$ID))))')))
  eval(parse(text=glue('print(paste0(\'oath-plink: \', length(setdiff(oath_res_{PC_num}pc_cut$ID,plink_res_{PC_num}pc_cut$ID))))')))
  eval(parse(text=glue('print(length(intersect(plink_res_{PC_num}pc_cut$ID,oath_res_{PC_num}pc_cut$ID)))')))
}


# cool --------------------------------------------------------------------




