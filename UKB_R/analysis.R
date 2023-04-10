library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)

cov_yy = as.matrix(read.table('./14M_cov_yy.table'))
cov_yy_df = as.data.frame(cov_yy)

# 绘制热图
ggplot(cov_yy_df, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank())
  
heatmap_object = Heatmap(cov_yy, col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
        , cluster_rows = TRUE
        , cluster_columns = TRUE,
        show_column_names =  F,
        show_row_names = F,name = 'Correlation of\n phenotypes')
pdf('./heatmap_of_cor.pdf',width = 9, height = 8)
heatmap_object = draw(heatmap_object)
dev.off()

row_order(heatmap_object)
col_order(heatmap_object)

# 获取聚类结果
row_order <- heatmap_object@row_order
col_order <- heatmap_object@column_order

# 根据聚类结果提取关系较强的区块
strong_relation_rows <- cov_yy[row_order, ]
strong_relation_cols <- strong_relation_rows[, col_order]

# 显示提取的区块
print(strong_relation_cols)


##### clustering
# 对协方差矩阵进行层次聚类
row_dist <- dist(cov_yy)
col_dist <- dist(t(cov_yy))
row_hclust <- hclust(row_dist)
col_hclust <- hclust(col_dist)


# 选择一个高度值
height_threshold <- 2

# 在树状图中，选择具有较强关系的行和列
strong_relation_row_clusters <- cutree(row_hclust, h = height_threshold)

# 绘制行和列的树状图
par(mfrow=c(1,2))
plot(row_hclust, main = "Row Dendrogram", sub = "", xlab = "", labels = F)
summary(factor(strong_relation_row_clusters))
hist(strong_relation_row_clusters,breaks = 200)

#  related
phe = 'X1349.0.0'
id = strong_relation_row_clusters[phe]
id = 5
summary(factor(strong_relation_row_clusters))[id]

# 找到具有较强关系的行和列的索引
strong_relation_rows <- which(strong_relation_row_clusters == id)
strong_relation_rows
# 提取具有较强关系的行和列
cov_yy[strong_relation_rows, phe]
strong_relation_block <- cov_yy[strong_relation_rows, strong_relation_rows]
# 显示提取的区块
print(strong_relation_block)

# Allel freq plot
freq = fread('./UKB_14M_MAF.afreq')
MAF = sapply(freq$ALT_FREQS,function(x) min(x,1-x))
pdf('./allele_freq.pdf')
hist(MAF,breaks=100)
dev.off()


## Analysis results

setwd('/public3/project_users/chengb/hjc/projects/UKBioCoin/analysis/')
color_plates = c("#0071C2","#D75615","#EDB11A")
library(glue)
library(ggplot2)
library(cowplot)
library(data.table)
phe = 'X4079.0.0'
for (PC_num in c(0,5,10)) {
  eval(parse(text=glue('plink_res_{PC_num}pc=fread(\'./plink_results/plink-{phe}-{PC_num}PC.{phe}.glm.linear\')')))
  eval(parse(text=glue('plink_res_{PC_num}pc=plink_res_{PC_num}pc[plink_res_{PC_num}pc$TEST==\'ADD\',]')))
  eval(parse(text=glue('oath_res_{PC_num}pc=fread(\'./oath_results/c_oath-{phe}-{PC_num}PC_results.table\')')))
  print(PC_num)
}


setwd(glue('/public3/project_users/chengb/hjc/projects/UKBioCoin/analysis/plots{phe}_c'))
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
    scale_x_continuous(name = 'PLINK estimates',limits = c(-1,1)) + 
    scale_y_continuous(name = 'UKC estimates',limits = c(-1,1)) + 
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
    scale_x_continuous(name = 'PLINK estimates',limits = c(0,50)) + 
    scale_y_continuous(name = 'UKC estimates',limits = c(0,50)) + 
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

# overlap analysis
cut = -log(0.5/dim(oath_res_0pc)[1],10)
for (i in c(1:3)) {
  PC_num = c(0,5,10)[i]
  eval(parse(text=glue('plink_res_{PC_num}pc_cut=plink_res_{PC_num}pc[-log(plink_res_{PC_num}pc$P,10)>cut,]')))
  eval(parse(text=glue('oath_res_{PC_num}pc_cut=oath_res_{PC_num}pc[oath_res_{PC_num}pc$`-log10_P`>cut,]')))
}
for (i in c(1:3)) {
  PC_num = c(0,5,10)[i]
  eval(parse(text=glue('print(paste0(\'plink-oath: \', length(setdiff(plink_res_{PC_num}pc_cut$ID,oath_res_{PC_num}pc_cut$ID))))')))
  eval(parse(text=glue('print(paste0(\'oath-plink: \', length(setdiff(oath_res_{PC_num}pc_cut$ID,plink_res_{PC_num}pc_cut$ID))))')))
  eval(parse(text=glue('print(length(intersect(plink_res_{PC_num}pc_cut$ID,oath_res_{PC_num}pc_cut$ID)))')))
}

# Mahatten plot
library(qqman)
oath_res_0pc$`-log10_P` = oath_res_0pc$`-log10_P`

# ids = 1:dim(oath_res_0pc)[1]
jpeg('./manhattan_oath.jpg', height = 1500, width = 1000, quality = 100, res = 100)
par(mfrow=c(3,1))
manhattan(data.frame(CHR = oath_res_0pc[ids]$`#CHROM`,
                     BP = oath_res_0pc[ids]$`POS`,
                     SNP = oath_res_0pc[ids]$`ID`,
                     P = 10^-oath_res_0pc[ids]$`-log10_P`),
          col = c('#30A9DE','#EFDC05','#E53A40','#090707'),#交替使用颜色展示
          suggestiveline = -log10(1e-05),#－log10(1e－5)处添加"suggestive"横线
          genomewideline = -log10(5e-08),#－log10(5e－10)处添加"genome-wide sigificant"横线
          # annotatePval = 0.05,#标记p值小于0.05的点
          # annotateTop = T,#如果为T，则仅批注低于注解阈值的每个染色体上的顶部点，为F则标记所有小于注解阈值的点。
          main = "oath 0pc",
          ylim=c(0,40))
manhattan(data.frame(CHR = oath_res_5pc[ids]$`#CHROM`,
                         BP = oath_res_5pc[ids]$`POS`,
                         SNP = oath_res_5pc[ids]$`ID`,
                         P = 10^-oath_res_5pc[ids]$`-log10_P`),
              col = c('#30A9DE','#EFDC05','#E53A40','#090707'),#交替使用颜色展示
              suggestiveline = -log10(1e-05),#－log10(1e－5)处添加"suggestive"横线
              genomewideline = -log10(5e-08),#－log10(5e－10)处添加"genome-wide sigificant"横线
              # annotatePval = 0.05,#标记p值小于0.05的点
              # annotateTop = T,#如果为T，则仅批注低于注解阈值的每个染色体上的顶部点，为F则标记所有小于注解阈值的点。
              main = "oath 5pc",
              ylim=c(0,40))
manhattan(data.frame(CHR = oath_res_10pc[ids]$`#CHROM`,
                     BP = oath_res_10pc[ids]$`POS`,
                     SNP = oath_res_10pc[ids]$`ID`,
                     P = 10^-oath_res_10pc[ids]$`-log10_P`),
          col = c('#30A9DE','#EFDC05','#E53A40','#090707'),#交替使用颜色展示
          suggestiveline = -log10(1e-05),#－log10(1e－5)处添加"suggestive"横线
          genomewideline = -log10(5e-08),#－log10(5e－10)处添加"genome-wide sigificant"横线
          # annotatePval = 0.05,#标记p值小于0.05的点
          # annotateTop = T,#如果为T，则仅批注低于注解阈值的每个染色体上的顶部点，为F则标记所有小于注解阈值的点。
          main = "oath 10pc",
          ylim=c(0,40))
dev.off()

jpeg('./manhattan_plink.jpg', height = 1500, width = 1000, quality = 100, res = 100)
par(mfrow=c(3,1))
manhattan(data.frame(CHR = plink_res_0pc[ids]$`#CHROM`,
                     BP = plink_res_0pc[ids]$`POS`,
                     SNP = plink_res_0pc[ids]$`ID`,
                     P = plink_res_0pc[ids]$`P`),
          col = c('#30A9DE','#EFDC05','#E53A40','#090707'),#交替使用颜色展示
          suggestiveline = -log10(1e-05),#－log10(1e－5)处添加"suggestive"横线
          genomewideline = -log10(5e-08),#－log10(5e－10)处添加"genome-wide sigificant"横线
          # annotatePval = 0.05,#标记p值小于0.05的点
          # annotateTop = T,#如果为T，则仅批注低于注解阈值的每个染色体上的顶部点，为F则标记所有小于注解阈值的点。
          main = "plink 0pc",
          ylim=c(0,40))
manhattan(data.frame(CHR = plink_res_5pc[ids]$`#CHROM`,
                     BP = plink_res_5pc[ids]$`POS`,
                     SNP = plink_res_5pc[ids]$`ID`,
                     P = plink_res_5pc[ids]$`P`),
          col = c('#30A9DE','#EFDC05','#E53A40','#090707'),#交替使用颜色展示
          suggestiveline = -log10(1e-05),#－log10(1e－5)处添加"suggestive"横线
          genomewideline = -log10(5e-08),#－log10(5e－10)处添加"genome-wide sigificant"横线
          # annotatePval = 0.05,#标记p值小于0.05的点
          # annotateTop = T,#如果为T，则仅批注低于注解阈值的每个染色体上的顶部点，为F则标记所有小于注解阈值的点。
          main = "plink 5pc",
          ylim=c(0,40))
manhattan(data.frame(CHR = plink_res_10pc[ids]$`#CHROM`,
                     BP = plink_res_10pc[ids]$`POS`,
                     SNP = plink_res_10pc[ids]$`ID`,
                     P = plink_res_10pc[ids]$`P`),
          col = c('#30A9DE','#EFDC05','#E53A40','#090707'),#交替使用颜色展示
          suggestiveline = -log10(1e-05),#－log10(1e－5)处添加"suggestive"横线
          genomewideline = -log10(5e-08),#－log10(5e－10)处添加"genome-wide sigificant"横线
          # annotatePval = 0.05,#标记p值小于0.05的点
          # annotateTop = T,#如果为T，则仅批注低于注解阈值的每个染色体上的顶部点，为F则标记所有小于注解阈值的点。
          main = "plink 10pc",
          ylim=c(0,40))
dev.off()

task_name=paste(task_prefix,'-',phe_name,'-',PC_num,'PC',sep='')

## analysis outlier SNPs

fwrite(phes[!all_ids$V1%in%setdiff(all_ids$V1,sub_ids$IID)],"/public3/project_users/chengb/hjc/projects/MR/pheno/valid_584pheno_20pcs_withcolnames_scaled.table",sep=' ',na='NA', col.names = T, quote=F)

out = (oa$SE>0.8)
out_snps = oa[out]$ID[1:100]

# 
plink   --bfile /public3/project_users/chengb/hjc/data/UKB_white_Impute_10M_BF \
--extract ./outlier_ids.txt \
--out outlier \
--recode 
out_ped = fread("/public3/project_users/chengb/hjc/projects/UKBioCoin/analysis/outlier.ped")
out_map = fread("/public3/project_users/chengb/hjc/projects/UKBioCoin/analysis/outlier.map")
meta = fread()
out_meta = meta[out][1:100]

phes = fread("/public3/project_users/chengb/hjc/projects/MR/pheno/584pheno_20pcs_withcolnames_scaled.table")[1:292216]
phe = phes[,c('X4079.0.0','X31.0.0')]

out_ped = as.matrix(out_ped)
data = matrix(NA,dim(out_ped)[1],100)
for (i in c(1:100)) {
  data[,i] = as.numeric(out_ped[,(6+2*i-1)]==out_meta[i]$ALT) + as.numeric(out_ped[,(6+2*i)]==out_meta[i]$ALT)
  print(i)
}
data[out_ped[,seq(7,206,2)]==0] = NA

X = cbind(data[,3],phes$X31.0.0)
Y = phes$X4079.0.0
cases = complete.cases(cbind(X,Y))
X = X[cases,]
Y = Y[cases]
solve(t(X)%*%X)%*%(t(X)%*%Y)

oa[out][1:5]

lm(formula = y ~ x+z, data = data.frame(y=phes$X4079.0.0,z=phes$X31.0.0,x=data[,1]))
