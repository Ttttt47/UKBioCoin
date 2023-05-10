library(data.table)
library(ggplot2)
library(glue)
library(cowplot)
library(data.table)
color_plates = c("#0071C2","#D75615","#EDB11A")

oath_res = fread("/public3/project_users/chengb/hjc/projects/UKBioCoin/0422/analysis/binary_out/oath-smoking_results.table")
plink_res = fread("/public3/project_users/chengb/hjc/projects/UKBioCoin/0422/analysis/binary_out/plink-X20160.0.0-5PC_nonscaled.X20160.0.0.glm.logistic.hybrid")
cov_yy = fread('/public3/project_users/chengb/hjc/projects/UKBioCoin/0422/2.matrix/10M_cov_yy.table')
phe = fread("/public3/project_users/chengb/hjc/projects/MR/pheno/129_filtered_pheno_20pcs_non_scaled_0421.table")
std_err = sd(phe$X20160.0.0,na.rm=T)
setwd('/public3/project_users/chengb/hjc/projects/UKBioCoin/0422/analysis/binary_out/')



jpeg('./OR.vs.oath.jpg', height = 900, width = 900, quality = 100, res = 300)
ggplot(data.frame(x=plink_res$OR, y=oath_res$BETA*std_err), aes(x = x, y = y)) +
  scale_x_continuous(name = 'PLINK OR',limits = c(-1,1)) + 
  scale_y_continuous(name = 'UKC BETA estimates',limits = c(-1,1)) + 
  geom_point(color = color_plates[i], size=0.3, alpha = 0.2)
# annotate("text", x = -1.5, y = 1.5, label = glue("{cor}"), size = 4) +
geom_abline(intercept=-1, slope=1, linewidth= 0.4, color='grey', alpha=0.7)
dev.off()
p1 = ggplot(data.frame(x=plink_points,y=oath_points), aes(x = x, y = y)) +
  scale_x_continuous(name = 'PLINK estimates',limits = c(-1,1)) + 
  scale_y_continuous(name = 'UKC estimates',limits = c(-1,1)) + 
  geom_point(color = color_plates[i], size=0.3, alpha = 0.2) +
  ggtitle(glue('Beta, {PC_num} PC', )) +
  # annotate("text", x = -1.5, y = 1.5, label = glue("{cor}"), size = 4) +
  geom_abline(intercept=0, slope=1, linewidth= 0.4, color='grey', alpha=0.7)

p2 = ggplot(data.frame(x=plink_res$OR, y=oath_res$BETA*std_err), aes(x = x, y = y)) +
  scale_x_continuous(name = 'PLINK OR',limits = c(-1,1)) + 
  scale_y_continuous(name = 'UKC BETA estimates',limits = c(-1,1)) + 
  geom_point(color = color_plates[i], size=0.3, alpha = 0.2) +
  ggtitle(glue('Beta, {PC_num} PC', )) +
  # annotate("text", x = -1.5, y = 1.5, label = glue("{cor}"), size = 4) +
  geom_abline(intercept=0, slope=1, linewidth= 0.4, color='grey', alpha=0.7)
# theme_minimal()
plot_list[[length(plot_list)+1]] = p

idx = sample(c(1:dim(plink_res)[1]),10^5)
jpeg('./OR.vs.oath.jpg', height = 900, width = 900, quality = 100, res = 300)
ggplot(data.frame(x=plink_res$OR[idx], y=(oath_res$BETA/std_err)[idx]), aes(x = x, y = y)) +
  scale_x_continuous(name = 'PLINK OR',limits = c(-1,1)) + 
  scale_y_continuous(name = 'UKC BETA estimates',limits = c(-1,1)) + 
  geom_point(color = color_plates[1], size=0.3, alpha = 0.2)+
  # annotate("text", x = -1.5, y = 1.5, label = glue("{cor}"), size = 4) +
  geom_abline(intercept=-1, slope=1, linewidth= 0.4, color='grey', alpha=0.7)
dev.off()

jpeg('./pvalue.jpg', height = 900, width = 900, quality = 100, res = 300)
ggplot(data.frame(x=-log(plink_res$P,10)[idx], y=(oath_res$`-log10_P`)[idx]), aes(x = x, y = y)) +
  scale_x_continuous(name = 'PLINK -log10 p-value',limits = c(0,50)) + 
  scale_y_continuous(name = 'UKC -log10 p-value',limits = c(0,50)) + 
  geom_point(color = color_plates[1], size=0.3, alpha = 0.2)+
# annotate("text", x = -1.5, y = 1.5, label = glue("{cor}"), size = 4) +
geom_abline(intercept=0, slope=1, linewidth= 0.4, color='grey', alpha=0.7)
dev.off()