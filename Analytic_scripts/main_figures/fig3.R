# Figure 3 (Pertubation under different covariates).

library(ggplot2)
library(reshape2)
library(circlize)
library(data.table)
library(glue)
library(cowplot)


phe = 'X21002.0.0'  # weight
threshold = 0.05/10531641

# pertubation under sex ascovariate------------------------------------------
ukc_5pc = fread(glue('../single_reg_with_PC_vif/ukc-{phe}-5ukbPC_results.table'), header=T)
ukc_5pc_sex = fread(glue('../adjusting_sex_vif/ukc-{phe}-basic_adj_results.table'), header=T)
ukc_5pc_sexage = fread(glue('../adjusting_sexage_vif/ukc-{phe}-basic_adj_results.table'), header=T)
ukc_5pc_BMI = fread(glue('../adjusting_BMI_vif/ukc-{phe}-basic_adj_results.table'), header=T)
ukc_raw = fread(glue('../single_reg_with_PC_vif/ukc-{phe}-0PC_results.table'), header=T)
ukc_sex = fread(glue('../adjusting_only_sex_or_BMI/ukc-{phe}-sex_results.table'), header=T)
ukc_BMI = fread(glue('../adjusting_only_sex_or_BMI/ukc-{phe}-BMI_results.table'), header=T)

# plot a scatter plot of the significant SNPs' -log10(p) in the two models
vif_threshold = 50
sig_snp1 = which(ukc_5pc$`-log10_P` > -log10(threshold) & ukc_5pc_sex$VIF < vif_threshold)
sig_snp2 = which(ukc_5pc_sex$`-log10_P` > -log10(threshold) & ukc_5pc_sex$VIF < vif_threshold)
sig_snp3 = which(ukc_5pc_sexage$`-log10_P` > -log10(threshold) & ukc_5pc_sexage$VIF < vif_threshold)
sig_snp_5pc_BMI = which(ukc_5pc_BMI$`-log10_P` > -log10(threshold) & ukc_5pc_BMI$VIF < vif_threshold)
sig_snp_raw = which(ukc_raw$`-log10_P` > -log10(threshold) & ukc_raw$VIF < vif_threshold)
sig_snp_sex = which(ukc_sex$`-log10_P` > -log10(threshold) & ukc_sex$VIF < vif_threshold)
sig_snp_BMI = which(ukc_BMI$`-log10_P` > -log10(threshold) & ukc_BMI$VIF < vif_threshold)

# save QTLs
system("mkdir -p QTLs")
write.table(ukc_5pc[sig_snp1], paste0(glue("./QTLs/1.QTL_NoClumping.ukc.{phe}_"), '5pc', ".table"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
write.table(ukc_5pc_sex[sig_snp2], paste0(glue("./QTLs/1.QTL_NoClumping.ukc.{phe}_"), '5pc_sex', ".table"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
write.table(ukc_5pc_sexage[sig_snp3], paste0(glue("./QTLs/1.QTL_NoClumping.ukc.{phe}_"), '5pc_sexage', ".table"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
write.table(ukc_raw[sig_snp_raw], paste0(glue("./QTLs/1.QTL_NoClumping.ukc.{phe}_"), 'raw', ".table"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
write.table(ukc_sex[sig_snp_sex], paste0(glue("./QTLs/1.QTL_NoClumping.ukc.{phe}_"), 'sex', ".table"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
write.table(ukc_BMI[sig_snp_BMI], paste0(glue("./QTLs/1.QTL_NoClumping.ukc.{phe}_"), 'BMI', ".table"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")


## raw v.s. sex
# p-values
sig_snp = union(sig_snp_raw, sig_snp_sex)
highlight = which(ukc_raw$`-log10_P`[sig_snp] < ukc_sex$`-log10_P`[sig_snp])
lr = lm(ukc_sex$`-log10_P` ~ ukc_raw$`-log10_P`)
cor = cor(ukc_sex$`-log10_P`,ukc_raw$`-log10_P`, use='complete.obs')
cor = round(cor,digits = 3)
p11 = ggplot() +
    geom_point(data = data.frame(x=ukc_raw$`-log10_P`[sig_snp][highlight], y=ukc_sex$`-log10_P`[sig_snp][highlight]), 
                aes(x=x, y=y), size = 0.5, alpha = 0.5, color = 'red') +
    labs(x = expression(-log[10]*'(p) without adjustment'), y = expression(-log[10]*'(p) adjust for sex')) +
    geom_point(data = data.frame(x=ukc_raw$`-log10_P`[sig_snp][-highlight], y=ukc_sex$`-log10_P`[sig_snp][-highlight]), 
                aes(x=x, y=y), size = 0.5, alpha = 0.5, color = 'black') +
    # reference line
    geom_abline(intercept=0, slope=1, linewidth= 0.8, color='grey', alpha=0.7) +
    geom_abline(intercept=lr$coeff[[1]], slope=lr$coeff[[2]], linewidth= 0.8, color='black', alpha=0.5) +
    theme(panel.background = element_rect(fill='transparent',colour = NA),
                        panel.grid = element_blank(),
                        panel.border = element_blank(),
                        axis.line = element_line(colour = "black",linewidth = 0.5),
                        axis.ticks = element_line(colour = "black",linewidth = 0.5))
ggsave(filename = glue("./{phe}_adjust_raw_vs_sex_pvalues.png"), plot = p11, width = 1000, height = 1000, units = 'px')
print(glue("y={lr$coeff[[2]]}x+{lr$coeff[[1]]},cor={cor}"))
# betas
lr = lm(ukc_sex$BETA ~ ukc_raw$BETA)
cor = cor(ukc_sex$BETA,ukc_raw$BETA, use='complete.obs')
cor = round(cor,digits = 3)
p12 = ggplot() +
    labs(x = expression(beta*' without adjustment'), y = expression(beta*' adjust for sex')) +
    geom_point(data = data.frame(x=ukc_raw$`BETA`[sig_snp], y=ukc_sex$`BETA`[sig_snp]), 
                aes(x=x, y=y), size = 0.5, alpha = 0.5, color = 'black') +
    # reference line
    geom_abline(intercept=0, slope=1, linewidth= 0.8, color='grey', alpha=0.7) +
    geom_abline(intercept=lr$coeff[[1]], slope=lr$coeff[[2]], linewidth= 0.8, color='black', alpha=0.5) +
    theme(panel.background = element_rect(fill='transparent',colour = NA),
                        panel.grid = element_blank(),
                        panel.border = element_blank(),
                        axis.line = element_line(colour = "black",linewidth = 0.5),
                        axis.ticks = element_line(colour = "black",linewidth = 0.5))
ggsave(filename = glue("./{phe}_adjust_raw_vs_sex_beta.png"), plot = p12, width = 1000, height = 1000, units = 'px')
print(glue("y={lr$coeff[[2]]}x+{lr$coeff[[1]]},cor={cor}"))

## raw v.s. BMI
sig_snp = union(sig_snp_raw, sig_snp_BMI)
highlight = which(ukc_raw$`-log10_P`[sig_snp] < ukc_BMI$`-log10_P`[sig_snp])
lr = lm(ukc_BMI$`-log10_P` ~ ukc_raw$`-log10_P`)
cor = cor(ukc_BMI$`-log10_P`,ukc_raw$`-log10_P`, use='complete.obs')
cor = round(cor,digits = 3)
p21 = ggplot() +
    geom_point(data = data.frame(x=ukc_raw$`-log10_P`[sig_snp][highlight], y=ukc_BMI$`-log10_P`[sig_snp][highlight]), 
                aes(x=x, y=y), size = 0.5, alpha = 0.5, color = 'red') +
    labs(x = expression(-log[10]*'(p) without adjustment'), y = expression(-log[10]*'(p) adjust for BMI')) +
    geom_point(data = data.frame(x=ukc_raw$`-log10_P`[sig_snp][-highlight], y=ukc_BMI$`-log10_P`[sig_snp][-highlight]), 
                aes(x=x, y=y), size = 0.5, alpha = 0.5, color = 'black') +
    # reference line
    geom_abline(intercept=0, slope=1, linewidth= 0.8, color='grey', alpha=0.7) +
    geom_abline(intercept=lr$coeff[[1]], slope=lr$coeff[[2]], linewidth= 0.8, color='black', alpha=0.5) +
    theme(panel.background = element_rect(fill='transparent',colour = NA),
                        panel.grid = element_blank(),
                        panel.border = element_blank(),
                        axis.line = element_line(colour = "black",linewidth = 0.5),
                        axis.ticks = element_line(colour = "black",linewidth = 0.5))
ggsave(filename = glue("./{phe}_adjust_raw_vs_BMI_pvalues.png"), plot = p21, width = 1000, height = 1000, units = 'px')
print(glue("y={lr$coeff[[2]]}x+{lr$coeff[[1]]},cor={cor}"))
# betas
lr = lm(ukc_BMI$BETA ~ ukc_raw$BETA)
cor = cor(ukc_BMI$BETA,ukc_raw$BETA, use='complete.obs')
cor = round(cor,digits = 3)
p22 = ggplot() +
    labs(x = expression(beta*' without adjustment'), y = expression(beta*' adjust for BMI')) +
    geom_point(data = data.frame(x=ukc_raw$`BETA`[sig_snp], y=ukc_BMI$`BETA`[sig_snp]), 
                aes(x=x, y=y), size = 0.5, alpha = 0.5, color = 'black') +
    # reference line
    geom_abline(intercept=0, slope=1, linewidth= 0.8, color='grey', alpha=0.7) +
    geom_abline(intercept=lr$coeff[[1]], slope=lr$coeff[[2]], linewidth= 0.8, color='black', alpha=0.5) +
    theme(panel.background = element_rect(fill='transparent',colour = NA),
                        panel.grid = element_blank(),
                        panel.border = element_blank(),
                        axis.line = element_line(colour = "black",linewidth = 0.5),
                        axis.ticks = element_line(colour = "black",linewidth = 0.5))
ggsave(filename = glue("./{phe}_adjust_raw_vs_BMI_beta.png"), plot = p22, width = 1000, height = 1000, units = 'px')
print(glue("y={lr$coeff[[2]]}x+{lr$coeff[[1]]},cor={cor}"))

## raw v.s. 5pc
sig_snp = union(sig_snp_raw, sig_snp1)
highlight = which(ukc_raw$`-log10_P`[sig_snp] < ukc_5pc$`-log10_P`[sig_snp])
lr = lm(ukc_5pc$`-log10_P` ~ ukc_raw$`-log10_P`)
cor = cor(ukc_5pc$`-log10_P`,ukc_raw$`-log10_P`, use='complete.obs')
cor = round(cor,digits = 3)
p31 = ggplot() +
    geom_point(data = data.frame(x=ukc_raw$`-log10_P`[sig_snp][highlight], y=ukc_5pc$`-log10_P`[sig_snp][highlight]), 
                aes(x=x, y=y), size = 0.5, alpha = 0.5, color = 'red') +
    labs(x = expression(-log[10]*'(p) without adjustment'), y = expression(-log[10]*'(p) adjust for 5 PCs')) +
    geom_point(data = data.frame(x=ukc_raw$`-log10_P`[sig_snp][-highlight], y=ukc_5pc$`-log10_P`[sig_snp][-highlight]), 
                aes(x=x, y=y), size = 0.5, alpha = 0.5, color = 'black') +
    # reference line
    geom_abline(intercept=0, slope=1, linewidth= 0.8, color='grey', alpha=0.7) +
    geom_abline(intercept=lr$coeff[[1]], slope=lr$coeff[[2]], linewidth= 0.8, color='black', alpha=0.5) +
    theme(panel.background = element_rect(fill='transparent',colour = NA),
                        panel.grid = element_blank(),
                        panel.border = element_blank(),
                        axis.line = element_line(colour = "black",linewidth = 0.5),
                        axis.ticks = element_line(colour = "black",linewidth = 0.5))
ggsave(filename = glue("./{phe}_adjust_raw_vs_5pc_pvalues.png"), plot = p31, width = 1000, height = 1000, units = 'px')
print(glue("y={lr$coeff[[2]]}x+{lr$coeff[[1]]},cor={cor}"))
# betas
lr = lm(ukc_5pc$BETA ~ ukc_raw$BETA)
cor = cor(ukc_5pc$BETA,ukc_raw$BETA, use='complete.obs')
cor = round(cor,digits = 3)
p32 = ggplot() +
    labs(x = expression(beta*' without adjustment'), y = expression(beta*' adjust for 5 PCs')) +
    geom_point(data = data.frame(x=ukc_raw$`BETA`[sig_snp], y=ukc_5pc$`BETA`[sig_snp]), 
                aes(x=x, y=y), size = 0.5, alpha = 0.5, color = 'black') +
    # reference line
    geom_abline(intercept=0, slope=1, linewidth= 0.8, color='grey', alpha=0.7) +
    geom_abline(intercept=lr$coeff[[1]], slope=lr$coeff[[2]], linewidth= 0.8, color='black', alpha=0.5) +
    theme(panel.background = element_rect(fill='transparent',colour = NA),
                        panel.grid = element_blank(),
                        panel.border = element_blank(),
                        axis.line = element_line(colour = "black",linewidth = 0.5),
                        axis.ticks = element_line(colour = "black",linewidth = 0.5))
ggsave(filename = glue("./{phe}_adjust_raw_vs_5pc_beta.png"), plot = p32, width = 1000, height = 1000, units = 'px')
print(glue("y={lr$coeff[[2]]}x+{lr$coeff[[1]]},cor={cor}"))

p = plot_grid(plotlist = list(p11, p12, p21, p22, p31, p32), nrow = 2, ncol = 3, align = 'h', axis = 'bt', byrow = F)
ggsave(filename = glue("./{phe}_adjust_all_5ukbpc_sex_BMI.png"), plot = p, width = 3000, height = 2000, units = 'px')



## 5pc v.s. 5pc+sex
sig_snp = union(sig_snp1, sig_snp2)
# highlight those points with smaller p-value in model adjust for other covariates
highlight = which(ukc_5pc$`-log10_P`[sig_snp] < ukc_5pc_sex$`-log10_P`[sig_snp])
lr = lm(ukc_5pc_sex$`-log10_P` ~ ukc_5pc$`-log10_P`)
cor = cor(ukc_5pc_sex$`-log10_P`,ukc_5pc$`-log10_P`, use='complete.obs')
cor = round(cor,digits = 3)
p = ggplot() +
    geom_point(data = data.frame(x=ukc_5pc$`-log10_P`[sig_snp][highlight], y=ukc_5pc_sex$`-log10_P`[sig_snp][highlight]), 
                aes(x=x, y=y), size = 0.5, alpha = 0.5, color = 'red') +
    labs(x = expression(-log[10]*'(p) adjust for 5 PCs'), y = expression(-log[10]*'(p) adjust for 5 PCs and sex')) +
    geom_point(data = data.frame(x=ukc_5pc$`-log10_P`[sig_snp][-highlight], y=ukc_5pc_sex$`-log10_P`[sig_snp][-highlight]), 
                aes(x=x, y=y), size = 0.5, alpha = 0.5, color = 'black') +
    # reference line
    geom_abline(intercept=0, slope=1, linewidth= 0.8, color='grey', alpha=0.7) +
    geom_abline(intercept=lr$coeff[[1]], slope=lr$coeff[[2]], linewidth= 0.8, color='black', alpha=0.5) +
    theme(panel.background = element_rect(fill='transparent',colour = NA),
                        panel.grid = element_blank(),
                        panel.border = element_blank(),
                        axis.line = element_line(colour = "black",linewidth = 0.5),
                        axis.ticks = element_line(colour = "black",linewidth = 0.5))
ggsave(filename = glue("./{phe}_adjust_sex_pvalues.png"), plot = p, width = 1000, height = 1000, units = 'px')
print(glue("y={lr$coeff[[2]]}x+{lr$coeff[[1]]},cor={cor}"))

## 5pc v.s. 5pc+BMI
sig_snp = union(sig_snp1, sig_snp_5pc_BMI)
# highlight those points with smaller p-value in model adjust for other covariates
highlight = which(ukc_5pc$`-log10_P`[sig_snp] < ukc_5pc_BMI$`-log10_P`[sig_snp])
lr = lm(ukc_5pc_BMI$`-log10_P` ~ ukc_5pc$`-log10_P`)
cor = cor(ukc_5pc_BMI$`-log10_P`,ukc_5pc$`-log10_P`, use='complete.obs')
cor = round(cor,digits = 3)
print(glue("y={lr$coeff[[2]]}x+{lr$coeff[[1]]},cor={cor}"))
p = ggplot() +
    geom_point(data = data.frame(x=ukc_5pc$`-log10_P`[sig_snp][highlight], y=ukc_5pc_BMI$`-log10_P`[sig_snp][highlight]), 
                aes(x=x, y=y), size = 0.5, alpha = 0.5, color = 'red') +
    labs(x = expression(-log[10]*'(p) adjust for 5 PCs'), y = expression(-log[10]*'(p) adjust for 5 PCs and BMI')) +
    geom_point(data = data.frame(x=ukc_5pc$`-log10_P`[sig_snp][-highlight], y=ukc_5pc_BMI$`-log10_P`[sig_snp][-highlight]), 
                aes(x=x, y=y), size = 0.5, alpha = 0.5, color = 'black') +
    # reference line
    geom_abline(intercept=0, slope=1, linewidth= 0.8, color='grey', alpha=0.7) +
    geom_abline(intercept=lr$coeff[[1]], slope=lr$coeff[[2]], linewidth= 0.8, color='black', alpha=0.5) +
    theme(panel.background = element_rect(fill='transparent',colour = NA),
                        panel.grid = element_blank(),
                        panel.border = element_blank(),
                        axis.line = element_line(colour = "black",linewidth = 0.5),
                        axis.ticks = element_line(colour = "black",linewidth = 0.5))
ggsave(filename = glue("./{phe}_adjust_BMI_pvalues.png"), plot = p, width = 1000, height = 1000, units = 'px')

# plot a histogram of the changes of -log10(p) in the two models
p = ggplot() +
    geom_histogram(data = data.frame(x=ukc_5pc_sex$`-log10_P`[sig_snp] - ukc_5pc$`-log10_P`[sig_snp]), 
                aes(x=x), binwidth = 0.1, fill = 'black', alpha = 0.8) +
    labs(x =  expression("Changes of "*-log[10]*'(p)'), y = 'Count') +
    theme(panel.background = element_rect(fill='transparent',colour = NA),
                        panel.grid = element_blank(),
                        panel.border = element_blank(),
                        axis.line = element_line(colour = "black",linewidth = 0.5),
                        axis.ticks = element_line(colour = "black",linewidth = 0.5))
ggsave(filename = glue("./{phe}_adjust_sex_pvalues_hist.png"), plot = p, width = 2000, height = 1000, units = 'px')
print(glue("y={lr$coeff[[2]]}x+{lr$coeff[[1]]},cor={cor}"))


library(ggVennDiagram)
library(ggplot2)
library(cowplot)
library(glue)
library(data.table)

# draw a histogram of the height(X50.0.0) and BMI(X21001.0.0) separately for female (X31.0.0==0) and male (X31.0.0==1)
# height separated into two groups
p1 = ggplot() + 
    geom_density(data=phes, aes(x=X50.0.0, group=sex, fill=sex),alpha=0.5, adjust=2) +
    labs(x =  expression("height"), y = 'density') +
    theme(panel.background = element_rect(fill='transparent',colour = NA),
                        panel.grid = element_blank(),
                        panel.border = element_blank(),
                        axis.line = element_line(colour = "black",linewidth = 0.5),
                        axis.ticks = element_line(colour = "black",linewidth = 0.5),
                        legend.position = c(0.2,0.7),legend.title=element_blank())
# BMI
p2 = ggplot() +
    geom_density(data=phes, aes(x=X21001.0.0, group=sex, fill=sex),alpha=0.5, adjust=2) +
    labs(x = expression('BMI'), y = 'density') +
    theme(panel.background = element_rect(fill='transparent',colour = NA),
                        panel.grid = element_blank(),
                        panel.border = element_blank(),
                        axis.line = element_line(colour = "black",linewidth = 0.5),
                        axis.ticks = element_line(colour = "black",linewidth = 0.5),
                        legend.position = c(0.9,0.7),legend.title=element_blank())
# scatter plot of height and BMI
p3 = ggplot(data = phes, aes(x=X21001.0.0, y=X50.0.0)) +
    geom_point(size = 0.5, alpha = 0.5, color = 'black') +
    labs(x = expression('BMI'), y = expression('height')) +
    geom_smooth(method = "lm", se = FALSE, color = "grey", alpha=0.5) +
    theme(panel.background = element_rect(fill='transparent',colour = NA),
                        panel.grid = element_blank(),
                        panel.border = element_blank(),
                        axis.line = element_line(colour = "black",linewidth = 0.5),
                        axis.ticks = element_line(colour = "black",linewidth = 0.5))

jpeg(glue('./Height_BMI_dist.jpg'), height = 1000, width = 1000*3, quality = 100, res = 300)
plot_grid(plotlist = list(p1, p2, p3), nrow = 1, ncol = 3, align = 'h', axis = 'bt', byrow = F)
dev.off()

p1 = ggplot() + 
    geom_density(data=phes, aes(x=X21002.0.0, group=sex, fill=sex),alpha=0.5, adjust=2) +
    labs(x =  expression("weight"), y = 'density') +
    theme(panel.background = element_rect(fill='transparent',colour = NA),
                        panel.grid = element_blank(),
                        panel.border = element_blank(),
                        axis.line = element_line(colour = "black",linewidth = 0.5),
                        axis.ticks = element_line(colour = "black",linewidth = 0.5),
                        legend.position = c(0.9,0.7),legend.title=element_blank())
# BMI
p2 = ggplot() +
    geom_density(data=phes, aes(x=X21001.0.0, group=sex, fill=sex),alpha=0.5, adjust=2) +
    labs(x = expression('BMI'), y = 'density') +
    theme(panel.background = element_rect(fill='transparent',colour = NA),
                        panel.grid = element_blank(),
                        panel.border = element_blank(),
                        axis.line = element_line(colour = "black",linewidth = 0.5),
                        axis.ticks = element_line(colour = "black",linewidth = 0.5),
                        legend.position = c(0.9,0.7),legend.title=element_blank())
# scatter plot of weight and BMI
p3 = ggplot(data = phes, aes(x=X21001.0.0, y=X21002.0.0)) +
    geom_point(size = 0.5, alpha = 0.5, color = 'black') +
    labs(x = expression('BMI'), y = expression('weight')) +
    geom_smooth(method = "lm", se = FALSE, color = "grey", alpha=0.5) +
    theme(panel.background = element_rect(fill='transparent',colour = NA),
                        panel.grid = element_blank(),
                        panel.border = element_blank(),
                        axis.line = element_line(colour = "black",linewidth = 0.5),
                        axis.ticks = element_line(colour = "black",linewidth = 0.5))

jpeg(glue('./Weight_BMI_dist.jpg'), height = 1000, width = 1000*3, quality = 100, res = 300)
plot_grid(plotlist = list(p1, p2, p3), nrow = 1, ncol = 3, align = 'h', axis = 'bt', byrow = F)
dev.off()