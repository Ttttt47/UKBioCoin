# Figure 4 (heritability via LDSC).

# Fig 4a and 4b
library(ggplot2)
library(ggrepel)
library(ggsci)
out1 = read.table("LDSC.summary.UKC2024.txt", sep = '\t', header=T)
out1$Correction = sub("raw","Raw",out1$Correction)
out1 = out1[order(out1$Heritability,decreasing = F),]
out1$Field.Code = factor(out1$Field.Code,levels = unique(out1$Field.Code))
out1$Correction = factor(out1$Correction,levels = c("Raw","PC5","PC10",paste0("G",1:5)))

out1_new = data.frame()
for(i in 1:length(unique(out1$Field.Code))){
  code = unique(out1$Field.Code)[i]
  
  temp = out1[which(out1$Field.Code==code),]
  rawH2 = temp$Heritability[which(temp$Correction=='Raw')]
  temp$DEV = temp$Heritability-temp$Heritability[which(temp$Correction=='Raw')]
  temp$Heritability[which(temp$Correction=='Raw')] = rawH2
  out1_new = rbind(out1_new,temp)
}
out1_new_target = out1_new[which(abs(out1_new$DEV)>0.1),]
targetCode = as.character(unique(out1_new_target$Field.Code))
textDT = data.frame()
for(i in 1:length(targetCode)){
  code = targetCode[i]
  temp = out1_new[which(out1_new$Field.Code==code),]
  
  temp = temp[which(temp$Heritability==max(temp$Heritability)),]
  textDT = rbind(textDT,temp)
}

unnorm = unique(out1$Field[which(out1$Heritability>=0.75)])
p = 
  ggplot(data = out1, mapping = aes(x = Field.Code, y = Heritability, fill = Correction, color = Correction, group = Correction))+
  geom_point()+
  geom_line()+
  geom_text_repel(data = textDT,
                  mapping = aes(label = Field),
                  nudge_y = textDT$Heritability*1.03,
                  size = 2,
                  box.padding = 0.1,
                  point.padding = 0.1,
                  force  = 100,
                  segment.size  = 0.1,
                  segment.color = "grey50",
                  direction = "x")+
  scale_color_npg(alpha = 0.6)+
  scale_fill_npg(alpha = 0.6)+
  xlab("505 Traits")+
  ylab("SNP-heritability")+
  guides(names="")+
  scale_y_continuous(expand = c(0,0),breaks = seq(0,0.7,0.1),limits = c(0,0.75),
                     labels = round(seq(0,0.7,0.1),digits = 2))+
  theme(plot.background = element_rect(fill = 'transparent', colour = 'transparent'),
        panel.background = element_rect(fill = 'transparent', colour = 'transparent'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(family = 'serif',size=16),
        axis.title = element_text(size = 16),
        #axis.text.x = element_text(size = 13,angle = 90,vjust = 0.5),
        legend.text = element_text(size = 12),
        legend.position = 'none',
        axis.line = element_line(colour = 'black'),
        legend.background = element_rect(fill = 'transparent', colour = 'transparent'),
        legend.key = element_blank(),
        legend.title = element_blank(),
        panel.grid = element_blank())

ggsave("Figure.4a.pdf",plot = p,width = 10,height = 3)


library(ggplot2)
library(ggsci)
out1 = read.table("LDSC.summary.UKC2024.txt", sep = '\t', header=T)
out1$Correction = sub("raw","Raw",out1$Correction)
out1$Correction = sub("PC5","5 PCs",out1$Correction)
out1$Correction = sub("PC10","10 PCs",out1$Correction)

out1$Correction = factor(out1$Correction,levels = c("Raw","5 PCs","10 PCs",paste0("G",1:5)))
p = 
  ggplot(data = out1, mapping = aes(x = Correction, y = Heritability, fill = Correction))+
  geom_boxplot(width = 0.6)+
  scale_fill_npg(alpha = 0.6)+
  ylab("SNP-heritability")+
  guides(names="")+
  scale_y_continuous(expand = c(0,0),breaks = seq(0,0.7,0.1),limits = c(0,0.75),
                     labels = round(seq(0,0.7,0.1),digits = 2))+
  # scale_y_continuous(expand = c(0,0),breaks = seq(0,0.25,0.05),limits = c(0,0.3),
  #                    labels = c("0.00",'0.05','0.10','0.15','0.20','0.25'))+
  theme(plot.background = element_rect(fill = 'transparent', colour = 'transparent'),
        panel.background = element_rect(fill = 'transparent', colour = 'transparent'),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        text = element_text(family = 'serif',size=18,colour = 'black'),
        axis.title = element_text(size = 18,colour = 'black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16,angle = 315,vjust = 0.5,colour = 'black'),
        axis.text.y = element_text(colour = 'black'),
        legend.text = element_text(size = 14),
        axis.line = element_line(colour = 'black'),
        axis.ticks.x=element_blank(),
        legend.position = 'none',
        legend.background = element_rect(fill = 'transparent', colour = 'transparent'),
        legend.key = element_blank(),
        legend.title = element_blank(),
        panel.grid = element_blank())
p
ggsave("Figure.4b.pdf",plot = p,width = 8,height = 4.8)

# Fig 4c compared with nealers

library(ggsci)
library(data.table)
library(ggplot2)
neales = fread("../Heritability/ukb31063_h2_topline.02Oct2019.tsv.gz",header=T,sep ="\t")
neales$phenotype = sub("_irnt","",neales$phenotype,fixed=T)

out1 = read.table("LDSC.summary.UKC2024.txt", sep = '\t', header=T)
ukc = out1[which(out1$Correction=='G5'),]

dt = merge(ukc,neales,by.x='Field.ID',by.y='phenotype',all.x = T)
dt = dt[-which(is.na(dt$h2_liability)),]
dt = dt[,c(4,11)]
dt$DEV = ifelse(dt$Heritability>=dt$h2_observed,'T1','T2')
p = 
  ggplot(data = dt, mapping = aes(x = h2_observed, y = Heritability, color = DEV))+
  geom_point(size=2)+
  scale_color_npg()+
  xlab("SNP-heritability (Neale's Lab)")+
  ylab("SNP-heritability (UKC)")+
  scale_y_continuous(expand = c(0.02,0.02))+
  geom_abline(slope = 1, intercept = 0)+
  annotate('text', x = 0.2, y = 0.44, label = "Pearson's correlation = 0.97",family='serif',size = 5)+
  theme(plot.background = element_rect(fill = 'transparent', colour = 'transparent'),
        panel.background = element_rect(fill = 'transparent', colour = 'transparent'),
        text = element_text(family = 'serif',size=16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.line = element_line(colour = 'black'),
        legend.position = 'none',
        legend.background = element_rect(fill = 'transparent', colour = 'transparent'),
        legend.key = element_blank(),
        legend.title = element_blank(),
        panel.grid = element_blank())
ggsave("Figure.4c.pdf",plot = p,width = 4,height = 4)
