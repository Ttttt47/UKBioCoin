# Figure 5 (PGS).

# Fig 5a and 5b

library(xlsx)
DT1 = read.xlsx("Summary.PGS.202412.xlsx",sheetIndex = 1)
DT1 = DT1[-which(DT1$FieldID==129 | DT1$FieldID==130 | DT1$FieldID==20074 | DT1$FieldID==20075),]
prsD = DT1

UKC_D = prsD[which(prsD$Correction.x=='UKC_PC'),]
correction = unique(DT1$Correction.x)
for(i in c(1,2,3,4)){
  UKB_D = prsD[which(prsD$Correction==correction[i]),]
  
  dt = merge(UKC_D,UKB_D[,c(2,5)],by = 'FieldIDLong.x')

  colnames(dt)[c(5,9)] = c("PRSC_UKC","PRSC_Target")
  dt$wPRSCS_UKC = (dt$PRSC_UKC*dt$QS)^2
  dt$wPRSCS_Target = (dt$PRSC_Target*dt$QS)^2
  dt$DEV = abs(dt$wPRSCS_UKC-dt$wPRSCS_Target)
  p =
    ggplot(data = dt,mapping = aes(x = wPRSCS_UKC, y = wPRSCS_Target, color=Path))+
    geom_point(size=dt$DEV*200)+
    geom_abline(slope = 1)+
    geom_hline(yintercept = mean(dt$wPRSCS_UKC),linetype=2,color = 'grey50')+
    geom_vline(xintercept = mean(dt$wPRSCS_Target),linetype=2,color = 'grey50')+
    scale_color_igv(palette = 'default')+
    scale_x_continuous(limits = c(0,0.2),breaks = seq(0,0.2,0.05),labels = round(seq(0,0.2,0.05),digits = 2))+
    scale_y_continuous(limits = c(0,0.2),breaks = seq(0,0.2,0.05),labels = round(seq(0,0.2,0.05),digits = 2))+
    theme_bw()+
    ylab(bquote(italic(R)^2 ~ " in " ~ .(correction[i]) ~ " adjusted"))+
    xlab(expression(paste(italic("R")^2," in UKC-PC adjusted")))+
    theme(plot.background = element_rect(color = 'transparent',fill = 'transparent'),
          legend.background = element_rect(color = 'transparent',fill = 'transparent'),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = 'black',linewidth = 0.5,linetype = 1),
          text = element_text(family = 'serif',size=16),
          axis.title = element_text(size = 16),
          axis.text.x = element_text(size = 13),
          axis.text.y = element_text(size = 13),
          legend.text = element_text(size = 12),
          legend.position = 'none',
          legend.title = element_blank(),
          legend.box.background = element_rect(colour = "black",fill = 'white'),
          legend.key = element_blank())
  ggsave(paste0("PRS_D_scatter", correction[i], ".pdf"),p,width = 6,height = 5.8) # Figure 5A when correction[i] = 'UKB_PC'
  
  p =
    ggplot(data = dt,mapping = aes(x = wPRSCS_UKC, y = wPRSCS_Target, color=Path))+
    geom_point(size=dt$DEV*200)+
    geom_abline(slope = 1)+
    geom_hline(yintercept = mean(dt$wPRSCS_Target),linetype=2,color = 'grey50')+
    geom_vline(xintercept = mean(dt$wPRSCS_UKC),linetype=2,color = 'grey50')+
    scale_color_igv(palette = 'default')+
    scale_x_continuous(limits = c(0,0.2),breaks = seq(0,0.2,0.05),labels = round(seq(0,0.2,0.05),digits = 2))+
    scale_y_continuous(limits = c(0,0.2),breaks = seq(0,0.2,0.05),labels = round(seq(0,0.2,0.05),digits = 2))+
    theme_bw()+
    ylab(bquote(italic(R)^2 ~ " in " ~ .(correction[i]) ~ " adjusted"))+
    xlab(expression(paste(italic("R")^2," in UKC-PC adjusted")))+
    theme(plot.background = element_rect(color = 'transparent',fill = 'transparent'),
          legend.background = element_rect(color = 'transparent',fill = 'transparent'),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = 'black',linewidth = 0.5,linetype = 1),
          text = element_text(family = 'serif',size=16),
          axis.title = element_text(size = 16),
          axis.text.x = element_text(size = 13),
          axis.text.y = element_text(size = 13),
          legend.text = element_text(size = 12),
          legend.title = element_blank(),
          legend.box.background = element_rect(colour = "black",fill = 'white'),
          legend.key = element_blank())
  ggsave(paste0("PRS_D_scatter.legend.", correction[i], ".pdf"),p,width = 6,height = 5.8)
  
  dt$DEV2 = dt$wPRSCS_UKC-dt$wPRSCS_Target
  p2 =
    ggplot(data = dt,mapping = aes(x = DEV2, fill=Path))+
    geom_density(color=NA)+
    scale_fill_igv(alpha = 0.6)+
    xlab(expression(paste(italic("R")^2," differences corrected using different schemes")))+
    ylab("Frequency")+
    theme_bw()+
    theme(plot.background = element_rect(color = 'transparent',fill = 'transparent'),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = 'black',linewidth = 0.5,linetype = 1),
          text = element_text(family = 'serif',size=16),
          axis.title = element_text(size = 16),
          axis.text.x = element_text(size = 13, angle = 315,hjust = 0.5,vjust = 0),
          axis.text.y = element_text(size = 13),
          legend.text = element_text(size = 12),
          legend.background = element_rect(color = 'transparent',fill = 'transparent'),
          legend.position = 'none',
          legend.title = element_blank(),
          legend.box.background = element_rect(colour = "black",fill = 'white'),
          legend.key = element_blank(),legend.box = element_blank())
  ggsave(paste0("PRS_D_density.", correction[i], ".pdf"),p2,width = 6,height = 6)
  
  p2 =
    ggplot(data = dt, mapping = aes(x = Path, y = DEV2, fill = Path))+
    geom_boxplot(width = 0.6)+
    scale_fill_igv(alpha = 0.9)+
    ylab(expression(paste(italic("R")^2," differences corrected using different schemes")))+
    xlab("")+
    coord_flip()+
    theme_bw()+
    geom_hline(yintercept = 0,color='grey30',linewidth = 0.5,linetype = 2)+
    theme(plot.background = element_rect(color = 'transparent',fill = 'transparent'),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = 'black',linewidth = 0.5,linetype = 1),
          text = element_text(family = 'serif',size=16,colour = 'black'),
          axis.title = element_text(size = 16),
          axis.text.x = element_text(size = 13),
          axis.text.y = element_text(size = 16,colour = 'black'),
          legend.text = element_text(size = 12),
          legend.background = element_rect(color = 'transparent',fill = 'transparent'),
          legend.position = 'none',
          legend.title = element_blank(),
          legend.box.background = element_rect(colour = "black",fill = 'white'),
          legend.key = element_blank(),legend.box = element_blank())
  ggsave(paste0("PRS_D_boxplot.", correction[i], ".pdf"),p2,width = 7,height = 6) # Figure 5B when correction[i] = 'UKB_PC'
}

# Fig 5c
library(ggplot2)
library(cowplot)
library(xlsx)
ukbPrs=read.table("./UKC.PRS.PLINKScoring.UKB_PC.txt", header = T, sep="\t")
ukcPrs=read.table("./UKC.PRS.PLINKScoring.UKC_PC.txt",header = T, sep="\t")

ukbPrs$PC = "UKB-PC"
ukcPrs$PC = "UKC-PC"
DT = rbind(ukbPrs,ukcPrs)
DT$Rsq = DT$Correlation^2
DT$Field[which(DT$Field=='Hair colour (natural, before greying)')] = 'Hair colour, natural before greying'

QSDT = plotDT2[,c(3,13)]
DT = merge(DT,QSDT,"FieldIDLong")
DT$wRsq = DT$Rsq*(DT$QS^2)
tagBody=c("21001", "50", "1747","5057", "26415","1727")
p.list = list()
for(i in 1:length(tagBody)){
  id = tagBody[i]
  dt = DT[DT$FieldID==id,]
  dt = dt[order(dt$PRange,decreasing = F),]
    
  dt$PRange = factor(dt$PRange, levels = as.character(unique(dt$PRange)))
  if(i==1){
    p =
      ggplot(data = dt, mapping = aes(x = PRange, y = wRsq, fill = PC))+
      geom_bar(stat="identity",position=position_dodge(0.9))+
      scale_fill_d3()+
      xlab("p-value cutoffs")+
      ylab(expression(paste(italic("R")^2)))+
      scale_y_continuous(limits = c(0,0.2), breaks = seq(0,0.2,0.05),labels = formatC(seq(0,0.2,0.05),digits = 2,format = 'f'))+
      labs(title = paste0(unique(dt$Field),"\n(",unique(prsD$Path[which(prsD$FieldID==id)]),")"))+
      theme_bw()+
      geom_hline(yintercept = 0)+
      theme(text = element_text(family = 'serif',size=16),
            axis.title = element_text(size = 16),
            axis.text.x = element_text(size = 12, angle = 315,hjust = 0.5,vjust = 0),
            axis.text.y = element_text(size = 12),
            legend.text = element_text(size = 12),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(colour = 'black',linewidth = 0.5,linetype = 1),
            plot.title=element_text(family = 'serif', face="bold", color="black",size=16, hjust=0.5,vjust=0.5, angle=360),
            plot.background = element_rect(color = 'transparent',fill = 'transparent'),
            legend.background = element_rect(color = 'transparent',fill = 'transparent'),
            legend.position = c(0.18,0.89),
            legend.title = element_blank(),
            legend.box.background = element_rect(colour = "transparent",fill = 'transparent'),
            legend.key = element_blank(),legend.box = element_blank())
  } else if(i==4){
    p =
      ggplot(data = dt, mapping = aes(x = PRange, y = wRsq, fill = PC))+
      geom_bar(stat="identity",position=position_dodge(0.9))+
      scale_fill_d3()+
      xlab("p-value cutoffs")+
      ylab(expression(paste(italic("R")^2)))+
      scale_y_continuous(limits = c(0,0.2), breaks = seq(0,0.2,0.05),labels = formatC(seq(0,0.2,0.05),digits = 2,format = 'f'))+
      labs(title = paste0(unique(dt$Field),"\n(",unique(prsD$Path[which(prsD$FieldID==id)]),")"))+
      theme_bw()+
      geom_hline(yintercept = 0)+
      theme(text = element_text(family = 'serif',size=16),
            axis.title = element_text(size = 16),
            axis.text.x = element_text(size = 12, angle = 315,hjust = 0.5,vjust = 0),
            axis.text.y = element_text(size = 12),
            legend.text = element_text(size = 12),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(colour = 'black',linewidth = 0.5,linetype = 1),
            plot.title=element_text(family = 'serif', face="bold", color="black",size=16, hjust=0.5,vjust=0.5, angle=360),
            plot.background = element_rect(color = 'transparent',fill = 'transparent'),
            legend.background = element_rect(color = 'transparent',fill = 'transparent'),
            legend.position = 'none',
            legend.title = element_blank(),
            legend.box.background = element_rect(colour = "transparent",fill = 'transparent'),
            legend.key = element_blank(),legend.box = element_blank())
  } else {
    p =
      ggplot(data = dt, mapping = aes(x = PRange, y = wRsq, fill = PC))+
      geom_bar(stat="identity",position=position_dodge(0.9))+
      scale_fill_d3()+
      xlab("p-value cutoffs")+
      ylab("")+
      scale_y_continuous(limits = c(0,0.2), breaks = seq(0,0.2,0.05),labels = formatC(seq(0,0.2,0.05),digits = 2,format = 'f'))+
      labs(title = paste0(unique(dt$Field),"\n(",unique(prsD$Path[which(prsD$FieldID==id)]),")"))+
      theme_bw()+
      geom_hline(yintercept = 0)+
      theme(text = element_text(family = 'serif',size=16),
            axis.title = element_text(size = 16),
            axis.text.x = element_text(size = 12, angle = 315,hjust = 0.5,vjust = 0),
            axis.text.y = element_text(size = 12),
            legend.text = element_text(size = 12),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(colour = 'black',linewidth = 0.5,linetype = 1),
            plot.title=element_text(family = 'serif', face="bold", color="black",size=16, hjust=0.5,vjust=0.5, angle=360),
            plot.background = element_rect(color = 'transparent',fill = 'transparent'),
            legend.background = element_rect(color = 'transparent',fill = 'transparent'),
            legend.position = 'none',
            legend.title = element_blank(),
            legend.box.background = element_rect(colour = "transparent",fill = 'transparent'),
            legend.key = element_blank(),legend.box = element_blank())
  }
  
  p.list[[i]] = p
}
P = plot_grid(plotlist = p.list, align = 'hv',ncol = 3)
ggsave(plot = P,filename = "Figure.5c.pdf",width = 12,height = 8)
