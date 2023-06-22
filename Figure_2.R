# Below script generates Figure 2A.
library(ggplot2)
library(RColorBrewer)
library(dplyr)
wat<-read.table("Fig2A_Data.txt", header = T, sep = "\t")
wat$Canonical.Pathways=factor(wat$Canonical.Pathways, levels=wat$Canonical.Pathways)
ggplot(wat, aes(x = Subcutaneous.Fat, y =  Canonical.Pathways)) + 
  geom_bar(stat = "identity", width=0.5, fill="#669999") + 
  xlab("-log10(adjusted p-values)") +
  ylab("Canonical Pathways") +
  theme_bw() +
  theme(axis.text = element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 3.3)) 
ggsave("Fig2A.pdf")



# Below script generates Figure 2C.
library(ggplot2)
library(RColorBrewer)
library(dplyr)
overlap<-read.table(textConnection("CPM
CTSS
FGR
GALNT16
GLCE
ITIH1
ITIH5
LEP
LGI2
LOC102115001
MMD
NEK6
RTKN
SCN7A
SLC2A1
SLC37A2
SUMF1
TMEM254"))
wat<-read.table("Subcutaneous_Fat_OBVsNonOB.txt", header = T, sep = "\t")
vat<-read.table("Visceral_Fat_OBVsNonOB.txt", header = T, sep = "\t")
bat<-read.table("Brown_Adipose_Tissue_OBVsNonOB.txt", header = T, sep = "\t")
wat_overlap<-mutate(overlap, logFC =  wat[match(overlap$V1, wat$Gene),4], CI.L =  wat[match(overlap$V1, wat$Gene),5], CI.R =  wat[match(overlap$V1, wat$Gene),6], Tissue = "Abdominal WAT")
vat_overlap<-mutate(overlap, logFC =  vat[match(overlap$V1, vat$Gene),4], CI.L =  vat[match(overlap$V1, vat$Gene),5], CI.R =  vat[match(overlap$V1, vat$Gene),6], Tissue = "Visceral WAT")
bat_overlap<-mutate(overlap, logFC =  bat[match(overlap$V1, bat$Gene),4], CI.L =  bat[match(overlap$V1, bat$Gene),5], CI.R =  bat[match(overlap$V1, bat$Gene),6], Tissue = "BAT")
overlap_combined<-rbind(wat_overlap,vat_overlap,bat_overlap)
colnames(overlap_combined)[1]="Gene"
overlap_combined<-arrange(overlap_combined, logFC)
overlap_combined$Gene=factor(overlap_combined$Gene, levels=unique(overlap_combined$Gene))
pd <- position_dodge(0.8)  
ggplot(overlap_combined, aes(y=Gene, x=logFC, col=Tissue, group=Tissue)) +
  geom_errorbar(aes(xmin=CI.L, xmax=CI.R), colour="black", width=0.7, size=0.4, position=pd) +
  geom_point(size=2.5, shape=16,  position=pd) + 
  ylab("Shared DEGs in Adipose Tissues of obesity") +
  xlab("Log2 fold change") + 
  expand_limits(x=0)+
  theme_bw() +
  scale_colour_manual(values = brewer.pal(8,"Dark2")[1:3]) +
  geom_vline(xintercept=0, linetype="dashed") +
  theme(axis.text=element_text(colour ="black"))
ggsave("Fig2C.pdf", height = 5, width = 8)



# Below script generates Figure 2D.
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)
tissue_rollup<-read.table("Tissue_Rollup_Terms.txt",head = T, sep="\t")
tissue_rollup<-arrange(tissue_rollup, Organ, Tissue.Rollup, Tissue)
tissue_rollup<-filter(tissue_rollup, Tissue!="Bone Marrow"  & Tissue!="Prostate Gland")
rownames(tissue_rollup)<-tissue_rollup$Tissue
fpkq<-read.table("QCed_primary_count_gene.FPKQ.Table.txt", header = T)
rownames(fpkq)<-fpkq$ID
head(fpkq)
dim(fpkq)
gene_list=c("LEP")
for(k in 1:length(gene_list)){
  gene=gene_list[k]
  gene_fpkq<-filter(fpkq, ID == gene)
  sample_anno<- read.table("NHP_Sample_Anno_CRO.txt",head = T, sep="\t")
  rownames(sample_anno)<-sample_anno$ID
  ObVsYc<-filter(sample_anno, Group == "NonObesity" | Group == "Obesity", Tissue!="Bone Marrow" & Tissue!="Prostate Gland", CRO == "Joinn", !str_detect(Subject,"OB4"))
  gene_fpkq<-select(gene_fpkq, ObVsYc$ID)
  ObVsYc<-mutate(ObVsYc,fpkq=t(gene_fpkq), Tissue = as.vector(Tissue), Group = as.vector(Group))
  ObVsYc<-select(ObVsYc, Tissue, Group, Subject, Sex, fpkq, CRO)
  ObVsYc$Tissue[ObVsYc$Tissue=="Subcutaneous Fat"]="Abdominal WAT"
  ObVsYc$Tissue[ObVsYc$Tissue=="Visceral Fat"]="Visceral WAT"
#tissue<-tissue_rollup$Tissue[order(tissue_rollup$Tissue)]
  tissue<-sort(unique(ObVsYc$Tissue))
  ObVsYc$Tissue=factor(ObVsYc$Tissue, levels = rev(tissue))
  pd <- position_dodge(0.8)
  p<-ggplot(ObVsYc,aes(x=Tissue, y=fpkq+0.1, col = Group, group = Group)) + geom_point(position=pd,alpha=0.5)
  for (i in 1:length(as.vector(tissue))) {
    cur_tissue_con<-filter(ObVsYc,Tissue == rev(tissue)[i], Group == "NonObesity")
    cur_tissue_di<-filter(ObVsYc,Tissue == rev(tissue)[i], Group == "Obesity")
    p<- p + geom_segment(x=i+0.2, y=log10(min(cur_tissue_con$fpkq)+0.1), xend=i+0.2, yend=log10(max(cur_tissue_con$fpkq)+0.1),col="#99CCCC",size=0.2) +
           geom_segment(x=i-0.2, y=log10(min(cur_tissue_di$fpkq)+0.1), xend=i-0.2, yend=log10(max(cur_tissue_di$fpkq)+0.1),col="#CC9752",size=0.2)
}
  fpkq_plot <- p+ geom_point(position=pd, size=2, alpha=0.7) + 
    xlab("Tissue") +
    ylab("FPKQ of Individuals") + 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(axis.text.x = element_text(colour="Black",angle = 0, hjust = 0.5, vjust=0.5,
        margin=unit(c(0.3,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(colour="Black",angle = 0, hjust = 1, vjust=0.5, margin=unit(c(0.5,0.3,0.5,0.3), "cm"))) +
    theme(legend.title=element_blank(), legend.position="top") +
    theme(axis.ticks.length=unit(-0.15, "cm")) +
    scale_y_log10(expand = c(0.03, 0), sec.axis = dup_axis(labels=NULL,name=NULL),breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", scales::math_format(10^.x)))  + scale_x_discrete(breaks = as.vector(tissue), expand = expand_scale(add = .25),labels = as.vector(tissue)) + coord_flip() + theme(axis.title.y = element_text(margin = margin(t = 0, r = -8, b = 0, l = 0))) + theme(axis.title.x = element_text(margin = margin(t = -8, r = 0, b = 0, l = 0))) +
    theme(legend.key.size = unit(0.6,"line"),legend.direction="horizontal",legend.spacing.x = unit(0.1, "cm")) + scale_color_manual(values=c("#99CCCC","#CC9752"),limits=c("Young_Control","Obesity"), labels=c("Young Control","Obesity")) + geom_hline(yintercept = 1, linetype="dashed", color = "darkred", size=0.5) + guides(color = guide_legend(order = 1),shape = guide_legend(order = 2))
## Plot log2foldchange and padj values for certain genes for Phase 3 individuals.
  tissue<-tissue_rollup$Tissue[order(tissue_rollup$Tissue)]
  fc_padj<-data.frame(tissue= as.vector(tissue), log2fc=c(rep(0,length(tissue))),     padj=c(rep(1,length(tissue))), CI.L=c(rep(0,length(tissue))), CI.R=c(rep(0,length(tissue))))
  fc_padj$tissue=factor(fc_padj$tissue, levels=rev(tissue))
  for (j in 1:length(as.vector(tissue))) {
    DEA<-read.table(paste(str_replace_all(tissue[j]," ", "_"), "_ObVsYc.txt", sep=""), header=T, sep="\t")
    if(gene %in% DEA$Gene) {
      fc_padj[j,2]= as.numeric(DEA[DEA$Gene==gene,]$logFC[1])
      fc_padj[j,3]= as.numeric(DEA[DEA$Gene==gene,]$adj.P.Val[1])
      fc_padj[j,4]= as.numeric(DEA[DEA$Gene==gene,]$CI.L[1])
      fc_padj[j,5]= as.numeric(DEA[DEA$Gene==gene,]$CI.R[1])
    }
  }
  fc_padj$tissue=as.vector(fc_padj$tissue)
  fc_padj$tissue[fc_padj$tissue=="Subcutaneous Fat"]="Abdominal WAT"
  fc_padj$tissue[fc_padj$tissue=="Visceral Fat"]="Visceral WAT"
  tissue<-sort(unique(ObVsYc$Tissue))
  fc_padj$tissue=factor(fc_padj$tissue, levels = tissue)
  DEGlabel <- rep("Non-DEG", nrow(fc_padj))
  DEGlabel[ (fc_padj$log2fc > log2(1.5) | fc_padj$log2fc < -log2(1.5)) & fc_padj$padj < 0.05 ]  <- "DEG"
  fc_padj<-mutate(fc_padj, DEG = DEGlabel)
  font <- rep("plain", nrow(fc_padj))
  font[ (fc_padj$log2fc > log2(1.5) | fc_padj$log2fc < -log2(1.5)) & fc_padj$padj < 0.05 ]  <- "bold"
  fc_padj<-mutate(fc_padj, font = font)
  fc_padj<-arrange(fc_padj, tissue)
  font=fc_padj$font
  if(length(DEGlabel[DEGlabel=="DEG"])<1){
    log2fc_plot<-ggplot(fc_padj,aes(x=tissue,y=log2fc,fill=DEG)) + geom_bar(stat ="identity") +
    ylab("Log2 Fold Change") + 
    geom_hline(yintercept = c(-log2(1.5),0,log2(1.5)),col= c("darkred","black","darkred"),linetype=c(2,1,2)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(axis.text.x = element_text(colour="Black",angle = 90, hjust = 1, vjust=0.5,
        margin=unit(c(0.3,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(colour="Black",margin=unit(c(0.5,0.3,0.5,0.3), "cm"))) +
    theme(legend.title=element_blank(), legend.position="top") +
    theme(axis.ticks.length=unit(-0.15, "cm")) +
    scale_y_continuous(expand = expand_scale(add = .25), sec.axis = dup_axis(labels=NULL,name=NULL)) +
    scale_x_discrete(breaks = tissue, expand = expand_scale(add = .2),labels = tissue) +
    coord_flip() +  
    theme(axis.text.y = element_text(face = font)) + 
    theme(axis.title.y = element_text(margin = margin(t = 0, r = -8, b = 0, l = 0))) +
    theme(axis.title.x = element_text(margin = margin(t = -8, r = 0, b = 0, l = 0))) +
    theme(axis.title.y = element_blank()) +
    theme(legend.key.height = unit(0.5, "line"), legend.key.width = unit(2, "line"), legend.direction="horizontal") + scale_fill_manual(breaks=c("Non-DEG"), values = c("lightgrey"))

  padj_plot <- ggplot(fc_padj, aes(x=tissue, y=-log10(padj),fill=DEG)) + 
    geom_bar(position=pd,stat="identity") +
    ylab("-Log10(adj p-value)") + 
    geom_hline(yintercept = c(-log10(0.05),0),col= c("darkred","black"),linetype=c(2,1)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(axis.text.x = element_text(face = rev(fc_padj$font), colour="Black",angle = 90, hjust = 1, vjust=0.5,
        margin=unit(c(0.3,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(colour="Black",margin=unit(c(0.5,0.3,0.5,0.3), "cm"))) +
    theme(legend.title=element_blank(), legend.position="top") +
    theme(axis.ticks.length=unit(-0.15, "cm")) +
    scale_y_continuous(expand = expand_scale(add = .25), sec.axis = dup_axis(labels=NULL,name=NULL)) +
    scale_x_discrete(breaks = tissue, expand = expand_scale(add = .2),labels = tissue) + 
    coord_flip() +  theme(axis.text.y = element_text(face = font)) +
    theme(axis.title.x = element_text(margin = margin(t = -8, r = 0, b = 0, l = 0))) +
    theme(legend.key.height = unit(0.5, "line"), legend.key.width = unit(2, "line") ,legend.direction="horizontal") + theme(axis.title.y = element_blank()) + scale_fill_manual(breaks=c("Non-DEG"), values = c("lightgrey"))
  } 
  else {
    log2fc_plot<-ggplot(fc_padj,aes(x=tissue,y=log2fc,fill=DEG)) + geom_bar(stat="identity",size=0.7) +
      ylab("Log2 Fold Change") + 
      geom_hline(yintercept = c(-log2(1.5),0,log2(1.5)),col= c("darkred","black","darkred"),linetype=c(2,1,2)) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
      theme(axis.text.x = element_text(colour="Black",angle = 90, hjust = 1, vjust=0.5,
        margin=unit(c(0.3,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(face =  rev(font), colour="Black",margin=unit(c(0.5,0.3,0.5,0.3), "cm"))) +
      theme(legend.title=element_blank(), legend.position="top") +
      theme(axis.ticks.length=unit(-0.15, "cm")) +
      scale_y_continuous(expand = expand_scale(add = .25), sec.axis = dup_axis(labels=NULL,name=NULL)) +
      scale_x_discrete(breaks = tissue, expand = expand_scale(add = .2),labels = tissue) + coord_flip() +  theme(axis.text.y = element_text(face = font)) + theme(axis.title.y = element_text(margin = margin(t = 0, r = -8, b = 0, l = 0))) + theme(axis.title.x = element_text(margin = margin(t = -8, r = 0, b = 0, l = 0))) + theme(axis.title.y = element_blank()) +
      theme(legend.key.height = unit(0.5, "line"), legend.key.width = unit(2, "line"), legend.direction="horizontal") + scale_fill_manual(breaks=c("DEG","Non-DEG"), values = c("#FB9A99", "lightgrey"))
#log2fc_plot
  padj_plot <- ggplot(fc_padj, aes(x=tissue, y=-log10(padj),fill=DEG)) + 
    geom_bar(position=pd,stat="identity",size=0.7) +
    ylab("-Log10(adj p-value)") + 
    geom_hline(yintercept = c(-log10(0.05),0),col= c("darkred","black"),linetype=c(2,1)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(axis.text.x = element_text(colour="Black",angle = 90, hjust = 1, vjust=0.5,
        margin=unit(c(0.3,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(face = rev(font), colour="Black",margin=unit(c(0.5,0.3,0.5,0.3), "cm"))) +
    theme(legend.title=element_blank(), legend.position="top") +
    theme(axis.ticks.length=unit(-0.15, "cm")) +
    scale_y_continuous(expand = expand_scale(add = .25), sec.axis = dup_axis(labels=NULL,name=NULL)) +
    scale_x_discrete(breaks = tissue, expand = expand_scale(add = .2),labels = tissue) + 
    coord_flip() +  theme(axis.text.y = element_text(face = font)) +
    theme(axis.title.x = element_text(margin = margin(t = -8, r = 0, b = 0, l = 0))) +
    theme(legend.key.height = unit(0.5, "line"), legend.key.width = unit(2, "line") ,legend.direction="horizontal") + theme(axis.title.y = element_blank()) + scale_fill_manual(breaks=c("DEG","Non-DEG"), values = c("#FB9A99", "lightgrey"))
  }
#padj_plot
  library(ggpubr)
  figure<-ggarrange(fpkq_plot, log2fc_plot, padj_plot, widths = c(1.5,1,1),
          labels = c("", "", ""),
          ncol = 3, nrow = 1, align = "h")
  annotate_figure(figure, top = text_grob(paste("Expression of", gene, sep=" "), color = "black", face = "bold", size = 14))
  ggsave(paste(gene, "_Fig2D.pdf", sep=""), width=12, height=5)
}



# Below script generates Figure 2E.
library(forestplot)
all<-read.table("Fig2E_Data.txt",header=TRUE, sep="\t")
pdf("Fig2E.pdf",width=7,height=9)
forestplot(as.vector(all$Cell_Type), subset(all,select=c("OR","X2.50.","X97.50.")),xlog = TRUE, xlab="OR of DEGs Enrichment", col=fpColors(box = "steelblue",line = "steelblue",zero = "gray"),zero = 1)
dev.off()



