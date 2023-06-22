# Below script generates Figure 5A.
library(ggplot2)
library(scatterpie)
library(dplyr)
library(stringr)
library(RColorBrewer)
DEG<-read.table("SameDir_DEG_DEPs_T2D.txt", header = T, sep = "\t")
DEG<-mutate(DEG, Disease_code=rep(1,21))
DEG<-mutate(DEG, Tissue_index=1:21)
DEG<-mutate(DEG, Fake=rep(0,21), Fake1=SameDir)
DEG_plot<-select(DEG, Tissue_index, Disease_code, SameDir, Fake, Fake1)
leg<-0.04*log2(as.vector(na.omit(DEG_plot$SameDir+1.1)))
leg<-as.numeric(leg)
DEG_NA<-filter(DEG_plot, is.na(DEG))
ggplot() + geom_scatterpie(aes(x=Tissue_index, y=Disease_code, group = Disease_code, r = 0.04 * log2(SameDir)), data = DEG_plot, cols = c("Fake1","Fake")) + geom_scatterpie_legend(leg, x=6, y=1, n=4, labeller=function(x) as.integer(2^(x/0.04))) + theme_bw()  + xlab("Tissues") + ylab("Disease") + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust =1, color = "black"), axis.text.y = element_text(color = "black"), legend.position = c(0.7,0.88), legend.title = element_blank(), legend.direction="horizontal", legend.key.size = unit(0.5, "cm"), )  +  theme(panel.grid.minor = element_blank(), legend.background = element_blank()) + scale_fill_manual(values = c( "#CCFFBB", "#005A04"), breaks = c("Fake1","Fake"), labels = c("Upregulated DEGs","Downregulated DEGs")) + coord_polar() + scale_x_continuous(limits=c(0.5,21.5), breaks = 1:21, labels = paste(DEG$Tissue,DEG$SameDir), expand = expand_scale(mult = c(0.01, 0.01))) + theme(axis.line.y = element_line(),panel.border = element_blank(),axis.text.y = element_text(size=9), axis.text.x = element_text(size=9,vjust=1,hjust=1,colour="Black", angle = 360/(2*pi)*rev( pi/2 + seq( pi/21, 2*pi-pi/21, len=21) + c( rep(135, as.integer(21/2)),rep(0,21-as.integer(21/2))))) ) 
ggsave("Fig5A.pdf")



# Below script generates Figure 5B.
library(limma)
library(dplyr)
library(edgeR)
library(stringr)
library(ggplot2)
library(mixOmics)
tissue=c("Pancreas")
for(i in 1:length(tissue)){
  protein_count<-read.table(paste("_T2DVsNonT2D_Norm_Count.txt", sep=""), sep="\t", header=T)
  gene_count<-read.table(paste(tissue[i], "_T2DVsNonT2D_Norm_Count.txt", sep=""), sep="\t", header=T)
  gene_count2<-gene_count
  protein_count1<-dplyr::select(protein_count, Gene, colnames(gene_count2))
  protein_count2<- protein_count1 %>% group_by(Gene) %>% summarise(across(.cols = is.numeric, .fns = list(Sum = mean), na.rm = TRUE, .names = "{col}"))
  protein_count3<-as.data.frame(protein_count2)
  rownames(protein_count3)<-protein_count3$Gene
  protein<-t(protein_count3[,-1])
  mrna<-t(gene_count2)
  pls.res = pls(protein, mrna, ncomp = 1)
  cor(pls.res$variates$X, pls.res$variates$Y)
  corre<-cor(pls.res$variates$X, pls.res$variates$Y)
  Y<-c(rep("NonT2D",table(str_detect(colnames(gene_count2),"NonT2D"))[2]),rep("T2D",table(str_detect(colnames(gene_count2),"T2D"))[2]))
  data = list(mRNA = mrna, proteomics = protein)
  design = matrix(corre[1,1], ncol = length(data), nrow = length(data), dimnames = list(names(data), names(data)))
  diag(design) = 0
  sgccda.res = block.splsda(X = data, Y = Y, ncomp = 5, design = design)
  set.seed(123)
  perf.diablo = perf(sgccda.res, validation = 'loo', nrepeat = 10)
  ncomp = 2
  list.keepX = list(mRNA = c(25,25), proteomics = c(5,5))
  pdf(paste(tissue[i], "_Fig5B.pdf", sep=""))
  cimDiablo(sgccda.res, size.legend = 0.6, margins = c(3, 6), legend.position = c(0.8,0.7), transpose = TRUE)
  cimDiablo(sgccda.res, size.legend = 0.6, margins = c(3, 15))
  dev.off()
}



# Below script generates Figure 5D.
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)
tissue_rollup<-read.table("Tissue_Rollup_Terms.txt",head = T, sep="\t")
tissue_rollup<-arrange(tissue_rollup, Organ, Tissue.Rollup, Tissue)
tissue_rollup<-filter(tissue_rollup, Tissue!="Uterus" & Tissue!="Brown Adipose Tissue" & Tissue!="Prostate Gland")
rownames(tissue_rollup)<-tissue_rollup$Tissue
fpkq<-read.table("QCed_primary_count_gene.FPKQ.Table.txt", header = T)
rownames(fpkq)<-fpkq$ID
gene_list=c("SLC30A8")
for(k in 1:length(gene_list)){
gene=gene_list[k]
gene_fpkq<-filter(fpkq, ID == gene)
sample_anno<- read.table("NHP_Sample_Anno_CRO.txt",head = T, sep="\t")
rownames(sample_anno)<-sample_anno$ID
DivsAC<-filter(sample_anno, Group == "Aged_Control" | Group == "Diabetes", Tissue!="Uterus" & Tissue!="Brown Adipose Tissue" & Tissue!="Prostate Gland", CRO == "Huazhen_Wuxi")
gene_fpkq<-select(gene_fpkq, DivsAC$ID)
DivsAC<-mutate(DivsAC,fpkq=t(gene_fpkq), Tissue = as.vector(Tissue), Group = as.vector(Group))
DivsAC<-select(DivsAC, Tissue, Group, Subject, Sex, fpkq, CRO)
DivsAC$Tissue[DivsAC$Tissue=="Subcutaneous Fat"]="Abdominal WAT"
DivsAC$Tissue[DivsAC$Tissue=="Visceral Fat"]="Visceral WAT"
tissue<-sort(unique(DivsAC$Tissue))
DivsAC$Tissue=factor(DivsAC$Tissue, levels = rev(tissue))
DivsAC$Group[DivsAC$Group=="Aged_Control"]="X_Aged_Control"
pd <- position_dodge(0.8)
p<-ggplot(DivsAC,aes(x=Tissue, y=fpkq+0.1, col = Group, group = Group)) + geom_point(position=pd,alpha=0.5)
for (i in 1:length(as.vector(tissue))) {
  cur_tissue_con<-filter(DivsAC,Tissue == rev(tissue)[i], Group == "X_Aged_Control")
  cur_tissue_di<-filter(DivsAC,Tissue == rev(tissue)[i], Group == "Diabetes")
   p<- p + geom_segment(x=i+0.2, y=log10(min(cur_tissue_con$fpkq)+0.1), xend=i+0.2, yend=log10(max(cur_tissue_con$fpkq)+0.1),col="#CABEE9",size=0.2) +
           geom_segment(x=i-0.2, y=log10(min(cur_tissue_di$fpkq)+0.1), xend=i-0.2, yend=log10(max(cur_tissue_di$fpkq)+0.1),col="#BC8E7D",size=0.2)
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
    theme(legend.key.size = unit(0.6,"line"),legend.direction="horizontal",legend.spacing.x = unit(0.1, "cm")) + scale_color_manual(values=c("#CABEE9","#BC8E7D"),limits=c("X_Aged_Control","Diabetes"), labels=c("Aged Control","Diabetes")) + geom_hline(yintercept = 1, linetype="dashed", color = "darkred", size=0.5) + guides(color = guide_legend(order = 1),shape = guide_legend(order = 2))
fpkq_plot
## Plot log2foldchange and padj values for certain genes for Phase 3 individuals.
tissue<-tissue_rollup$Tissue[order(tissue_rollup$Tissue)]
fc_padj<-data.frame(tissue= as.vector(tissue), log2fc=c(rep(0,length(tissue))), padj=c(rep(1,length(tissue))), CI.L=c(rep(0,length(tissue))), CI.R=c(rep(0,length(tissue))))
fc_padj$tissue=factor(fc_padj$tissue, levels=rev(tissue))
for (j in 1:length(as.vector(tissue))) {
  DEA<-read.table(paste(str_replace_all(tissue[j]," ", "_"), "_DivsAC_DM_vs_AC.txt", sep=""), header=T, sep="\t")
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
tissue<-sort(unique(DivsAC$Tissue))
fc_padj$tissue=factor(fc_padj$tissue, levels = tissue)
DEGlabel <- rep("Non-DEG", nrow(fc_padj))
DEGlabel[ (fc_padj$log2fc > log2(1.5) | fc_padj$log2fc < -log2(1.5)) & fc_padj$padj < 0.01 ]  <- "DEG"
fc_padj<-mutate(fc_padj, DEG = DEGlabel)
font <- rep("plain", nrow(fc_padj))
font[ (fc_padj$log2fc > log2(1.5) | fc_padj$log2fc < -log2(1.5)) & fc_padj$padj < 0.01 ]  <- "bold"
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
  scale_x_discrete(breaks = tissue, expand = expand_scale(add = .2),labels = tissue) + coord_flip() +  theme(axis.text.y = element_text(face = font)) + theme(axis.title.y = element_text(margin = margin(t = 0, r = -8, b = 0, l = 0))) + theme(axis.title.x = element_text(margin = margin(t = -8, r = 0, b = 0, l = 0))) + theme(axis.title.y = element_blank()) +
    theme(legend.key.height = unit(0.5, "line"), legend.key.width = unit(2, "line"), legend.direction="horizontal") + scale_fill_manual(breaks=c("Non-DEG"), values = c("lightgrey"))
#log2fc_plot
padj_plot <- ggplot(fc_padj, aes(x=tissue, y=-log10(padj),fill=DEG)) + 
    geom_bar(position=pd,stat="identity") +
  ylab("-Log10(adj p-value)") + 
  geom_hline(yintercept = c(-log10(0.01),0),col= c("darkred","black"),linetype=c(2,1)) +
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
} else {
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
  geom_hline(yintercept = c(-log10(0.01),0),col= c("darkred","black"),linetype=c(2,1)) +
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
library(ggpubr)
figure<-ggarrange(fpkq_plot, log2fc_plot, padj_plot, widths = c(1.5,1,1),
          labels = c("", "", ""),
          ncol = 3, nrow = 1, align = "h")
annotate_figure(figure, top = text_grob(paste("Expression of", gene, sep=" "), color = "black", face = "bold", size = 14))
ggsave("Fig5D.pdf", width=12, height=5)
}
