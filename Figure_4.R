# Below script generates Figure 4A.
library(dplyr)
library(stringr)
library(ggplot2)
library("RColorBrewer")
protein_gene_ncbi<-read.table("Protein_Gene_Match.txt", sep = "\t", header = T)
tissue=c("Subcutaneous_Fat")
pdf("Fig4A.pdf")
for(i in 1:length(tissue)){
  genes<-read.table(paste(tissue[i],"_OBVsNonOB_RNA.txt", sep=""), sep = "\t", header = T)
  proteins<-read.table(paste( tissue[i], "_OBVsNonOB_Protein.txt", sep=""), sep="\t", header=T)
  proteins<-mutate(proteins, protein_ID = matrix(unlist(str_split(proteins$Accession, "\\.")), ncol=2, byrow=TRUE)[,1])
  proteins<-mutate(proteins, gene_ncbi = protein_gene_ncbi[match(proteins$protein_ID,protein_gene_ncbi$Protein),][,2])
  all_overlap<-merge.data.frame(genes, proteins, by.x = "Gene", by.y="gene_ncbi")
  DEGs<-read.table(paste(tissue[i],"_SigDEGsOBVsNonOB.txt", sep=""), sep = "\t", header = T)
  DEPs<-read.table(paste(tissue[i], "_SigDEPsOBVsNonOB.txt", sep=""), sep="\t", header=T)
  if(dim(DEPs)[1]>0){
    DEPs<-mutate(DEPs, protein_ID = matrix(unlist(str_split(DEPs$Accession, "\\.")), ncol=2, byrow=TRUE)[,1])
  }
  DEPs<-mutate(DEPs, gene_ncbi = protein_gene_ncbi[match(DEPs$protein_ID,protein_gene_ncbi$Protein),][,2])
gene_to_plot<-unique(c(as.vector(DEGs$Gene),as.vector(DEPs$gene_ncbi)))
  gene_to_plot_mat<-filter(all_overlap, Gene %in% gene_to_plot)
  gene_fre<-c(as.vector(DEGs$Gene),as.vector(DEPs$gene_ncbi))
  gene_fre_data<-as.data.frame(table(gene_fre))
  DEGsDEPs<-filter(gene_fre_data, Freq == 2)
  DEGsDEPs<-mutate(DEGsDEPs, Group = "DEGs & DEPs")
  DEGsOnly<-filter(gene_fre_data, Freq == 1 & gene_fre %in% DEGs$Gene)
  DEGsOnly<-mutate(DEGsOnly, Group = "DEGs")
  DEPsOnly<-filter(gene_fre_data, Freq == 1 & gene_fre %in% DEPs$gene_ncbi)
  DEPsOnly<-mutate(DEPsOnly, Group = "DEPs")
  group_info<-rbind(DEGsDEPs,DEGsOnly,DEPsOnly)
  if(dim(gene_to_plot_mat)[1]>0){
    gene_to_plot_mat<-mutate(gene_to_plot_mat, Group = group_info[match(gene_to_plot_mat$Gene,group_info$gene_fre),][,3])
    gene_to_plot_mat<-filter(gene_to_plot_mat, Group == "DEGs & DEPs" & abs(logFC.x) > 0.585 & abs(logFC.y) > 0.585)
    plots<-ggplot(gene_to_plot_mat, aes(x =  logFC.x, y = logFC.y)) + geom_vline(xintercept = c(-0.585,0,0.585), linetype=c(2,1,2), color = c("black","grey","black"), size=0.3)  + geom_hline(yintercept = c(-0.585,0,0.585), linetype=c(2,1,2), color = c("black","grey","black"), size=0.3)  + geom_point(aes(col=Group),alpha=0.5) + coord_fixed()  + xlab("Log2 Fold Change of RNA") + ylab("Log2 Fold Change of Protein") + labs(title = str_replace_all(tissue[i], "_", " ")) + geom_smooth(method='lm', formula= y~x, size =0.5, col="darkred") + stat_cor(label.x = -3, label.y = 3) + stat_regline_equation(label.x = -3, label.y = 2.5) + xlim(-3,3) + ylim(-3,3) + theme(plot.title = element_text(hjust = 0.5)) + scale_color_manual(values=brewer.pal(3,"Dark2")[2], breaks = c("DEGs & DEPs"), labels = c("DEGs & DEPs"))
    print(plots)
  }
}
dev.off()



# Below script generates Figure 4B.
library(ggplot2)
library(RColorBrewer)
library(dplyr)
chrx<-read.table("Fig4B_Data.txt",header=T)
ax <- length(unique(chrx$Gene)) 
pd <- position_dodge(0.4) 
ggplot(chrx, aes(x=Gene, y=logfc,col=DataType, group=DataType))  + 
    geom_errorbar(aes(ymin=CIL, ymax=CIR, col=DataType), width=0, size=0.2, position=pd,alpha=0.8) + geom_polygon(fill=NA, color = "grey",size=0.5) +
    geom_point(size=1, position=pd, alpha=0.8) +  coord_polar(theta = "x") + 
    xlab("Overlapping DEGs and DEPs") +
    ylab("Log2 fold change") +
    expand_limits(y=0)    +
    theme_bw()  +
    theme(axis.line.y = element_line(),panel.border = element_blank(),axis.text.y = element_text(size=9), axis.text.x = element_text(size=9,vjust=1,hjust=1,colour="Black", angle = 360/(2*pi)*rev( pi/2 + seq( pi/ax, 2*pi-pi/ax, len=ax) + c( rep(135, as.integer(ax/2)),rep(0,ax-as.integer(ax/2))))) )    + 
    scale_x_discrete(limits=chrx$Gene[1:ax], labels=c(paste0(rep("",as.integer(ax/2)),chrx$Gene[1:as.integer(ax/2)]),paste0(chrx$Gene[as.integer(ax/2)+1:ax],rep("",ax-as.integer(ax/2))))) + geom_hline(yintercept=0, linetype="dashed", color = "black", size = 0.6) + scale_y_continuous(expand = c(0, 0),breaks=c(-4, -2,0,2,4),labels=c(4, -2,0,2,4),limits=c(-5.8,11.5)) +   
scale_colour_manual(name = "Omics Type",
                     labels = c("DEGs", "DEPs"),
                      values = c("#7d9fc2","#C582B2"))
ggsave("Fig4B.pdf")



# Below script generates Figure 4C.
library(ggplot2)
library(RColorBrewer)
library(dplyr)
chrx<-read.table("Fig4C_Data.txt",header=T)
ax <- length(unique(chrx$Gene)) 
pd <- position_dodge(0.4) 
ggplot(chrx, aes(x=Gene, y=logfc,col=DataType, group=DataType))  + 
    geom_errorbar(aes(ymin=CIL, ymax=CIR, col=DataType), width=0, size=0.3, position=pd,alpha=1) + geom_polygon(fill=NA, color = "grey",size=0.5) +
    geom_point(size=1.5, position=pd, alpha=1) +  coord_polar(theta = "x") + 
    xlab("Overlapping DEGs and DEPs") +
    ylab("Log2 fold change") +
    expand_limits(y=0)    +
    theme_bw()  +
    theme(axis.line.y = element_line(),panel.border = element_blank(),axis.text.y = element_text(size=9), axis.text.x = element_text(size=9,vjust=1,hjust=1,colour="Black", angle = 360/(2*pi)*rev( pi/2 + seq( pi/ax, 2*pi-pi/ax, len=ax) + c( rep(135, as.integer(ax/2)),rep(0,ax-as.integer(ax/2))))) )    + 
    scale_x_discrete(limits=chrx$Gene[1:ax],
    labels=c(paste0(rep("",as.integer(ax/2)),chrx$Gene[1:as.integer(ax/2)]),paste0(chrx$Gene[as.integer(ax/2)+1:ax],rep("",ax-as.integer(ax/2))))) + geom_hline(yintercept=0, linetype="dashed", color = "black", size = 0.6) + scale_y_continuous(expand = c(0, 0),breaks=c(-4, -2,0,2,4),labels=c(4, -2,0,2,4),limits=c(-10.4,4.9)) +   
scale_colour_manual(name = "Omics Type",
                     labels = c("DEGs", "DEPs"),
                      values = c("#7d9fc2","#C582B2"))
ggsave("Fig4C.pdf")



# Below script generates Figure 4D-E.
library(dplyr)
library(stringr)
library(circlize)
library(pals)
gene_anno<- read.table("feature.Gene.anno",head = T, sep="\t")
rownames(gene_anno)<-gene_anno$Id
protein_gene_ncbi<-read.table("Protein_Gene_Match.txt", sep = "\t", header = T)
tissue=c("Subcutaneous_Fat","Visceral_Fat")
for(i in 1:length(tissue)){
  DEGs<-read.table(paste(tissue[i],"_OBVsNonOB_DEGs.txt", sep=""), sep = "\t", header = T)
  DEPs<-read.table(paste(tissue[i], "_OBVsNonOB_DEPs.txt", sep=""), sep="\t", header=T)
  if (dim(DEPs)[1]==0 | dim(DEGs)[1]==0) {
  next
  }
  if(dim(DEPs)[1]>0){
    DEPs<-mutate(DEPs, protein_ID = matrix(unlist(str_split(DEPs$Accession, "\\.")), ncol=2, byrow=TRUE)[,1])
  }
  DEPs<-mutate(DEPs, gene_ncbi = protein_gene_ncbi[match(DEPs$protein_ID,protein_gene_ncbi$Protein),][,2])
  DEPs<-mutate(DEPs, Chr=gene_anno[match(DEPs$gene_ncbi,gene_anno$Id),][,9], Start=gene_anno[match(DEPs$gene_ncbi,gene_anno$Id),][,11], End=gene_anno[match(DEPs$gene_ncbi,gene_anno$Id),][,12],)
  DEGs_circos<-select(DEGs, Chr, Start, End, logFC, CI.L, CI.R)
  DEPs_circos<-select(DEPs, Chr, Start, End, logFC, CI.L, CI.R)
  DEGs_circos$Chr = paste0("chr", DEGs_circos$Chr)
  DEPs_circos$Chr = paste0("chr", DEPs_circos$Chr)
  DEGs_circos<-arrange(DEGs_circos, Chr, Start)
  DEPs_circos<-arrange(DEPs_circos, Chr, Start)
  pdf(paste(tissue[i], "_circos.pdf", sep=""))
  circos.par("start.degree" = 90)
  circos.initializeWithIdeogram(species = "macFas5")
  circos.genomicTrack(ylim = c(-4, 4), DEGs_circos[,1:4], numeric.column = 4, bg.col = c(kelly()[c(-2,-22)],"#2BCE48"),
    panel.fun = function(region, value, ...) {
        circos.segments(0, 0, CELL_META$xlim[2], 0, col = "#FFFFFF", lwd = 1, lty = 1)
        circos.genomicLines(region, value, type = 'h', baseline = 0, col = ifelse(value>0, "#565656", "#565656"))
    })
  circos.yaxis("right", sector.index = get.current.sector.index(), labels.cex = 0.5)
  circos.genomicTrack(ylim = c(-4, 4), DEPs_circos[,1:4], numeric.column = 4, bg.col = c(kelly()[c(-2,-22)],"#2BCE48"),
    panel.fun = function(region, value, ...) {
        circos.segments(0, 0, CELL_META$xlim[2], 0, col = "#FFFFFF", lwd = 1, lty = 1)
        circos.genomicLines(region, value, type = 'h', baseline = 0, col = ifelse(value>0, "#565656", "#565656"))
  })
  circos.yaxis("right", sector.index = get.current.sector.index(), labels.cex = 0.5)
  circos.genomicDensity(DEGs_circos[,1:4], window.size = 5e6, col = c("#99CC99"), track.height = 0.12)
  circos.yaxis("right", at = c(0, 0.1, 0.2), sector.index = get.current.sector.index(), labels.cex = 0.5)
  circos.genomicDensity(DEPs_circos[,1:4], window.size = 5e6, col = c("#993300"), track.height = 0.12)
  circos.yaxis("right", at = c(0, 0.04, 0.08), sector.index = get.current.sector.index(), labels.cex = 0.5)
  circos.clear()
  dev.off()
}



# Below script generates Figure 4F-H.
library(limma)
library(dplyr)
library(edgeR)
library(stringr)
library(ggplot2)
library(mixOmics)
tissue=c("Subcutaneous_Fat")
for(i in 1:length(tissue)){
  protein_count<-read.table(paste(tissue[i], "_OBVsNonOB_Norm_Count.txt", sep=""), sep="\t", header=T)
  gene_count<-read.table(paste(tissue[i], "_OBVsNonOB_Norm_Count.txt", sep=""), sep="\t", header=T)
  gene_count2<-gene_count
  protein_count1<-dplyr::select(protein_count, Gene, colnames(gene_count2))
  protein_count2<- protein_count1 %>% group_by(Gene) %>% summarise(across(.cols = where(is.numeric), .fns = list(Sum = mean), na.rm = TRUE, .names = "{col}"))
  protein_count3<-as.data.frame(protein_count2)
  rownames(protein_count3)<-protein_count3$Gene
  protein<-t(protein_count3[,-1])
  mrna<-t(gene_count2)
  pls.res = pls(protein, mrna, ncomp = 1)
  corre<-cor(pls.res$variates$X, pls.res$variates$Y)
  Y<-c(rep("YC",table(str_detect(colnames(gene_count2),"YC|OW2"))[2]),rep("OB",table(str_detect(colnames(gene_count2),"OB|PD|OW1"))[2]))
  data = list(mRNA = mrna, proteomics = protein)
  design = matrix(corre[1,1], ncol = length(data), nrow = length(data), dimnames = list(names(data), names(data)))
  diag(design) = 0
  ncomp = 2
  list.keepX = list(mRNA = c(5,5), proteomics = c(6,30))
  sgccda.res = block.splsda(X = data, Y = Y, ncomp = ncomp, keepX = list.keepX, design = design)
  sgccda.res$design
  pdf(paste(tissue[i], "_Fig4F.pdf", sep=""))
  cimDiablo(sgccda.res, size.legend = 0.6, margins = c(3, 6), legend.position = c(0.8,0.7), transpose = TRUE)
  cimDiablo(sgccda.res, size.legend = 0.6, margins = c(3, 15))
  dev.off()  
  varDat.train <- do.call(rbind, sgccda.res$variates[1:length(data)])
  varDat.train <- as.data.frame(varDat.train) %>% mutate(subj = rownames(varDat.train))
  varDat.train$Dataset <- rep(names(data), e = length(Y))
  varDat.train$Class <- Y
  varDat.train <- varDat.train %>% group_by(Class, subj) %>% 
    summarise_all(funs(mean)) %>% 
    mutate(Dataset = "Consensus") %>% 
    dplyr::select(`comp1`, `comp2`, subj, Dataset, Class) %>% 
    as.data.frame() %>% 
    rbind(., varDat.train) %>% 
    mutate(Dataset = factor(Dataset, c("Consensus", "mRNA", "Proteins"))) 
  filter(varDat.train, Dataset == "Consensus") %>% 
    ggplot(aes(x = `comp1`, y = `comp2`, color = Class)) + 
    geom_point(size = 2) +
    facet_wrap(~Dataset, scales = "free", ncol = 5) + 
   stat_ellipse(data = filter(varDat.train, Dataset == "Consensus"), size = 1) +
    xlab("Component 1") + ylab("Component 2") + 
    theme(strip.text.x = element_text(size=26, face = "bold")) + 
    scale_color_manual(values=color.mixo(1:4))
  ggsave(paste(tissue[i], "_Fig4G.pdf", sep=""))  
  list.keepX = list(mRNA = c(5), proteomics = c(6))
  ncomp=1
  protein_count4<-protein_count3
  edit_protein_names<-paste0("P_",rownames(protein_count4))
  rownames(protein_count4)<-edit_protein_names
  protein<-t(protein_count4[,-1])
  data = list(mRNA = mrna, proteomics = protein)
  design = matrix(corre[1,1], ncol = length(data), nrow = length(data), dimnames = list(names(data), names(data)))
  diag(design) = 0
  sgccda.res = block.splsda(X = data, Y = Y, ncomp = ncomp, keepX = list.keepX, design = design)
  library(igraph)
  my.network = network(sgccda.res, blocks = c(1,2), color.node = c('darkorchid', 'brown1'), cutoff = 0.4, save = 'pdf') 
  write.graph(my.network$gR, file = paste(cur_dir, tissue[i], "_myNetwork.gml", sep=""), format = "gml")
  corMat.P <- circosPlot(sgccda.res, cutoff = 0.7, ncol.legend = 2, size.legend = 0.8)
  corrplot::corrplot(corMat.P, order="hclust", method = 'circle', col=brewer.pal(n=8, name="RdYlBu"), addrect = 2 ,tl.cex = 0.55, type = 'lower', col = COL2('BrBG'))
  pdf("Fig4H.pdf")
  corrplot::corrplot(corMat.P, order="hclust", method = 'circle', col=colorRampPalette(brewer.pal(3, "PuOr"))(50), addrect = 2, tl.cex = 0.55, type = 'lower')
  dev.off()
}








