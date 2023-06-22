# Below script describes the workflow of differential RNA expression analysis between obese vs. non-obese monkeys for all tissues.
# Load libraries.
library(limma)
library(dplyr)
library(edgeR)
library(stringr)
library(ggplot2)
# Read gene annotation info.
gene_anno<- read.table("feature.Gene.anno",head = T, sep="\t")
rownames(gene_anno)<-gene_anno$Id
head(gene_anno)
# Read sample annotation info.
sample_anno<- read.table("NHP_Sample_Anno_CRO.txt",head = T, sep="\t")
rownames(sample_anno)<-sample_anno$ID
head(sample_anno)
tissue<-sort(unique(sample_anno$Tissue))
for(i in c(1:2,4:17,19:27)){
  tissue_sample<-read.table(paste(tissue[i], "_OBVsNonOB_Sample_Anno.txt", sep=""), sep="\t", header=T)
  tissue_count<-read.table(paste(tissue[i], "_OBVsNonOB_Raw_Count.txt", sep=""), sep="\t", header=T)
  fpkq<-read.table(paste(tissue[i], "_OBVsNonOB_FPKQ.txt", sep=""),sep="\t",header=T)
  row.names(tissue_count)=tissue_count$ID
  table(tissue_sample$ID==colnames(tissue_count[,-1]))
  if(i==26){
    y <- DGEList(counts = tissue_count[,-1])
    y <- calcNormFactors(y)
    design <- model.matrix(~ Disease, data = tissue_sample)
    v <- voom(y, design, plot=TRUE, normalize="quantile")
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
  } 
  else {
    y <- DGEList(counts = tissue_count[,-1])
    y <- calcNormFactors(y)
    design <- model.matrix(~ Disease + Sex, data = tissue_sample)
    v <- voom(y, design, plot=TRUE, normalize="quantile")
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
  }
  summary(decideTests(fit))  
  ObVsYc_Ob_vs_Yc<-topTable(fit, coef="DiseaseObesity", confint=TRUE, n = Inf)
  ObVsYc_Ob_vs_Yc_Gene<-mutate(as.data.frame(ObVsYc_Ob_vs_Yc), Gene=row.names(ObVsYc_Ob_vs_Yc))
  ObVsYc_Ob_vs_Yc_Gene<-select(ObVsYc_Ob_vs_Yc_Gene, Gene, everything())
  ObVsYc_Ob_vs_Yc_ID<-mutate(ObVsYc_Ob_vs_Yc_Gene, ID=gene_anno[match(ObVsYc_Ob_vs_Yc_Gene$Gene,gene_anno$Id),][,4], Chr=gene_anno[match(ObVsYc_Ob_vs_Yc_Gene$Gene,gene_anno$Id),][,9], Start=gene_anno[match(ObVsYc_Ob_vs_Yc_Gene$Gene,gene_anno$Id),][,11], End=gene_anno[match(ObVsYc_Ob_vs_Yc_Gene$Gene,gene_anno$Id),][,12], Obesity_Ave_FPKQ=fpkq[match(ObVsYc_Ob_vs_Yc_Gene$Gene,fpkq$Gene),][,ncol(fpkq)-1], YoungCon_Ave_FPKQ=fpkq[match(ObVsYc_Ob_vs_Yc_Gene$Gene,fpkq$Gene),][,ncol(fpkq)])
  ObVsYc_Ob_vs_Yc_ID<-select(ObVsYc_Ob_vs_Yc_ID, ID, everything())
  ObVsYc_Ob_vs_Yc_ID<-mutate(ObVsYc_Ob_vs_Yc_ID, ENSEMBL=substring(ObVsYc_Ob_vs_Yc_ID$ID, regexpr("EN", ObVsYc_Ob_vs_Yc_ID$ID)))
  ObVsYc_Ob_vs_Yc_ID<-select(ObVsYc_Ob_vs_Yc_ID, ENSEMBL, everything())
  write.table(ObVsYc_Ob_vs_Yc_ID %>% mutate(ID=str_replace(ID, "'", "")), paste( str_replace(tissue[i]," ", "_"),"_OBVsNonOB.txt", sep=""),row.names = F, sep="\t",quote=F)
}




# Below script describes the workflow of differential RNA expression analysis between T2D vs. non-diabetic monkeys for all tissues.
# Load libraries.
library(limma)
library(dplyr)
library(edgeR)
library(stringr)
library(ggplot2)
# Read gene annotation info.
gene_anno<- read.table("feature.Gene.anno",head = T, sep="\t")
rownames(gene_anno)<-gene_anno$Id
head(gene_anno)
# Read sample annotation info.
sample_anno<- read.table("NHP_Sample_Anno_CRO.txt",head = T, sep="\t")
rownames(sample_anno)<-sample_anno$ID
head(sample_anno)
tissue<-sort(unique(sample_anno$Tissue))
for(i in c(1:4,6:17,19:25,27)){
  tissue_sample<-read.table(paste(tissue[i], "_T2DVsNonT2D_Sample_Anno.txt", sep=""),sep="\t",header=T)
  tissue_count<-read.table(paste(tissue[i], "_T2DVsNonT2D_Raw_Count.txt", sep=""),sep="\t",header=T)
  fpkq<-read.table(paste(tissue[i], "_T2DVsNonT2D_FPKQ.txt", sep=""),sep="\t",header=T)
  row.names(tissue_count)=tissue_count$ID
  table(tissue_sample$ID==colnames(tissue_count[,-1]))
  y <- DGEList(counts = tissue_count[,-1])
  y <- calcNormFactors(y)
  design <- model.matrix(~ Group + Sex, data = tissue_sample)
  v <- voom(y, design, plot=TRUE, normalize="quantile")
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  summary(decideTests(fit))
  DiVsAc_DM_vs_AC<-topTable(fit, coef="GroupDiabetes", confint=TRUE, n = Inf)
  DiVsAc_DM_vs_AC_Gene<-mutate(as.data.frame(DiVsAc_DM_vs_AC), Gene=row.names(DiVsAc_DM_vs_AC))
  DiVsAc_DM_vs_AC_Gene<-select(DiVsAc_DM_vs_AC_Gene, Gene, everything())
  DiVsAc_DM_vs_AC_ID<-mutate(DiVsAc_DM_vs_AC_Gene, ID=gene_anno[match(DiVsAc_DM_vs_AC_Gene$Gene,gene_anno$Id),][,4], Chr=gene_anno[match(DiVsAc_DM_vs_AC_Gene$Gene,gene_anno$Id),][,9], Start=gene_anno[match(DiVsAc_DM_vs_AC_Gene$Gene,gene_anno$Id),][,11], End=gene_anno[match(DiVsAc_DM_vs_AC_Gene$Gene,gene_anno$Id),][,12], Diabetes_Ave_FPKQ=fpkq[match(DiVsAc_DM_vs_AC_Gene$Gene,fpkq$Gene),][,ncol(fpkq)-1], AgedCon_Ave_FPKQ=fpkq[match(DiVsAc_DM_vs_AC_Gene$Gene,fpkq$Gene),][,ncol(fpkq)])
  DiVsAc_DM_vs_AC_ID<-select(DiVsAc_DM_vs_AC_ID, ID, everything())
  DiVsAc_DM_vs_AC_ID<-mutate(DiVsAc_DM_vs_AC_ID, ENSEMBL=substring(DiVsAc_DM_vs_AC_ID$ID, regexpr("EN", DiVsAc_DM_vs_AC_ID$ID)))
  DiVsAc_DM_vs_AC_ID<-select(DiVsAc_DM_vs_AC_ID, ENSEMBL, everything())
  write.table(DiVsAc_DM_vs_AC_ID %>% mutate(ID=str_replace(ID, "'", "")), paste( str_replace(tissue[i]," ", "_"),"_T2DVsNonT2D.txt", sep=""),row.names = F, sep="\t",quote=F)
}




