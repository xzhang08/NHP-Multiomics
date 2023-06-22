# Below script describes the workflow of differential protein expression analysis between obese vs. non-obese monkeys for all tissues.
library(limma)
library(dplyr)
library(edgeR)
library(stringr)
library(ggplot2)
## Read gene annotation info.
gene_anno<- read.table("Protein_Gene_ID.txt",head = T, sep="\t")
rownames(gene_anno)<-gene_anno$Protein
sample_anno<- read.table("Obesity_Proteomics_Sample_Anno.txt",head = T, sep="\t")
rownames(sample_anno)<-sample_anno$ID
tissue=c("Adrenal_Gland","Biceps_Femoris","Brain_Stem","Brown_Adipose_Tissue","Cerebellum","Esophagus","Left_Atrium","Left_Ventricle","Liver","Lumbar_Spinal_Cord","Lung","Occipital_Lobe","Optic_Nerve","Pancreas","Pituitary_Gland","Renal_Cortex","Renal_Medulla","Right_Atrium","Right_Ventricle","Salivary_Gland","Subcutaneous_Fat","Thalamus","Visceral_Fat", "Uterus")
for(i in 1:length(tissue)){
  tissue_count<-read.table(paste("Norm_Obesity_", tissue[i], ".txt", sep=""),sep="\t",header=T)
  row.names(tissue_count)=tissue_count$Accession
  if(tissue[i] == "Uterus"){
    tissue_sample<- read.table("Obesity_Proteomics_Uterus_Sample_Anno.txt",head = T, sep="\t")
    tissue_count<-select(tissue_count, Accession, Human_id, tissue_sample$ID)
    table(tissue_sample$ID==colnames(tissue_count[,c(-1,-2)]))
    design <- model.matrix(~ Disease, data = tissue_sample)
    y <- tissue_count[,c(-1,-2)]
    E <- new("EList", list(E=y))
    fit <- lmFit(E, design)
    fit <- eBayes(fit)
    summary(decideTests(fit))
  } 
  else {
    tissue_sample<-sample_anno
    tissue_count<-select(tissue_count, Accession, Human_id, sample_anno$ID)
    table(tissue_sample$ID==colnames(tissue_count[,c(-1,-2)]))
    design <- model.matrix(~ Disease + Sex, data = tissue_sample)
    y <- tissue_count[,c(-1,-2)]
    E <- new("EList", list(E=y))
    fit <- lmFit(E, design)
    fit <- eBayes(fit)
    summary(decideTests(fit))
  }
ObVsYc_Ob_vs_Yc<-topTable(fit, coef="DiseaseObesity", confint=TRUE, n = Inf)
ObVsYc_Ob_vs_Yc<-mutate(ObVsYc_Ob_vs_Yc, Accession = rownames(ObVsYc_Ob_vs_Yc))
ObVsYc_Ob_vs_Yc<-mutate(ObVsYc_Ob_vs_Yc, Gene = gene_anno[match(ObVsYc_Ob_vs_Yc$Accession,gene_anno$Protein),2], Protein = gene_anno[match(ObVsYc_Ob_vs_Yc$Accession,gene_anno$Protein),1])
write.table(ObVsYc_Ob_vs_Yc, paste(tissue[i], "_OBVsNonOB.txt", sep=""), row.names = T, sep="\t",quote=F)
}



# Below script describes the workflow of differential protein expression analysis between T2D vs. non-diabetic monkeys for all tissues.
library(limma)
library(dplyr)
library(edgeR)
library(stringr)
library(ggplot2)
## Read gene annotation info.
gene_anno<- read.table("Protein_Gene_ID.txt",head = T, sep="\t")
rownames(gene_anno)<-gene_anno$Protein
head(gene_anno)
sample_anno<- read.table("T2D_Proteomics_Sample_Anno.txt",head = T, sep="\t")
rownames(sample_anno)<-sample_anno$ID
head(sample_anno)
tissue=c("Adrenal_Gland","Bicep_Femoris","Brain_Stem","Cerebellum","Esophagus","Left_Atrium","Left_Ventricle","Liver","Lumbar_Spinal_Cord","Lung","Occipital_Lobe","Optic_Nerve","Pancreas","Pituitary_Gland","Renal_Cortex","Renal_Medulla","Right_Atrium","Right_Ventricle","Salivary_Gland","Subcutaneous_Fat","Visceral_Fat","Plasma")
for(i in 1:length(tissue)){
  tissue_count<-read.table(paste(tissue[i], "_T2D_Proteomics.txt", sep=""), sep="\t", header=T)
  row.names(tissue_count)=tissue_count$Accession
  tissue_sample<-sample_anno
  tissue_count<-select(tissue_count, Accession, sample_anno$ID)
  table(tissue_sample$ID==colnames(tissue_count[,-1]))
  design <- model.matrix(~ Disease + Sex, data = tissue_sample)
  y <- tissue_count[,-1]
  E <- new("EList", list(E=y))
  fit <- lmFit(E, design)
  fit <- eBayes(fit)
  summary(decideTests(fit))
  T2DVsAc_T2D_vs_Ac<-topTable(fit, coef="DiseaseT2D", confint=TRUE, n = Inf)
  T2DVsAc_T2D_vs_Ac<-mutate(T2DVsAc_T2D_vs_Ac, Accession = rownames(T2DVsAc_T2D_vs_Ac))
  T2DVsAc_T2D_vs_Ac<-mutate(T2DVsAc_T2D_vs_Ac, Gene = gene_anno[match(T2DVsAc_T2D_vs_Ac$Accession,gene_anno$Protein),2], Protein = gene_anno[match(T2DVsAc_T2D_vs_Ac$Accession,gene_anno$Protein),1])
  write.table(T2DVsAc_T2D_vs_Ac, paste(tissue[i], "_T2DVsNonT2D.txt", sep=""), row.names = T, sep="\t",quote=F)
}




