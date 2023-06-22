# Below script describes the workflow of WGCNA of RNA-seq data for the subcutaneous abdominal white adipose tissue.
library(dplyr)
library(WGCNA)
library(stringr)
library(ggplot2)
library(limma)
library(edgeR)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(reshape2)
### Read Sample Anno file.
sample_anno<- read.table("Sample_Anno_T2D_NonT2D_OB_NonOB_Same_Tis.txt",head = T, sep="\t")
rownames(sample_anno)<-sample_anno$ID
head(sample_anno)
tissue<-sort(unique(sample_anno$Tissue))
### Select the tissue to be analyzed.
sample_anno<-filter(sample_anno, Tissue == "Subcutaneous Fat")
### Read gene count
gene_count<-read.table("QCed_primary_count_gene.Counts.Table.txt",sep="\t",header=T)
row.names(gene_count)=gene_count$ID
### Select the gene count for the tissue to be analyzed.
gene_count<-select(gene_count, ID, sample_anno$ID)
### Check sample order consistency.
table(colnames(gene_count[,-1])==sample_anno$ID)
### Quantile normalize the data.
y <- DGEList(counts = gene_count[,-1])
dim(y)
## Keep genes with total counts more than 50.
A <- rowSums(y$counts) 
isexpr <- A > 100
y <- y[isexpr, keep.lib.size = FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~ Sex + Group, data = sample_anno)
v <- voom(y, design, plot=FALSE, normalize="quantile")
write.table(v[["E"]],"Quantile_Norm_T2D_NonT2D_OB_NonOB_count.txt", sep = "\t", quote = F)
### WGCNA starts
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
Data = read.table("Quantile_Norm_T2D_NonT2D_OB_NonOB_count.txt",header=T)
# We now transpose the expression data for further analysis.
datExpr = as.data.frame(t(Data));
## Loading clinical trait data
traitData = read.table("Sample_Anno_T2D_NonT2D_OB_NonOB_Same_Tis.txt", sep = "\t", header = T)
# Form a data frame analogous to expression data that will hold the clinical traits.
Samples = rownames(datExpr)
traitRows = match(Samples, traitData$ID)
datTraits = traitData[traitRows, -1]
rownames(datTraits) = traitData[traitRows, 1]
collectGarbage()
## Choosing the soft-thresholding power: analysis of network topology
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=40, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
pdf("Soft_thrshold_power.pdf", width = 9, height = 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
############################################
##### Choose power of 5 based on plot.
############################################
tissue<-sort(unique(sample_anno$Tissue))
### Select the tissue to be analyzed.
sample_anno<-filter(sample_anno, Tissue == "Subcutaneous Fat")
net = blockwiseModules(datExpr, maxBlockSize = 20000, power = 5, networkType = "signed",TOMType = "signed", nThreads=20, minModuleSize = 20,reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE,saveTOMs = TRUE,saveTOMFileBase = "T2D_AC_OB_YC_Same_Tis_TOM",verbose = 3)
### save the module assignment and module eigengene information necessary for subsequent analysis
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,file = "T2D_NonT2D_OB_NonOB_Same_Tis-02-networkConstruction-auto.RData")
datME=moduleEigengenes(datExpr,moduleColors)$eigengenes
signif(cor(datME, use="p"), 2)
print(signif(cor(datME, use="p"), 2))
dissimME=(1-t(cor(datME, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )
#### Get the number of genes for each module
novak<-unique(moduleColors)
module_gene_count<-matrix(rep(0,length(novak)*2),length(novak),2)
for(i in 1:length(novak)){
#print(novak[i])
#print(ncol(datExpr[,moduleColors==novak[i] ]))
module_gene_count[i,]=c(novak[i],ncol(datExpr[,moduleColors==novak[i] ]))
}
### Correlate the module eigengenes with the trait
phenotype<-read.table("Phenotypes.txt", header = T)
vars <- colsplit(datTraits$Subject, "-", c("Donor", "Tis_Index"))
vars
datTraits2<-mutate(datTraits, Donor = vars$Donor)
datTraits3<-merge(datTraits2,phenotype, by = "Donor")
rownames(datTraits3)<-rownames(datTraits)
datTraits3$SR_ID==datTraits$SR_ID
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
datTraits4=datTraits3[,12:25]
moduleTraitCor = cor(MEs, datTraits4, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
pdf("Trait_Eigen_Cor.pdf",width=10,height=7)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot 
labeledHeatmap(Matrix = moduleTraitCor,
             xLabels = names(datTraits4),
             yLabels = names(MEs),
             ySymbols = names(MEs),
             colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5,zlim = c(-1,1),
main = paste("Module-trait relationships"))
dev.off()
### Correlate the module eigengenes with the trait
signif(cor(datTraits3$Phenotype,datME, use="p"),2)
print(signif(cor(datTraits3$Phenotype,datME, use="p"),2))
##### The following code can be used to get all p-values:
p.values = corPvalueStudent(cor(datTraits3$Phenotype,datME, use="p"), nSamples = length(datTraits3$Phenotype))
print(p.values)
### Measure of module significance as average gene significance
GS1=as.numeric(cor(datTraits3$Phenotype,datExpr, use="p"))
GeneSignificance=abs(GS1)
write.table(GeneSignificance,"GeneSignificance.txt")
# Next module significance is defined as average gene significance.
ModuleSignificance=tapply(GeneSignificance, moduleColors, mean, na.rm=T)
print(ModuleSignificance)
################## Plot the modules of interest only.
datME=moduleEigengenes(datExpr,moduleColors)$eigengenes
novak<-unique(moduleColors)
df<-as.data.frame(sample_anno[,10:11])
df<-arrange(df, match(Group, c("NonOB", "NonT2D", "OB", "T2D")))
tmp<-data.matrix(t(scale(datExpr)))
tmp<-dplyr::select(as.data.frame(tmp),rownames(df))
plot_add<-mutate(as.data.frame(tmp), module = moduleColors)
rownames(plot_add)<-rownames(tmp)
plot_add<-filter(plot_add, module == "green" | module == "brown")
ann_colors = list(module = c(brown = "brown", green = "green"))
plot_order<-arrange(plot_add, module)
plot_order_gene<-mutate(plot_order, gene=rownames(plot_order))
anno_row<-as.data.frame(plot_order_gene[,(ncol(plot_order_gene)-1):ncol(plot_order_gene)])
tmp_rownames<-anno_row$gene
anno_row<-dplyr::select(anno_row, module)
df<-dplyr::select(df, Group)
rownames(anno_row)<-tmp_rownames
plot_order<-dplyr::select(plot_order, rownames(df))
pdf("plasma_ModuleOfInterest_heatmap.pdf")
pheatmap(plot_order, cluster_rows=F, cluster_cols=F, scale = "none", show_rownames=FALSE, show_colnames=FALSE, annotation_col=df, annotation_row = anno_row, color = plasma(200),na_col = "black", annotation_colors = ann_colors[1])
dev.off()





