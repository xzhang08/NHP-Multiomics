# Below script generates Figure 7A.
library(networkD3)
library(dplyr)
library(RColorBrewer)
links<-read.table("SubFat_OB_T2D_DEPs_Sankey.txt", sep="\t", header=T)
nodes <- data.frame(
  name=c(as.character(links$Source), 
  as.character(links$Target)) %>% unique()
)
links$IDsource <- match(links$Source, nodes$name)-1 
links$IDtarget <- match(links$Target, nodes$name)-1
links$group <- as.factor(c(rep("type_a",3),rep("type_b",3),rep("type_c",2),rep("type_d",3)))
my_color <- 'd3.scaleOrdinal() .domain(["Control","Upregulated_In_Obesity","Unchanged_In_Obesity","Downregulated_In_Obesity","Upregulated_In_T2D","Unchanged_In_T2D","Downregulated_In_T2D","type_a","type_b","type_c","type_d"]) .range(["#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#66C2A5","#FC8D62","#8DA0CB","#E78AC3"])'
# Make the Network
sankeyNetwork(Links = links, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "Value", NodeID = "name", sinksRight=FALSE, colourScale = my_color, LinkGroup="group", NodeGroup="name", nodeWidth=15, fontSize=13, fontFamily="Arial", nodePadding=20)



# Below script generates Figure 7B.
library(networkD3)
library(dplyr)
library(RColorBrewer)
links<-read.table("VisFat_OB_T2D_DEPs_Sankey.txt", sep="\t", header=T)
nodes <- data.frame(
  name=c(as.character(links$Source), 
  as.character(links$Target)) %>% unique()
)
links$IDsource <- match(links$Source, nodes$name)-1 
links$IDtarget <- match(links$Target, nodes$name)-1
links$group <- as.factor(c(rep("type_a",3),rep("type_b",2),rep("type_c",2),rep("type_d",3)))
my_color <- 'd3.scaleOrdinal() .domain(["Control","Upregulated_In_Obesity","Unchanged_In_Obesity","Downregulated_In_Obesity","Upregulated_In_T2D","Unchanged_In_T2D","Downregulated_In_T2D","type_a","type_b","type_c","type_d"]) .range(["#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#66C2A5","#FC8D62","#8DA0CB","#E78AC3"])'
# Make the Network
sankeyNetwork(Links = links, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "Value", NodeID = "name", sinksRight=FALSE, colourScale = my_color, LinkGroup="group", NodeGroup="name", nodeWidth=15, fontSize=13, fontFamily="Arial", nodePadding=20)



# Below script generates Figure 7C.
library(dplyr)
library(stringr)
library(ggplot2)
library(pheatmap)
library(viridis)
library(RColorBrewer)
samples<-c("YC1","YC2","YC3","OW2","AC1","AC2","AC3","AC4","OB3","PD2","PD3","OB1","OB2","OW1","PD1","PD4","DM1","DM2","DM3","DM4","DM5","DM6")
samples_2<-paste0(samples,".17")
sample_anno<-data.frame(ID=samples, Group=c(rep("NonObesity",4), rep("NonT2D",4), rep("Obesity",8),rep("T2D",6)), Gender=c(rep("F",3), "M", rep("M",2), rep("F",2), rep("F",3), rep("M",5), rep("M",3), rep("F",3)))
rownames(sample_anno)=samples_2
normalized_pro<- read.table("Normalized_data_All_samples.txt",head = T, quote = "", sep="\t")
sub_norm<-select(normalized_pro, Accession, all_of(samples_2))
genes_to_plot<-read.table("Subcutaneous_FatOB&T2D_SameDir_overlapped_DEPs.txt", sep = "\t", quote = "", header = T)
sub_norm_genes<-filter(sub_norm, Accession %in% genes_to_plot$Accession)
keep <- apply(sub_norm_genes, 1, function(x) length(unique(x[!is.na(x)])) != 1)
sub_norm_filter<-sub_norm_genes[keep, ]
rownames(sub_norm_filter)<-genes_to_plot[match(sub_norm_filter$Accession,genes_to_plot$Accession),10]
Var1 <- c(brewer.pal(4, "Paired"))
names(Var1) <- c("NonObesity","Obesity","NonT2D","T2D")
anno_colors <- list(Group = Var1)
pdf("Fig7C.pdf")
pheatmap(sub_norm_filter[,-1], cluster_rows=TRUE, show_rownames=TRUE, show_colnames=FALSE, cluster_cols=FALSE, annotation_col=select(sample_anno,2), scale = "row", treeheight_row = 0, color = inferno(256), annotation_colors = anno_colors, fontsize_row = 8)
dev.off()



# Below script generates Figure 7D.
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(UpSetR)
ob_ab<-read.table("/Users/xzhang08/Desktop/Amgen_China_Cyno/Proteomics/Obesity/limma/Subcutaneous_Fat_ObVsYc_SigOb_vs_Yc.txt",sep="\t", header = T)
ob_vis<-read.table("/Users/xzhang08/Desktop/Amgen_China_Cyno/Proteomics/Obesity/limma/Visceral_Fat_ObVsYc_SigOb_vs_Yc.txt",sep="\t",header = T)
t2d_ab<-read.table("/Users/xzhang08/Desktop/Amgen_China_Cyno/Proteomics/T2D/limma/Subcutaneous_Fat_T2DVsAc_SigT2D_vs_Ac.txt", sep="\t",header = T)
t2d_vis<-read.table("/Users/xzhang08/Desktop/Amgen_China_Cyno/Proteomics/T2D/limma/Visceral_Fat_T2DVsAc_SigT2D_vs_Ac.txt", sep="\t",header = T)
input <- c(
  "Ob_Ab_WAT&Ob_Vis_WAT" = 15,
  "Ob_Ab_WAT&T2D_Ab_WAT" = 71,
  "Ob_Ab_WAT&T2D_Vis_WAT" = 132,
  "Ob_Vis_WAT&T2D_Ab_WAT" = 24,
  "Ob_Vis_WAT&T2D_Vis_WAT" = 21,
  "T2D_Ab_WAT&T2D_Vis_WAT" = 135,
  "Ob_Vis_WAT&T2D_Ab_WAT&T2D_Vis_WAT" = 3,
  "Ob_Ab_WAT&T2D_Ab_WAT&T2D_Vis_WAT" = 41,
  "Ob_Ab_WAT&Ob_Vis_WAT&T2D_Vis_WAT" = 5,
  "Ob_Ab_WAT&Ob_Vis_WAT&T2D_Ab_WAT" = 5,
  "Ob_Ab_WAT&Ob_Vis_WAT&T2D_Ab_WAT&T2D_Vis_WAT" = 3
)
upset(fromExpression(input), 
      nintersects = 11, 
      nsets = 4, 
      order.by = "freq",
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 1.1, 
      point.size = 2.8, 
      line.size = 1
)


