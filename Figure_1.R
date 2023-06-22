# Below script generates the scatterpie plot in Figure 1C.
library(ggplot2)
library(scatterpie)
library(dplyr)
library(stringr)
library(RColorBrewer)
DEG<-read.table("DEGs_Count_OB&T2D.txt", header = T, sep = "\t")
DEG<-mutate(DEG, Disease_code=c(rep(1,26),rep(3,26)))
DEG<-arrange(DEG, Disease_code, Organ, Tissue_Rollup, Tissue)
DEG<-mutate(DEG, Tissue_index=rep(1:26,2))
DEG_NA<-filter(DEG, is.na(DEGs))
dim(DEG)
dim(DEG_NA)
DEG_plot<-select(DEG, Tissue_index, Disease_code, DEGs,Up_DEGs,Down_DEGs)
leg<-0.04*log2(as.vector(na.omit(DEG_plot$DEGs+1.1)))
leg<-as.numeric(leg)
DEG_NA<-filter(DEG_plot, is.na(DEGs))
ggplot() + geom_scatterpie(aes(x=Tissue_index, y=Disease_code, group = Disease_code, r = 0.04 * log2(DEGs)), data = DEG_plot, cols = c("Up_DEGs","Down_DEGs"))  + geom_scatterpie_legend(leg, x=5, y=5, n=6, labeller=function(x) as.integer(2^(x/0.04))) + theme_bw() + scale_y_continuous(breaks = c(1,3), labels=c("Obesity","T2D")) + xlab("Tissues") + ylab("Disease") + scale_x_continuous(limits=c(0.5,26.5), breaks = 1:26, labels = str_replace_all(as.vector(DEG$Tissue[1:26]), "_", " "), expand = expand_scale(mult = c(0.01, 0.01))) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust =1, color = "black"), axis.text.y = element_text(color = "black"), legend.position = c(0.7,0.88), legend.title = element_blank(), legend.direction="horizontal", legend.key.size = unit(0.5, "cm"), )  + coord_fixed() + theme(panel.grid.minor = element_blank(), legend.background = element_blank()) + scale_fill_manual(values = c( "#CCFFBB", "#005A04"), breaks = c("Up_DEGs","Down_DEGs"), labels = c("Upregulated DEGs","Downregulated DEGs")) + annotate(geom="text", x=DEG_NA$Tissue_index, y=DEG_NA$Disease_code, label="NA", color="black", size = 3)
ggsave("DEGs_Count_OB&T2D.pdf")



# Below script generates the scatterpie plot in Figure 1D.
library(ggplot2)
library(scatterpie)
library(dplyr)
library(stringr)
library(RColorBrewer)

DEP<-read.table("OB_T2D_DEPs_Count.txt", header = T, sep = "\t")
DEP<-mutate(DEP, Disease_code=c(rep(1,25),rep(3,25)))
DEP<-arrange(DEP, Disease_code, Organ, Tissue_Rollup, Tissue)
DEP<-mutate(DEP, Tissue_index=rep(1:25,2))

DEP_NA<-filter(DEP, is.na(DEPs))
dim(DEP)
dim(DEP_NA)

DEP_plot<-select(DEP, Tissue_index, Disease_code, DEPs,Up_DEPs,Down_DEPs)
leg<-0.04*log2(as.vector(na.omit(DEP_plot$DEPs+1.1)))
leg<-as.numeric(leg)
DEP_NA<-filter(DEP_plot, is.na(DEPs))

ggplot() + geom_scatterpie(aes(x=Tissue_index, y=Disease_code, group = Disease_code, r = 0.04 * log2(DEPs)), data = DEP_plot, cols = c("Up_DEPs","Down_DEPs"))  + geom_scatterpie_legend(leg, x=5, y=4.5, n=3, labeller=function(x) as.integer(2^(x/0.04))) + theme_bw() + scale_y_continuous(breaks = c(1,3), labels=c("Obesity","T2D")) + xlab("Tissues") + ylab("Disease") + scale_x_continuous(limits=c(0.5,25.5), breaks = 1:25, labels = str_replace_all(as.vector(DEP$Tissue[1:25]), "_", " "), expand = expand_scale(mult = c(0.01, 0.01))) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust =1, color = "black"), axis.text.y = element_text(color = "black"), legend.position = c(0.7,0.88), legend.title = element_blank(), legend.direction="horizontal", legend.key.size = unit(0.5, "cm"), )  + coord_fixed() + theme(panel.grid.minor = element_blank(), legend.background = element_blank()) + scale_fill_manual(values = brewer.pal(8, "Paired")[3:4], breaks = c("Up_DEPs","Down_DEPs"), labels = c("Upregulated DEPs","Downregulated DEPs")) + annotate(geom="text", x=DEP_NA$Tissue_index, y=DEP_NA$Disease_code, label="NA", color="black", size = 3)
ggsave("OB_T2D_DEPs_Count.pdf")


























