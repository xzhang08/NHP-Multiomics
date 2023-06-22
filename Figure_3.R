# Below script generates Figure 3A.
library(ggplot2)
library(RColorBrewer)
library(dplyr)
ab_wat<-read.table("Fig3A_Data.txt", header = T, sep = "\t")
ab_wat$Canonical.Pathways=factor(ab_wat$Canonical.Pathways, levels=rev(ab_wat$Canonical.Pathways))
ggplot(ab_wat, aes(x = Abdominal.WAT, y =  Canonical.Pathways)) + 
  geom_bar(stat = "identity", width=0.5, fill="#99CC99") + 
  xlab("-log10(adjusted p-values)") +
  ylab("Canonical Pathways") +
  theme_bw() +
  theme(axis.text = element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 10)) 
ggsave("Fig3A.pdf")



# Below script generates Figure 3C.
library(ggplot2)
library(RColorBrewer)
library(dplyr)
ab_wat<-read.table("Fig3C_Data.txt", header = T, sep = "\t")
ab_wat$Canonical.Pathways=factor(ab_wat$Canonical.Pathways, levels=rev(ab_wat$Canonical.Pathways))
ggplot(ab_wat, aes(x = Abdominal.WAT, y =  Canonical.Pathways)) + 
  geom_bar(stat = "identity", width=0.5, fill="#B3B3D7") + 
  xlab("-log10(adjusted p-values)") +
  ylab("Canonical Pathways") +
  theme_bw() +
  theme(axis.text = element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 18)) 
ggsave("Fig3C.pdf")




