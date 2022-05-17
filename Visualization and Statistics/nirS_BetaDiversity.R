############################################################
##################### Beta-diversity  ######################
############################################################

BiocManager::install("microbiome")
library(phyloseq)
library(ggplot2)
library(plyr)
library(gridExtra)
library(BiocManager)
library(microbiome)
library(ggplot2)
library(gridExtra)
library(grid)
library(ggpubr)

col_Symb=c("#CECACC", "#E0A800", "#006C67") #APO, SymbA, SymbB
col_host=c("#AA4388", "#167777", "#777611")

setwd("~/WP3_Aiptasia_nirS/nirS_analysis/")
load("nirS_ordination.Rdata")
save.image("nirS_ordination.Rdata")

asv=read.table("~/WP3_Aiptasia_nirS/nirS_analysis/R_output/ASV_Tax_noConta.txt", header = TRUE, sep ='\t')
met=read.table("~/WP3_Aiptasia_nirS/nirS_analysis/R_input/Aip_nirS_metadata.txt", header = TRUE, row.names = 1, sep ='\t')
met[is.na(met)] <- "None"

met$host=factor(met$host, levels = c("None", "CC7", "H2"))
met$symbiont=factor(met$symbiont, levels = c("None","SSA01", "SSB01"))
asv$CA3=NULL # remove CA3,  samples without DNA input, 7704 reads
asv$CA5=NULL # remove CA5,  samples without DNA input, 11666 reads

otu.t=otu_table(asv[, 1:41], taxa_are_rows=TRUE)
sam.t= sample_data(data.frame(met))
tax.t= tax_table(as.matrix(asv[, 43:ncol(asv)]))
phy.all= phyloseq(otu.t, tax.t,  sam.t)

####################
#### Ordination ####
####################

#### PCA by Subsetting host ############
phy.t=microbiome::transform(phy.all, transform = "clr", target = "OTU", shift = 0, scale = 1) 
PCA = ordinate(phy.t, method = "RDA", distance = "euclidean")

phy_CC7=subset_samples(phy.t, host == "CC7")
phy_H2=subset_samples(phy.t, host == "H2")
phy_AF=subset_samples(phy.t, host == "None")

PCA_CC7 = ordinate(phy_CC7, method = "RDA", distance = "euclidean")
PCA_H2 = ordinate(phy_H2, method = "RDA", distance = "euclidean")
PCA_AF = ordinate(phy_AF, method = "RDA", distance = "euclidean")

plot_CC7=plot_ordination(phy_CC7, PCA_CC7, color = "symbiont", shape = "host") + 
  geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("") + 
  theme_classic() + scale_colour_manual(values=c(col_Symb)) +
  #facet_wrap(~host, ncol=3, scales="free") +
  ggtitle("PCA on Euclidean distance")

plot_H2=plot_ordination(phy_H2, PCA_H2, color = "symbiont", shape = "host") + 
  geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("") + 
  theme_classic() + scale_colour_manual(values=c(col_Symb)) +
  #facet_wrap(~host, ncol=3, scales="free") +
  ggtitle("PCA on Euclidean distance")

plot_AF=plot_ordination(phy_AF, PCA_AF, color = "symbiont", shape = "host") + 
  geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("") + 
  theme_classic() + scale_colour_manual(values=c(col_Symb)) +
  #facet_wrap(~host, ncol=3, scales="free") +
  ggtitle("PCA on Euclidean distance")


pdf("/Users/nanxiang/Desktop/PhD/WP3_Aiptasia_nirS/nirS_analysis/R_output/nirS_PCA_new.pdf", width=6,height=2.5, pointsize = 12)
ggarrange(plot_AF, plot_CC7, plot_H2 + rremove("x.text"),
          legend = "none",
          labels = c("Algae and food", "CC7", "H2"),
          ncol = 3, nrow = 1)
dev.off()






