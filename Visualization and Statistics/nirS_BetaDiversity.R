############################################################
#####################   NMDS  #############################
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


############ PCA ####################
phy.t=microbiome::transform(phy.all, transform = "clr", target = "OTU", shift = 0, scale = 1)
PCA = ordinate(phy.t, method = "RDA", distance = "euclidean")

pdf("./R_output/new_PCA_ordination.pdf", width=6.5,height=3, pointsize = 12)
plot_ordination(phy.t, PCA, color = "symbiont", shape = "host") +
  geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("") +
  theme_classic() + scale_colour_manual(values=c(col_Symb)) +
  facet_wrap(~host, ncol=3, scales="free") +
  ggtitle("PCA on Euclidean distance")
dev.off()


############ NMDS ####################
phy.t=microbiome::transform(phy.all, transform = "relative.abundance", target = "OTU", shift = 0, scale = 1) #NMDS,
NMDS = ordinate(phy.t, method = "NMDS", distance =  "bray")
pdf("./R_output/new_NMDS_ordination.pdf", width=6.5,height=3, pointsize = 12)
plot_ordination(phy.t,NMDS, color = "symbiont", shape = "host") +
  geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("") +
  theme_classic() + scale_colour_manual(values=c(col_Symb)) +
  facet_wrap(~host, ncol=3, scales="free") +
  ggtitle("NMDS on Bray Curtis distance")
dev.off()


######## PcoA ###############
phy.t=microbiome::transform(phy.all, transform = "relative.abundance", target = "OTU", shift = 0, scale = 1) #PcoA,
PcoA = ordinate(phy.t, method = "PCoA", distance ="bray", weighted=TRUE)
pdf("./R_output/new_PcoA_ordination.pdf", width=6.5,height=3, pointsize = 12)
plot_ordination(phy.t,PcoA, color = "symbiont", shape = "host") +
  geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("") +
  theme_classic() + scale_colour_manual(values=c(col_Symb)) +
  facet_wrap(~host, ncol=3, scales="free") +
  ggtitle("PCoA on Bray Curtis distance")
dev.off()





  geom_point(size = 5, alpha = 1) + theme_bw()  + ggtitle("") +
  theme_classic() + scale_colour_manual(values=c(col_Symb)) +
  #facet_wrap(~host, ncol=3, scales="free") +
  ggtitle("PCoA on Bray Curtis distance")
dev.off()
