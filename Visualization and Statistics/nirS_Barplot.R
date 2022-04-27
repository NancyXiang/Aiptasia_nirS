#############################################################
#####   Taxonomic profiles of the 10 most abundant ASVs #####
#############################################################

library(reshape2)
library(ggplot2)
library(scales)
library(dplyr)
library(ggpubr)


setwd("~/WP3_Aiptasia_nirS/nirS_analysis/")
load("nirS_Barplot.Rdata")
save.image("nirS_Barplot.Rdata")

asv = read.table("~/WP3_Aiptasia_nirS/nirS_analysis/R_output/ASV_Tax_noConta_formated.txt", header = T, row.names = 1, sep = "\t")
map = read.table("~/WP3_Aiptasia_nirS/nirS_analysis/R_input/Aip_nirS_metadata.txt", header = T, sep = "\t", row.names = 1)
map[is.na(map)] <- "None"
asv[is.na(asv)] <- "unclassified"

# keep in parallel with alphaDiversity and Venn, remove CA3, CA5 (no DNA input during library prep)
asv$CA3=NULL
asv$CA5=NULL
names(asv)

asv$ASV <- rownames(asv)
asv$ASV=gsub("sq", "ASV", asv$ASV)

asv.tax.ag=aggregate(asv[, 1:41], by = list(asv[, 51]), FUN =  sum) #define sample range and group factor

topASV=asv.tax.ag[order(rowSums(asv.tax.ag[,2:ncol(asv.tax.ag)]),decreasing = TRUE),][1:10,1]
ASV.top=subset(asv.tax.ag, asv.tax.ag$Group.1 %in% topASV)
ASV.bot=subset(asv.tax.ag, !asv.tax.ag$Group.1 %in% topASV)
ASV.bot$Group.1=gsub(".*","Others", ASV.bot$Group.1)
others=aggregate(ASV.bot[, 2:ncol(ASV.bot)], by = list(ASV.bot[, 1]), FUN =  sum)

all.2 =rbind(ASV.top, others)
all.l=reshape2::melt(all.2, id.vars=c("Group.1"), variable.name = "ASV", value.name = "Abundance")
colnames(all.l)=c("ASV","Sample","Abundance")

## Add sample information
all.l$host=map$host[match(all.l$Sample, rownames(map))]
all.l$symbiont=map$symbiont[match(all.l$Sample, rownames(map))]
all.l$other=map$other[match(all.l$Sample, rownames(map))]

P11=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#CCCCCC")

all.l$Taxa=paste(asv$Phylum, asv$Family, asv$Genus,sep = "_")[match(all.l$ASV, asv$ASV)]

CC7_all=ggplot() +geom_bar(aes(y = Abundance, x = Sample, fill = ASV), data = subset(all.l, host == "CC7" ), stat="identity", position = "fill") +  theme(axis.text.x=element_text(angle=90,hjust=1)) + labs( y= "Percentage of nirS sequences", x="Host-symbiont combinations", main = "Spring") + scale_fill_manual(values=P11) + facet_wrap(~symbiont, ncol=3, scales="free") + theme( legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'bottom') +  theme_classic() + guides(fill=guide_legend(ncol=1))
H2_all=ggplot() +geom_bar(aes(y = Abundance, x = Sample, fill = ASV), data = subset(all.l, host == "H2" ), stat="identity", position = "fill") +  theme(axis.text.x=element_text(angle=90,hjust=1)) + labs( y= "Percentage of nirS sequences", x="Host-symbiont combinations", main = "Spring") + scale_fill_manual(values=P11) + facet_wrap(~symbiont, ncol=3, scales="free") + theme( legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'bottom') +  theme_classic() + guides(fill=guide_legend(ncol=1))

# food,algae
FooAlgae=subset(all.l, host == "None")
write.table(FooAlgae, "R_output/FooAlgae_barplot", quote =F, sep = "\t")

asv.FooAlgae = read.table("~/WP3_Aiptasia_nirS/nirS_analysis/R_output/FooAlgae_barplot.txt", header = T, sep = "\t")
asv.FooAlgae[is.na(asv.FooAlgae)] <- "others"

Food_Algae=ggplot() +geom_bar(aes(y = AbundanceMerged, x = Sample, fill = ASV), data = asv.FooAlgae, stat="identity", position = "fill") +  theme(axis.text.x=element_text(angle=90,hjust=1)) + labs( y= "Percentage of nirS sequences", x="Algal strains") + scale_fill_manual(values=P11) + facet_wrap(~symbiont, ncol=3, scales="free") + theme( legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'bottom') +  theme_classic() + guides(fill=guide_legend(ncol=1))

pdf("/Users/nanxiang/Desktop/PhD/WP3_Aiptasia_nirS/nirS_analysis/R_output/nirS_barplot_ASV10_merged.pdf", width=6.5,height=6, pointsize = 12)
ggarrange(Food_Algae, CC7_all, H2_all + rremove("x.text"),
          #legend = "none",
          labels = c("A", "B", "C"),
          ncol = 1, nrow = 3)
dev.off()


