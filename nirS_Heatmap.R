#######################################
###   Heatmap after ANCOM analysis  ###
#######################################

library(ANCOMBC)
library(phyloseq)
library("ape")
library(tibble)
library(ComplexHeatmap)
library(RColorBrewer)

setwd("~/WP3_Aiptasia_nirS/nirS_analysis/")

## input
ASV.enriched=read.table("~/WP3_Aiptasia_nirS/nirS_analysis/R_input/Ancom_enriched_ASV.txt", header = TRUE, sep ='\t')
met=read.table("~/WP3_Aiptasia_nirS/nirS_analysis/R_input/Aip_nirS_metadata.txt", header = TRUE, row.names = 1, sep ='\t')
asv=read.table("~/WP3_Aiptasia_nirS/nirS_analysis/R_output/ASV_Tax_noConta_formated.txt", header = TRUE, row.names = 1, sep ='\t')

sam.t= sample_data(data.frame(met))
sam.t[is.na(sam.t)] <- "None"
sam.t$comb=paste(sam.t$host, sam.t$symbiont, sam.t$other, sep = "_")
sam.t$group1=ifelse(sam.t$comb %in% c("CC7_None_None","H2_None_None"), "APO", ifelse(sam.t$comb == "None_None_food", "Food",  ifelse(sam.t$symbiont == "SSA01", "SSA01", "SSB01" )))
row_names_to_remove<-c("NEC1","NEC2","NEC3","NPC1","NPC2", "CA3", "CA5")
sam.t=sam.t[!(row.names(sam.t) %in% row_names_to_remove),]

# phyloseq
otu.t= otu_table(as.matrix(asv[, 1:43]), taxa_are_rows=TRUE)
tax.t= tax_table(as.matrix(asv[, 45:ncol(asv)]))
physeq = phyloseq(otu.t,tax.t)

random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
plot(random_tree)

physeq1 = merge_phyloseq(physeq, sam.t, random_tree)

ps <- physeq1
ps.taxa <- tax_glom(ps, taxrank = 'ASV', NArm = FALSE)
ps.taxa.rel <- transform_sample_counts(ps, function(x) x/sum(x)*100)

####################### Both hosts ########################
# subset sample
ps.taxa.sub.both <- subset_samples(ps.taxa, host %in% c("CC7", "H2"))

# enriched ASVs
taxa_sig.both=ASV.enriched$both

ps.taxa.rel.sig.both <- prune_taxa(taxa_sig.both, ps.taxa.rel)

ps.taxa.rel.sig.both <- prune_samples(colnames(otu_table(ps.taxa.sub.both)), ps.taxa.rel.sig.both)

## heatmap
matrix.both <- as.matrix(data.frame(otu_table(ps.taxa.rel.sig.both)))
rownames(matrix.both) <- as.character(tax_table(ps.taxa.rel.sig.both)[, "ASV"])
metadata_sub.both <- data.frame(sample_data(ps.taxa.rel.sig.both))

# Define the annotation color for columns and rows

annotation_col = data.frame(
  Host = as.factor(metadata_sub.both$host),
  Symbiont = as.factor(metadata_sub.both$symbiont),
  check.names = FALSE
)
rownames(annotation_col) = rownames(metadata_sub.both)


annotation_row = data.frame(
  Genus = as.factor(tax_table(ps.taxa.rel.sig.both)[, "Genus"])
)
rownames(annotation_row) = rownames(matrix.both)

# ann_color should be named vectors
col=c("#771155","#CC99BB","#4477AA","#117777","#77CCCC","#44AA77","#777711","#DDDD77","#AA7744","#771122")
names(col) = levels(annotation_row$Genus)

ann_colors = list(
  Host = c(CC7 = "#771155", H2 = "#114477"),
  Symbiont = c(None = "#CECACC", SSA01 = "#DFA800", SSB01 = "#006C67"),
  Genus = col
)


pdf("~/WP3_Aiptasia_nirS/nirS_analysis/R_output/heatmap_genus.pdf", width=10,height=5, pointsize = 7)
ComplexHeatmap::pheatmap(matrix.both,
                         color = colorRampPalette(c("#011936", "#FFFFFF", "#C9492C"))(50),
                         scale= "row",
                         border_color = "grey90",
                         column_split = annotation_col$Host,
                         annotation_col = annotation_col,
                         annotation_row = annotation_row,
                         annotation_colors = ann_colors)
dev.off()
