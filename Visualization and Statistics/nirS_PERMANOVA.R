#####################################################
######## Stats on community composition ############
####################################################

install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
install.packages('compositions')

library(vegan)
library(devtools)
library(pairwiseAdonis)
library(compositions)
library(ggplot2)
library(ggpubr)

setwd("/Users/nanxiang/Desktop/PhD/WP3_Aiptasia_nirS/nirS_analysis/")

#####################
#### PERMANOVA   ####
#####################

asv=read.table("/Users/nanxiang/Desktop/PhD/WP3_Aiptasia_nirS/nirS_analysis/R_output/ASV_Tax_noConta.txt", header = TRUE, sep ='\t')
met=read.table("/Users/nanxiang/Desktop/PhD/WP3_Aiptasia_nirS/nirS_analysis/R_input/Aip_nirS_metadata.txt", header = TRUE, row.names = 1, sep ='\t')
asv$CA3=NULL # remove CA3, CA5
asv$CA5=NULL

met$factor=ifelse(met$host %in% c("CC7","H2"), "animal", ifelse(met$symbiont %in% c("SSA01","SSB01"), "algal culture",  "food" ))#define new column
met$symbiont[is.na(met$symbiont)] <- "APO"
met$host[is.na(met$host)] <- "None"

asv.n=as.data.frame(t(apply(asv[, 1:41],2,clr))) # clr-transformed counts
asv.n$host=met$host[match(rownames(asv.n), rownames(met))]
asv.n$symbiont=met$symbiont[match(rownames(asv.n), rownames(met))]
asv.n$other=met$other[match(rownames(asv.n), rownames(met))]
asv.n$factor=met$factor[match(rownames(asv.n), rownames(met))]

####################################################################
### comparison1:  overall model test animal，algal cultures,food  ##
####################################################################

## betadispersion
bdisper_overall=betadisper(dist(asv.n[,1:427]), group = asv.n$factor, sqrt.dist = T) # dist: euclidean distance; veg: bray-curtis distance
boxplot(bdisper_overall)
## PERMANOVA
asv_adonis=adonis(asv.n[,1:427]~ asv.n$factor, method = "euclidean")
pairwise.adonis(asv.n[,1:427], asv.n$factor,  sim.method = "euclidean", p.adjust.m = "fdr", perm = 999)

###################################################################
###     comparison2:  subset Six combs - APO vs A vs B         ###
##################################################################

## subset six combinations
comb_six=subset(asv.n, factor == "animal" )
## betadispersion
bdisper_comb_six=betadisper(dist(comb_six[,1:427]), group = comb_six$symbiont, sqrt.dist = T) # dist: euclidean distance; veg: bray-curtis distance
boxplot(bdisper_comb_six)
asv.n$symbiont[is.na(asv.n$symbiont)] <- "APO"
asv.n$host[is.na(asv.n$host)] <- "None"
## PERMANOVA
asv_adonis=adonis(comb_six[,1:427]~ comb_six$symbiont, method = "euclidean")
pairwise.adonis(comb_six[,1:427], comb_six$symbiont,  sim.method = "euclidean", p.adjust.m = "fdr", perm = 999)

# compare CC7-APO vs H2-APO
apo_animal=subset(asv.n, symbiont == "APO" & host != "None")
asv_adonis=adonis(apo_animal[,1:427]~ apo_animal$host, method = "euclidean")


####################################################
#####     Comparison 3: between C, CA, CB     ######
####################################################

CAB=subset(asv.n,host == "CC7")

bdisper_CAB=betadisper(dist(CAB[,1:427]), group = CAB$symbiont, sqrt.dist = T) # dist: euclidean distance; veg: bray-curtis distance
boxplot(bdisper_CAB)

adonis(CAB[,1:427]~ CAB$symbiont, method = "euclidean" )
pairwise.adonis(CAB[,1:427], CAB$symbiont,  sim.method = "euclidean", p.adjust.m = "fdr", perm = 999)

####################################################
#####      Comparison 4: between H, HA, HB     #####
####################################################

HAB=subset(asv.n,host == "H2")

bdisper_HAB=betadisper(dist(HAB[,1:427]), group = HAB$symbiont, sqrt.dist = T) # dist: euclidean distance; veg: bray-curtis distance
boxplot(bdisper_HAB)

adonis(HAB[,1:427]~ HAB$symbiont, method = "euclidean" )
pairwise.adonis(HAB[,1:427], HAB$symbiont,  sim.method = "euclidean", p.adjust.m = "fdr", perm = 999)


##### beta-dispersion ANOVA ######
# overall
bp_overall <- data.frame(bdisper_overall$distances)
bp_overall$factor=met$factor[match(rownames(bp_overall), rownames(met))]
bp_overall$factor=as.factor(bp_overall$factor)
bp_overall_ANOVA <- aov(bp_overall$bdisper_overall.distances ~ bp_overall$factor, data = bp_overall)
shapiro.test(residuals(bp_overall_ANOVA)) # P= 0,846
summary(bp_overall_ANOVA)
TukeyHSD(bp_overall_ANOVA)


# APO vs A vs B
bp_comb6 <- data.frame(bdisper_comb_six$distances)
bp_comb6$symbiont=met$symbiont[match(rownames(bp_comb6), rownames(met))]
bp_comb6$symbiont=as.factor(bp_comb6$symbiont)
bp_comb6_ANOVA <- aov(bp_comb6$bdisper_comb_six.distances ~ bp_comb6$symbiont, data = bp_comb6)
shapiro.test(residuals(bp_comb6_ANOVA)) # P=0,830
summary(bp_comb6_ANOVA)
TukeyHSD(bp_comb6_ANOVA)


# CA vs CB vs C-APO
bp_CAB <- data.frame(bdisper_CAB$distances)
bp_CAB$symbiont=met$symbiont[match(rownames(bp_CAB), rownames(met))]
bp_CAB$symbiont=as.factor(bp_CAB$symbiont)
bp_CAB_ANOVA <- aov(bp_CAB$bdisper_CAB.distances ~ bp_CAB$symbiont, data = bp_CAB)
shapiro.test(residuals(bp_CAB_ANOVA)) # P=0.7115
summary(bp_CAB_ANOVA)
TukeyHSD(bp_CAB_ANOVA)


# HA vs HB vs H-APO
bp_HAB <- data.frame(bdisper_HAB$distances)
bp_HAB$symbiont=met$symbiont[match(rownames(bp_HAB), rownames(met))]
bp_HAB$symbiont=as.factor(bp_HAB$symbiont)
bp_HAB_ANOVA <- aov(bp_HAB$bdisper_HAB.distances ~ bp_HAB$symbiont, data = bp_HAB)
shapiro.test(residuals(bp_HAB_ANOVA)) # P=0.751
summary(bp_HAB_ANOVA)
TukeyHSD(bp_HAB_ANOVA)




##############################       Old Analysis    ################################
asv=read.table("/Users/nanxiang/Desktop/PhD/WP3_Aiptasia_nirS/nirS_analysis/R_output/ASV_Tax_noConta.txt", header = TRUE, sep ='\t')
met=read.table("/Users/nanxiang/Desktop/PhD/WP3_Aiptasia_nirS/nirS_analysis/R_input/Aip_nirS_metadata.txt", header = TRUE, row.names = 1, sep ='\t')
asv$CA3=NULL # remove CA3, CA5
asv$CA5=NULL

#asv.n=as.data.frame(t(sweep(asv[, 1:95],2,colSums(asv[, 1:95]),`/`))) # relative abundances
asv.n=as.data.frame(t(apply(asv[, 1:41],2,clr))) # clr-transformed counts
asv.n$host=met$host[match(rownames(asv.n), rownames(met))]
asv.n$symbiont=met$symbiont[match(rownames(asv.n), rownames(met))]
asv.n$other=met$other[match(rownames(asv.n), rownames(met))]

##overal model
asv.n[is.na(asv.n)] <- "None"
asv_adonis=adonis(asv.n[,1:427]~ asv.n$host * asv.n$symbiont * asv.n$other, method = "euclidean")
asv_adonis_df=as.data.frame(asv_adonis[["aov.tab"]])
# write.table(asv_adonis_df, "outputs/lightPollution_overal_adonis_ASVs.txt", sep = "\t", row.names = T, quote = F)

## plot
asv_adonis_df$Factor=gsub(".*\\$", "", rownames(asv_adonis_df))
asv_adonis_df[4,7]="host-symbiont"
asv_adonis_df$Factor=factor(asv_adonis_df$Factor, levels= c("host", "symbiont", "other","host-symbiont"))

b=ggplot() +geom_bar(aes(y = R2, x = reorder(Factor,R2)), data = asv_adonis_df[1:4,], stat="identity", fill="black") +
  labs(x="", y="PERMANOVA R2") +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "none") +
  annotate("text", x = 1, y = 0.014, label = "*") +
  annotate("text", x =2, y = 0.016, label = "**") +
  annotate("text", x = 3, y = 0.041, label = "***") +
  annotate("text", x =4, y = 0.085, label = "***")

##### Comparison 1: between hosts regardless of algal identity #####
asv.n.host=subset(asv.n, host %in% c("CC7", "H2"))
adonis(asv.n.host[,1:427]~ asv.n.host$host, method = "euclidean" )
pairwise.adonis(asv.n.host[,1:427], asv.n.host$host,  sim.method = "euclidean", p.adjust.m = "fdr", perm = 999)


##### Comparison 2: between algae regardless of host identity #####
asv.n.algae=subset(asv.n, host %in% c("CC7", "H2") & symbiont %in% c("SSA01", "SSB01"))
adonis(asv.n.algae[,1:427]~ asv.n.algae$symbiont, method = "euclidean" )

##### Comparison 3: between C, CA, CB #####
CAB=subset(asv.n,host == "CC7")
adonis(CAB[,1:427]~ CAB$symbiont, method = "euclidean" )

##### Comparison 4: between H, HA, HB #####
HAB=subset(asv.n,host == "H2")
adonis(HAB[,1:427]~ HAB$symbiont, method = "euclidean" )

##### Comparison 5: comparing C, CA, CB pairwise #####
CAB_df=pairwise.adonis(CAB[,1:427], CAB$symbiont,  sim.method = "euclidean", p.adjust.m = "fdr", perm = 999)

##### Comparison 6: comparing H, HA, HB pairwise #####
HAB_df=pairwise.adonis(HAB[,1:427], HAB$symbiont,  sim.method = "euclidean", p.adjust.m = "fdr", perm = 999)

##### Comparison 7: comparing between cultures A and b #####
adonis(AB[,1:427]~ AB$symbiont,
       method = "euclidean")
anosim(AB[,1:427], AB$symbiont, permutations = 999, distance = "euclidean", strata = NULL,
       parallel = getOption("mc.cores"))
