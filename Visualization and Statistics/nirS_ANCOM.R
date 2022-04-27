#######################
###  ANCOM analysis ###
#######################

## Differential abundance analysis between samples

library(ANCOMBC)
library(phyloseq)
library("ape")
library(tibble)
library(ComplexHeatmap)

setwd("~/WP3_Aiptasia_nirS/nirS_analysis/")

###################
###  ASV level ###
##################

met=read.table("~/WP3_Aiptasia_nirS/nirS_analysis/R_input/Aip_nirS_metadata.txt", header = TRUE, row.names = 1, sep ='\t')

asv=read.table("~/WP3_Aiptasia_nirS/nirS_analysis/R_output/ASV_Tax_noConta_formated.txt", header = TRUE, row.names = 1, sep ='\t')

sam.t= sample_data(data.frame(met))
sam.t[is.na(sam.t)] <- "None"
sam.t$comb=paste(sam.t$host, sam.t$symbiont, sam.t$other, sep = "_")
sam.t$group1=ifelse(sam.t$comb %in% c("CC7_None_None","H2_None_None"), "APO", ifelse(sam.t$comb == "None_None_food", "Food",  ifelse(sam.t$symbiont == "SSA01", "SSA01", "SSB01" )))
row_names_to_remove<-c("NEC1","NEC2","NEC3","NPC1","NPC2", "CA3", "CA5")
sam.t=sam.t[!(row.names(sam.t) %in% row_names_to_remove),]

otu.t= otu_table(as.matrix(asv[, 1:43]), taxa_are_rows=TRUE)
tax.t= tax_table(as.matrix(asv[, 45:ncol(asv)]))
phy.all= phyloseq(otu.t, tax.t,  sam.t)

########## Comparison 1, between host CC7 and H2, regardless of the algal identity ##########
com_host=subset_samples(phy.all, !host == "None" )
res1_com_host=ancombc(phyloseq=com_host,formula="host",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "host",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res1_com_host_df=data.frame(res1_com_host[["res"]])
colnames(res1_com_host_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res1_com_host_df_sig=subset(res1_com_host_df, Diff_abundant == "TRUE")#[,c(1,2,6,8,10)]
res1_com_host_df_sig$Diff_more_abundant=ifelse(res1_com_host_df_sig$W<0, "CC7", "H2")
res1_com_host_df_sig$Taxa=paste(asv$Phylum, asv$Family, asv$Genus,sep = "_")[match(rownames(res1_com_host_df_sig),rownames(asv))]
message("Number of DA ASVs: ", nrow(res1_com_host_df_sig), "\nNumber of DA ASVs enriched in CC7: ", nrow(subset(res1_com_host_df_sig, Diff_more_abundant == "CC7" )), "\nNumber of DA ASVs enriched in H2: ", nrow(subset(res1_com_host_df_sig, Diff_more_abundant == "H2" )))
# Number of DA ASVs enriched in CC7: 4; Number of DA ASVs enriched in H2: 0
write.table(res1_com_host_df_sig,  "R_output/ANCOMBC_ASVs_Comp1_host.txt", sep = "\t", quote = F, row.names = T )


########## Comparison 2, between algae, regardless of the host ##########
# APO and SSA01
apo_ssa=subset_samples(phy.all, !group1 == "Food" & group1 %in% c( "APO", "SSA01"))
res1_apo_ssa=ancombc(phyloseq=apo_ssa,formula="group1",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "group1",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res1_apo_ssa_df=data.frame(res1_apo_ssa[["res"]])
colnames(res1_apo_ssa_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res1_apo_ssa_df_sig=subset(res1_apo_ssa_df, Diff_abundant == "TRUE")#[,c(1,2,6,8,10)]
res1_apo_ssa_df_sig$Diff_more_abundant=ifelse(res1_apo_ssa_df_sig$W<0, "APO", "A")
res1_apo_ssa_df_sig$Taxa=paste(asv$Phylum, asv$Family, asv$Genus,sep = "_")[match(rownames(res1_apo_ssa_df_sig),rownames(asv))]
message("Number of DA ASVs: ", nrow(res1_apo_ssa_df_sig), "\nNumber of DA ASVs enriched in APO: ", nrow(subset(res1_apo_ssa_df_sig, Diff_more_abundant == "APO" )), "\nNumber of DA ASVs enriched in A: ", nrow(subset(res1_apo_ssa_df_sig, Diff_more_abundant == "A" )))
# Number of DA ASVs enriched in APO: 0;  Number of DA ASVs enriched in A: 2
write.table(res1_apo_ssa_df_sig,  "R_output/ANCOMBC_ASVs_Comp2_host_APO_A.txt", sep = "\t", quote = F, row.names = T )

# APO and SSB01
apo_ssb=subset_samples(phy.all, !group1 == "Food" & group1 %in% c( "APO", "SSB01"))
res1_apo_ssb=ancombc(phyloseq=apo_ssb,formula="group1",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "group1",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res1_apo_ssb_df=data.frame(res1_apo_ssb[["res"]])
colnames(res1_apo_ssb_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res1_apo_ssb_df_sig=subset(res1_apo_ssb_df, Diff_abundant == "TRUE")#[,c(1,2,6,8,10)]
res1_apo_ssb_df_sig$Diff_more_abundant=ifelse(res1_apo_ssb_df_sig$W<0, "APO", "B")
res1_apo_ssb_df_sig$Taxa=paste(asv$Phylum, asv$Family, asv$Genus,sep = "_")[match(rownames(res1_apo_ssb_df_sig),rownames(asv))]
message("Number of DA ASVs: ", nrow(res1_apo_ssb_df_sig), "\nNumber of DA ASVs enriched in APO: ", nrow(subset(res1_apo_ssb_df_sig, Diff_more_abundant == "APO" )), "\nNumber of DA ASVs enriched in B: ", nrow(subset(res1_apo_ssb_df_sig, Diff_more_abundant == "B" )))
# Number of DA ASVs enriched in APO: 0;  Number of DA ASVs enriched in A: 2
write.table(res1_apo_ssb_df_sig,  "R_output/ANCOMBC_ASVs_Comp2_host_APO_B.txt", sep = "\t", quote = F, row.names = T )

# SSA01 and SSB01
ssa_ssb=subset_samples(phy.all, !group1 == "Food" & group1 %in% c( "SSA01", "SSB01"))
res1_ssa_ssb=ancombc(phyloseq=ssa_ssb,formula="group1",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "group1",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res1_ssa_ssb_df=data.frame(res1_ssa_ssb[["res"]])
colnames(res1_ssa_ssb_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res1_ssa_ssb_df_sig=subset(res1_ssa_ssb_df, Diff_abundant == "TRUE")#[,c(1,2,6,8,10)]
res1_ssa_ssb_df_sig$Diff_more_abundant=ifelse(res1_ssa_ssb_df_sig$W<0, "A", "B")
res1_ssa_ssb_df_sig$Taxa=paste(asv$Phylum, asv$Family, asv$Genus,sep = "_")[match(rownames(res1_ssa_ssb_df_sig),rownames(asv))]
message("Number of DA ASVs: ", nrow(res1_ssa_ssb_df_sig), "\nNumber of DA ASVs enriched in SSA01: ", nrow(subset(res1_ssa_ssb_df_sig, Diff_more_abundant == "A" )), "\nNumber of DA ASVs enriched in SSB01: ", nrow(subset(res1_ssa_ssb_df_sig, Diff_more_abundant == "B" )))
# Number of DA ASVs enriched in APO: 0;  Number of DA ASVs enriched in A: 2
write.table(res1_ssa_ssb_df_sig,  "R_output/ANCOMBC_ASVs_Comp2_host_A_B.txt", sep = "\t", quote = F, row.names = T )


########## Comparison 3  subset CC7 -APO, A, B ##########
# C and CA
cc7_phy_ssa=subset_samples(phy.all, host == "CC7" & comb %in% c("CC7_None_None", "CC7_SSA01_None") ) #subset C-APO and C-A
res1_cc7_ssa=ancombc(phyloseq=cc7_phy_ssa,formula="comb",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "comb",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res1_cc7_ssa_df=data.frame(res1_cc7_ssa[["res"]])

colnames(res1_cc7_ssa_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res1_cc7_ssa_df_sig=subset(res1_cc7_ssa_df, Diff_abundant == "TRUE")# different taxa between C-APO and CA
res1_cc7_ssa_df_sig$Diff_more_abundant=ifelse(res1_cc7_ssa_df_sig$W<0, "APO", "A")

res1_cc7_ssa_df_sig$Taxa=paste(asv$Phylum, asv$Family, asv$Genus,sep = "_")[match(rownames(res1_cc7_ssa_df_sig),rownames(asv))]
message("\nNumber of DA ASVs: ", nrow(res1_cc7_ssa_df_sig), "\nNumber of DA ASVs enriched in CC7-APO: ", nrow(subset(res1_cc7_ssa_df_sig, Diff_more_abundant == "APO" )), "\nNumber of DA ASVs enriched in CC7 + A: ", nrow(subset(res1_cc7_ssa_df_sig, Diff_more_abundant == "A" )))
# Number of DA ASVs enriched in CC7-APO: 2 ; Number of DA ASVs enriched in CC7 + A: 3
write.table(res1_cc7_ssa_df_sig,  "R_output/ANCOMBC_ASVs_Comp3_C_CA.txt", sep = "\t", quote = F, row.names = T )


# C and CB
cc7_phy_ssb=subset_samples(phy.all, host == "CC7" & comb %in% c("CC7_None_None", "CC7_SSB01_None") )
res1_cc7_ssb=ancombc(phyloseq=cc7_phy_ssb,formula="comb",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "comb",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res1_cc7_ssb_df=data.frame(res1_cc7_ssb[["res"]])
colnames(res1_cc7_ssb_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res1_cc7_ssb_df_sig=subset(res1_cc7_ssb_df, Diff_abundant == "TRUE")#[,c(1,2,6,8,10)]
res1_cc7_ssb_df_sig$Diff_more_abundant=ifelse(res1_cc7_ssb_df_sig$W<0, "APO", "B")
res1_cc7_ssb_df_sig$Taxa=paste(asv$Phylum, asv$Family, asv$Genus,sep = "_")[match(rownames(res1_cc7_ssb_df_sig),rownames(asv))]
message("Number of DA ASVs: ", nrow(res1_cc7_ssb_df_sig), "\nNumber of DA ASVs enriched in CC7 APO: ", nrow(subset(res1_cc7_ssb_df_sig, Diff_more_abundant == "APO" )), "\nNumber of DA ASVs enriched in CC7 + B: ", nrow(subset(res1_cc7_ssb_df_sig, Diff_more_abundant == "B" )))
# Number of DA ASVs enriched in CC7 APO: 1; Number of DA ASVs enriched in CC7 + B: 7
write.table(res1_cc7_ssb_df_sig,  "R_output/ANCOMBC_ASVs_Comp3_C_CB.txt", sep = "\t", quote = F, row.names = T )

# CA and CB
cc7_phy_CAB=subset_samples(phy.all, host == "CC7" & comb %in% c("CC7_SSA01_None", "CC7_SSB01_None") )
res1_CAB=ancombc(phyloseq=cc7_phy_CAB,formula="comb",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "comb",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res1_CAB_df=data.frame(res1_CAB[["res"]])
colnames(res1_CAB_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res1_CAB_df_sig=subset(res1_CAB_df, Diff_abundant == "TRUE")#[,c(1,2,6,8,10)]
res1_CAB_df_sig$Diff_more_abundant=ifelse(res1_CAB_df_sig$W<0, "A", "B")
res1_CAB_df_sig$Taxa=paste(asv$Phylum, asv$Family, asv$Genus,sep = "_")[match(rownames(res1_CAB_df_sig),rownames(asv))]
message("Number of DA ASVs: ", nrow(res1_CAB_df_sig), "\nNumber of DA ASVs enriched in CA: ", nrow(subset(res1_CAB_df_sig, Diff_more_abundant == "A" )), "\nNumber of DA ASVs enriched in CB: ", nrow(subset(res1_CAB_df_sig, Diff_more_abundant == "B" )))
# Number of DA ASVs enriched in CC7 APO: 1; Number of DA ASVs enriched in CC7 + B: 7
write.table(res1_CAB_df_sig,  "R_output/ANCOMBC_ASVs_Comp3_CA_CB.txt", sep = "\t", quote = F, row.names = T )


########## Comparison 4  H2 -APO, A, B ##########
# H and HA
H2_phy_ssa=subset_samples(phy.all, host == "H2" & comb %in% c("H2_None_None", "H2_SSA01_None") ) #subset
res1_H2_ssa=ancombc(phyloseq=H2_phy_ssa,formula="comb",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "comb",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res1_H2_ssa_df=data.frame(res1_H2_ssa[["res"]])
colnames(res1_H2_ssa_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res1_H2_ssa_df_sig=subset(res1_H2_ssa_df, Diff_abundant == "TRUE")# different taxa between C-APO and CA
res1_H2_ssa_df_sig$Diff_more_abundant=ifelse(res1_H2_ssa_df_sig$W<0, "APO", "A")

res1_H2_ssa_df_sig$Taxa=paste(asv$Phylum, asv$Family, asv$Genus,sep = "_")[match(rownames(res1_H2_ssa_df_sig),rownames(asv))]
message("\nNumber of DA ASVs: ", nrow(res1_H2_ssa_df_sig), "\nNumber of DA ASVs enriched in H2-APO: ", nrow(subset(res1_H2_ssa_df_sig, Diff_more_abundant == "APO" )), "\nNumber of DA ASVs enriched in H2 + A: ", nrow(subset(res1_H2_ssa_df_sig, Diff_more_abundant == "A" )))
# Number of DA ASVs enriched in H2-APO: 2; Number of DA ASVs enriched in H2 + A: 4
write.table(res1_H2_ssa_df_sig,  "R_output/ANCOMBC_ASVs_Comp4_H_HA.txt", sep = "\t", quote = F, row.names = T )


# H and HB
H2_phy_ssb=subset_samples(phy.all, host == "H2" & comb %in% c("H2_None_None", "H2_SSB01_None") )
res1_H2_ssb=ancombc(phyloseq=H2_phy_ssb,formula="comb",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "comb",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res1_H2_ssb_df=data.frame(res1_H2_ssb[["res"]])
colnames(res1_H2_ssb_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res1_H2_ssb_df_sig=subset(res1_H2_ssb_df, Diff_abundant == "TRUE")#[,c(1,2,6,8,10)]
res1_H2_ssb_df_sig$Diff_more_abundant=ifelse(res1_H2_ssb_df_sig$W<0, "APO", "B")
res1_H2_ssb_df_sig$Taxa=paste(asv$Phylum, asv$Family, asv$Genus,sep = "_")[match(rownames(res1_H2_ssb_df_sig),rownames(asv))]
message("Number of DA ASVs: ", nrow(res1_H2_ssb_df_sig), "\nNumber of DA ASVs enriched in H2 APO: ", nrow(subset(res1_H2_ssb_df_sig, Diff_more_abundant == "APO" )), "\nNumber of DA ASVs enriched in H2 + B: ", nrow(subset(res1_H2_ssb_df_sig, Diff_more_abundant == "B" )))
# Number of DA ASVs enriched in H2 APO: 2; Number of DA ASVs enriched in H2 + B: 2
write.table(res1_H2_ssb_df_sig,  "R_output/ANCOMBC_ASVs_Comp4_H_HB.txt", sep = "\t", quote = F, row.names = T )


# HA and HB
H2_HAB=subset_samples(phy.all, host == "H2" & comb %in% c("H2_SSA01_None", "H2_SSB01_None") )
res1_H2_HAB=ancombc(phyloseq=H2_HAB,formula="comb",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "comb",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res1_H2_HAB_df=data.frame(res1_H2_HAB[["res"]])
colnames(res1_H2_HAB_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res1_H2_HAB_df_sig=subset(res1_H2_HAB_df, Diff_abundant == "TRUE")#[,c(1,2,6,8,10)]
res1_H2_HAB_df_sig$Diff_more_abundant=ifelse(res1_H2_HAB_df_sig$W<0, "A", "B")
res1_H2_HAB_df_sig$Taxa=paste(asv$Phylum, asv$Family, asv$Genus,sep = "_")[match(rownames(res1_H2_HAB_df_sig),rownames(asv))]
message("Number of DA ASVs: ", nrow(res1_H2_HAB_df_sig), "\nNumber of DA ASVs enriched in HA: ", nrow(subset(res1_H2_HAB_df_sig, Diff_more_abundant == "A" )), "\nNumber of DA ASVs enriched in HB: ", nrow(subset(res1_H2_HAB_df_sig, Diff_more_abundant == "B" )))
# ASVs enriched in HA: 1; Number of DA ASVs enriched in HB: 1;
write.table(res1_H2_HAB_df_sig,  "R_output/ANCOMBC_ASVs_Comp4_HA_HB.txt", sep = "\t", quote = F, row.names = T )

