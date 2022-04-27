###################################################################
#####         Identifying and Removing contaminant ASVs       #####
###################################################################

library(tidyverse)

setwd("/Users/nanxiang/Desktop/PhD/WP3_Aiptasia_nirS/nirS_analysis")

##########################################################################
##### Step One: Determine Contamination from negative DNA extraction #####
##########################################################################

###     use  preASV table as input (DADA2 output table)
# Read Table
pre_asv=read.table("~/nirS_dada2/NirS_full_ASV_table.txt", header = T)
pre_asv$Old_ID=paste("sq",1:nrow(pre_asv), sep="" )
names(pre_asv)
pre_asv2 <- pre_asv[,-1]
rownames(pre_asv2) <- pre_asv[,1]


# Remove samples with low read number
names(pre_asv2)
pre_asv2.n=apply(pre_asv2[,c(1:48)], 2, as.numeric) # put into numeric
pre_asv2.o=pre_asv2.n[, colSums(pre_asv2.n) > 1000] # subset with sample > 1000 reads
dim(pre_asv2.o)
rownames(pre_asv2.o)=rownames(pre_asv2)
message(ncol(pre_asv2.o)," samples with > 1000 reads were retained out of ", ncol(pre_asv2.o), " total samples")


# Identify and removing contaminant ASVs from raw data
pre_asv2.r=as.data.frame(sweep(pre_asv2.o,2,colSums(pre_asv2.o),`/`)) # proportion of per ASV for per sample
pre_asv2.r=transform(pre_asv2.r,  Sum = rowSums(pre_asv2.r[,1:ncol(pre_asv2.r)])) # Sum for per ASV, add new column "Sum"
names(pre_asv2.r)

pre_asv2.r=transform(pre_asv2.r,  SumNegs = pre_asv2.r[,c(44:48)]) # define negative controls here!!
names(pre_asv2.r)

pre_asv2.r=transform(pre_asv2.r,  contaFactor1=(pre_asv2.r$SumNegs.NEC1/pre_asv2.r$Sum)*100) # proportion of NC_ASV/sum_ASV
pre_asv2.r=transform(pre_asv2.r,  contaFactor2=(pre_asv2.r$SumNegs.NEC2/pre_asv2.r$Sum)*100) #
pre_asv2.r=transform(pre_asv2.r,  contaFactor3=(pre_asv2.r$SumNegs.NEC3/pre_asv2.r$Sum)*100) #
pre_asv2.r=transform(pre_asv2.r,  contaFactor4=(pre_asv2.r$SumNegs.NPC1/pre_asv2.r$Sum)*100) #
pre_asv2.r=transform(pre_asv2.r,  contaFactor5=(pre_asv2.r$SumNegs.NPC2/pre_asv2.r$Sum)*100) #
names(pre_asv2.r)
rownames(pre_asv2.r)=rownames(pre_asv2)


# define Conta ASV using preASV table
pre_Conta=subset(pre_asv2.r, pre_asv2.r$contaFactor1 > 10 | pre_asv2.r$contaFactor2 > 10 | pre_asv2.r$contaFactor3 > 10 | pre_asv2.r$contaFactor4 > 10 | pre_asv2.r$contaFactor5 > 10)
pre_Conta$Family=pre_asv2$Family[match(rownames(pre_Conta), rownames(pre_asv2))]
message("Number of total ASVs: ", nrow(pre_asv2))
message("Number of identified contaminant ASVs removed from the analysis: ", length(rownames(pre_Conta)), "\n", pre_Conta$Family[1],"\n", pre_Conta$Family[2],"\n", pre_Conta$Family[3],"\n", pre_Conta$Family[4],"\n", pre_Conta$Family[5])


# check how much read from  conta ASV
pre_Conta$Old_ID=rownames(pre_Conta)
pre_Conta$Old_ID=gsub("preASV000", "sq",  pre_Conta$Old_ID)
pre_Conta$Old_ID=gsub("preASV00", "sq",  pre_Conta$Old_ID)
pre_Conta$Old_ID=gsub("preASV0", "sq",  pre_Conta$Old_ID)
pre_Conta$Old_ID=gsub("preASV", "sq",  pre_Conta$Old_ID)
pre_Conta$Old_ID


# subset a conta ASV table
asv_pre_conta=subset(pre_asv, pre_asv$Old_ID %in% pre_Conta$Old_ID)
sum(asv_pre_conta$sum)
sum(asv_pre_conta$sum)/sum(pre_asv$sum)

# subset a ASV table with contamination removed, using filtered 558_ASV table (very important!) (ASV table passed quality checks)

# read nirS filtered 558_ASV table
asv = read.table("/Users/nanxiang/Desktop/PhD/WP3_Aiptasia_nirS/nirS_analysis/tax_nirS/filtered_asv.txt", sep = "\t", header = T, row.names = 1)
met = read.table("/Users/nanxiang/Desktop/PhD/WP3_Aiptasia_nirS/nirS_analysis/R_input/Aip_nirS_metadata.txt", header = T, sep = "\t", row.names = 1) ## metadata with your sample varaibles
names(asv)
asv.n=apply(asv[,c(1:48)], 2, as.numeric) # put into numeric
asv.o=asv.n[, colSums(asv.n) > 1000] # subset with sample  >1000 reads
rownames(asv.o)=rownames(asv)


preasv.noConta=subset(asv.o, !rownames(asv.o) %in% rownames(pre_Conta))[,-c(44:48)] # define negative control
colnames(preasv.noConta)
dim(preasv.noConta) # 516 ASVs, 43 samples

preasv.noConta=as.data.frame(preasv.noConta)
preasv.noConta$sum = rowSums(preasv.noConta[1-43])
sum(preasv.noConta$sum)
sum(preasv.noConta$sum)/sum(pre_asv$sum)


##########################################################################
##### Step Two: Determine Contamination from genus "Cupriavidus"     #####
##########################################################################


## we used it as positive control for PCR development, presumbly there will be some contaminants of it
## it is known not associated with Aiptasia

### input file: noConta_ASV table, filtered conta by preASV, and subset from the filtered 558 ASV table
preasv.noConta$Old_ID=rownames(preasv.noConta)
preasv.noConta$Old_ID=gsub("preASV000", "sq",  preasv.noConta$Old_ID)
preasv.noConta$Old_ID=gsub("preASV00", "sq",  preasv.noConta$Old_ID)
preasv.noConta$Old_ID=gsub("preASV0", "sq",  preasv.noConta$Old_ID)
preasv.noConta$Old_ID=gsub("preASV", "sq",  preasv.noConta$Old_ID)
preasv.noConta$Old_ID
preasv.noConta = preasv.noConta %>% remove_rownames %>% column_to_rownames(var="Old_ID")

### read TAX table
Tax <- read.table("/Users/nanxiang/Desktop/PhD/WP3_Aiptasia_nirS/nirS_analysis/tax_nirS/nirS_ASV_Tax.txt", sep = "\t", header = T, row.names = 1, na.strings=c("","NA"))
Tax$sum_reads=NULL
Tax = separate(Tax, tax_blastp_fungene_550_without_uncultured, into = c("1", "2", "3","4", "5", "6","7", "8"), sep = ";", remove = F)
names(Tax)[names(Tax) == '6'] <- 'Genus'

### combine TAX and ASV
preasv.noConta.f=merge(preasv.noConta, Tax, by="row.names")
write.table(preasv.noConta.f, "/Users/nanxiang/Desktop/PhD/WP3_Aiptasia_nirS/nirS_analysis/R_output/ASV_Tax_noConta_preASVtxt",  quote = FALSE, row.names=F, sep = "\t") #define sample range

####   remove genus "Cupriavidus"  #####
nrow(subset(preasv.noConta.f, Genus == "Cupriavidus"))
subset(preasv.noConta.f, Genus == "Cupriavidus")$sum %>% sum ()
subset(preasv.noConta.f, Genus == "Cupriavidus")$Row.names

preasv.noConta$sum=NULL
str(preasv.noConta)

# sum reads for each ASV
preasv.noConta$sum=rowSums(preasv.noConta)
sum(preasv.noConta$sum)

# combine ASV and TAX together
preasv.noConta.f=merge(preasv.noConta, Tax, by="row.names")

# remove ASVs belong to genus "Cupriavidus"
Cupriavidus=subset(preasv.noConta.f, Genus == "Cupriavidus")
preasv.noConta.final=(subset(preasv.noConta.f, !(preasv.noConta.f$Row.names %in% Cupriavidus$Row.names)))

# formatted
row.names(preasv.noConta.final) = preasv.noConta.final$Row.names
preasv.noConta.final$Row.names=NULL
sum(preasv.noConta.final$sum)
sum(preasv.noConta.final$sum)/sum(pre_asv$sum)

### check the taxonomic information
final_ASV_TAX<-preasv.noConta.final[, c(44-53)]
final_ASV_TAX[is.na(final_ASV_TAX)] <- "unclassified"

### check the statistical information
final_ASV_TAX_stat <- aggregate(final_ASV_TAX$sum, by=list(Category=final_ASV_TAX$Genus), FUN=sum)
names(final_ASV_TAX_stat)[2] <- "tax_FunGene_noConta"
colnames(final_ASV_TAX_stat)[1]<-"Genus"
sum(final_ASV_TAX_stat$tax_FunGene_noConta) # sum reads
nrow(subset(final_ASV_TAX, Genus == "unclassified")) # 122 ASVs are unclassified bacteria
subset(final_ASV_TAX_stat, Genus == "unclassified")$tax_FunGene_noConta/sum(preasv.noConta.final$sum)

###  Output  final Tax, ASV, and stats Table
write.table(preasv.noConta.final, "/Users/nanxiang/Desktop/PhD/WP3_Aiptasia_nirS/nirS_analysis/R_output/ASV_Tax_noConta.txt",  quote = FALSE, row.names=F, sep = "\t")
message("Number of ASVs used in the analysis: ", length(rownames(preasv.noConta.final)))

# check stats at the genus level 
write.table(final_ASV_TAX_stat, "/Users/nanxiang/Desktop/PhD/WP3_Aiptasia_nirS/nirS_analysis/R_output/ASV_Tax_noConta_stat.txt",  quote = FALSE, row.names=F, sep = "\t")
