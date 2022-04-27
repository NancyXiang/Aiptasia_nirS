### ### ### ### ### ### ### ###
###   filter positive nirS  ###
### ### ### ### ### ### ### ###

setwd("~/Desktop/PhD/WP3_Aiptasia_nirS/nirS_analysis/")
save.image("filter_positive_nirS.Rdata")
load("filter_positive_nirS.Rdata")

#####################################
############ after DADA2 ############
#####################################

sum(pre_asv$sum)
# in total 7,414,102 reads represented by 1104 ASVs
pre_asv=read.table("/Users/nanxiang/Desktop/PhD/WP3_Aiptasia_nirS/nirS_analysis/bash_output/nirS_full/NirS_full_ASV_table.txt", header = T)
pre_asv$Old_ID=paste("sq",1:nrow(pre_asv), sep="" )
boxplot(nchar(pre_asv$Sequence))
hist(nchar(pre_asv$Sequence))
sum(pre_asv$sum)
# 7,414,102 reads represented by 1,104 ASVs

########## after Uniprot blastx ############

blasx_uniprot=read.table("/Users/nanxiang/Desktop/PhD/WP3_Aiptasia_nirS/nirS_analysis/bash_output/nirS_full/nirS_uniprot_positive_deduplicated.txt", header = T, sep = "\t")
nirS_uniport_asv=subset(pre_asv, pre_asv$Old_ID %in% blasx_uniprot$nirS)
sum(nirS_uniport_asv$sum)
# in total 6,907,743 reads represented by 616 ASVs
sum(nirS_uniport_asv$sum)/sum(pre_asv$sum)
# keep 93,17% of sequences
hist(nchar(nirS_uniport_asv$Sequence))


# plot seq_length distribution
attach(mtcars)
par(mfrow=c(3,1))
hist(nchar(pre_asv$Sequence))
hist(nchar(nirS_uniport_asv$Sequence))

# plot hist in a narrow window
attach(mtcars)
par(mfrow=c(3,1))
hist(nchar(nirS_uniport_asv$Sequence))
hist(nchar(nirS_uniport_asv$Sequence),ylim=c(0,20))
hist(nchar(nirS_uniport_asv$Sequence),xlim=c(200,270))

### histogram
library(ggplot2)
ggplot(nirS_uniprot_length_asv, aes(nchar(Sequence))) +            # ggplot2 histogram with default bins
  geom_histogram(bins = 20) +
  scale_x_continuous(breaks = c(210,215,220,225,230,235,240,245,250,255,260,265,270,275,280,285,290,295,300))


###########  subset nirS fasta after blastx #############
# Reading fasta
library(devtools)
fasta<-as.data.frame(phylotools::read.fasta("/Users/nanxiang/Desktop/PhD/WP3_Aiptasia_nirS/nirS_analysis/bash_output/nirS_full/dada2_nirS_full.fasta"))

library(reshape2)
fasta_new = cbind(fasta,colsplit(fasta$seq.name,";", c("ASV", "size")))

blasx_uniprot2 = blasx_uniprot[,-1]
rownames(blasx_uniprot2) = blasx_uniprot[,1]
blasx_uniprot2$ASV = blasx_uniprot$nirS
blasx_uniprot = blasx_uniprot2
rm(blasx_uniprot2)

fasta_uniprot_filtered=subset(fasta_new, fasta_new$ASV %in% blasx_uniprot$ASV)
fasta_uniprot_filtered=subset(fasta_uniprot_filtered, select = -c(ASV,size))

install.packages('seqinr')
library(seqinr)
write.fasta(as.list(fasta_uniprot_filtered$seq.text), fasta_uniprot_filtered$seq.name, "~/Desktop/fasta_uniprot_filtered.fasta", open = "w", nbchar = 60, as.string = FALSE)

###############################################################
###################  after length  selection ##################
###############################################################

###       run in Bash       ###

# filter uniport sequence by length 220-240 bp
# bioawk -c fastx '{ if(length($seq) > 220) { print ">"$name; print $seq }}' dada2_nirS_full_uniprot_filtered.fasta >output_above220.fasta
# bioawk -c fastx '{ if(length($seq) < 240) { print ">"$name; print $seq }}' output_above220.fasta >dada2_nirS_uniprot_filtered_220_240.fasta
# 565 ASVs

### subset fasta by length  ###
library(devtools)
nirS_fasta_filter_uniprot_length<-as.data.frame(phylotools::read.fasta("/Users/nanxiang/Desktop/PhD/WP3_Aiptasia_nirS/nirS_analysis/bash_output/nirS_full/dada2_nirS_uniprot_filtered_220_240.fasta"))

library(reshape2)
nirS_fasta_filter_uniprot_length = cbind(nirS_fasta_filter_uniprot_length,colsplit(nirS_fasta_filter_uniprot_length$seq.name,";", c("ASV", "size")))
nirS_uniprot_length_asv=subset(pre_asv, pre_asv$Old_ID %in% nirS_fasta_filter_uniprot_length$ASV)

sum(nirS_uniprot_length_asv$sum)
# in total 6,907,542 reads represented by 565 ASVs
sum(nirS_uniprot_length_asv$sum)/sum(pre_asv$sum)
# keep 93,17% of sequences
hist(nchar(nirS_uniprot_length_asv$Sequence))

nirS_fasta_filter_uniprot_length=subset(nirS_fasta_filter_uniprot_length, select = -c(ASV,size))
nirS_uniprot_length_asv$length = nchar(nirS_uniprot_length_asv$Sequence)
write.fasta(as.list(nirS_uniprot_length_asv$Sequence),
            paste(nirS_uniprot_length_asv$Old_ID, ",
                  length=",nirS_uniprot_length_asv$length,",
                  size=" ,nirS_uniprot_length_asv$sum ,sep= ""),
            "after_blastx_length_selection_nirS.fasta",
            open = "w", nbchar = 60, as.string = FALSE)


#########################################################################
###################  after translation to correct ORF  ##################
#########################################################################

after_translation=read.table("/Users/nanxiang/Desktop/PhD/WP3_Aiptasia_nirS/nirS_analysis/Aiptasia_nirS_seq_fasta/ASV_blastx_length_translated.txt", header = F, sep = "\t")
nirS_after_translation_asv=subset(pre_asv, pre_asv$Old_ID %in% after_translation$V1)
sum(nirS_after_translation_asv$sum)

sum(nirS_after_translation_asv$sum)/sum(pre_asv$sum)

hist(nchar(nirS_after_translation_asv$Sequence))
