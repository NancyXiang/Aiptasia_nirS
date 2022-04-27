### ### ### ### ### ### ### ### ###
###   nirS  prot database prep  ###
### ### ### ### ### ### ### ### ###


# filter database excluding unculture organsims

setwd("/Users/nanxiang/Desktop/PhD/WP3_Aiptasia_nirS/nirS_analysis")
save.image("database_nirS_excl_uncul.Rdata")

library(seqinr)
library(dplyr)

fasta<-as.data.frame(phylotools::read.fasta("/Users/nanxiang/Desktop/PhD/WP3_Aiptasia_nirS/nirS_analysis/database/fungene_repository_nirS_prot_ref_25330.fasta"))

fasta$species<-gsub(".*organism=", "", fasta$seq.name)

fasta$genus<-sub(" .*", "", fasta$species)

fasta$genus= as.factor(fasta$genus)
levels(fasta$genus)

# subset by removing genus = "uncultured"
fasta_subset=subset(fasta, fasta$genus != "uncultured")

# check specie numbers
fasta_subset$species_tax=gsub(",.*", "", fasta_subset$species)
species_final=fasta_subset %>% group_by(species_tax) %>%  sample_n(1)

write.fasta(as.list(fasta_subset$seq.text), fasta_subset$seq.name, "/Users/nanxiang/Desktop/PhD/WP3_Aiptasia_nirS/nirS_analysis/database/nirS_prot_dp_excludingUnc_3504.fasta", open = "w", nbchar = 60, as.string = FALSE)
