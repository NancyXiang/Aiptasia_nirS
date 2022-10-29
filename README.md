# Denitrifiers in the Cnidarian Holobiont

This repository contains the R and Bash scripts for nirS amplicon sequence analysis, data visualization, and statistics. 
Code Developper: Anny Cárdenas, Nan Xiang.

Nan Xiang, Nils Rädecker, Claudia Pogoreutz, Anny Cárdenas, Anders Meibom, Christian Wild, Astrid Gärdes, and Christian R. Voolstra. Presence of algal symbionts affects denitrifying bacterial communities in the sea anemone Aiptasia coral model. ISME COMMUN. 2, 105 (2022). https://doi.org/10.1038/s43705-022-00190-9.

Raw sequencing data are deposited in the NCBI Sequence Read Archive (SRA) under BioProject PRJNA836569. 

# Workflow

## Part A: nirS Sequence Analysis

A1: nirS Amplicon Sequence variants (ASV) were inferred using dada2 (link to github) using the script (nirS_Dada2.R).

A2: Quality checks: a) removal of sequences not annotated to “nirS” protein after Blastx and removal of sequences out of 220-240 bp using the script (nirS_Filter_positives.R). 

A3: ASVs were translated to the nirS protein in a correct open reading frame (ORF) using TranslatorX, ORFfinder in NCBI. 


## Part B: nirS database and Taxonomic Assignment 

B1: nirS database was prepared using the script (nirS_Database_accession_to_tax.sh) and (nirS_Prot_database_prep.R) based on the input file (fungene_repository_nirS_prot_ref_25330.fasta). 

B2: Taxonomic assignment for nirS sequences were done using the script (nirS_Tax_assign.sh). 


## Part C: Data visualization and Statistics 
Scripts from https://github.com/ajcardenasb?tab=repositories

C1: Removal of putative contaminant ASVs using the script (nirS_ContaRemoval.R).

C2: Alpha diversity and statistical analysis were done using the script (nirS_AlphaDiversity.R).

C3: Betadiversity ordination was done using the script (nirS_BetaDiversity.R). 

C4: Barplots showing most abundant bacterial ASVs were visualized using the script (nirS_Barplot.R).

C5: PERMANOVA and Betadispersion were done using the script (nirS_PERMANOVA.R)

C6: ASV enrichment analysis was done using the script (nirS_ANCOM.R).

C7: Heatmap showing differential abundances of ASVs was plotted using the script (nirS_Heatmap.R).
