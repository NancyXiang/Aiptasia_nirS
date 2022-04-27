### ### ### ### ### ### ### ### ###
###   nirS  Tax Assignment      ###
### ### ### ### ### ### ### ### ###

# use nirS Database excluding uncultured organisms

# limited unclutured bacteria to the genus level

# fungene_nirS_prot_dp_excludingUnc_3504.fasta   #3,504 seqs

# blastp against nirS database

# create database
makeblastdb -in fungene_nirS_prot_dp_excludingUnc_3504.fasta -parse_seqids -dbtype prot
# run Blast using the created DB  (database option:  -db)
blastp -query ~/project/nirS/nirS_tax_assignment/dada2_blastx_length_translated_nirS.fasta -db fungene_nirS_prot_dp_excludingUnc_3504.fasta -outfmt 6 -evalue 1e-30 >nirS_tax_fungene_subset_bp

scp xiangn@XXX:/home/xiangn/database/nirS/nirS_tax_fungene_subset_bp .

### subset the best hit for output from blastp in  R #### "filter_besthit_blastp.Rdata".

cut -f 2 blastp_fg_subset_output.txt >acn_nirS_bp_fg_subset_550 #NCBI prot accession
screen -S tax_fungene_subset_550

# step1 protAccession_to_TaxID
cat acn_nirS_bp_fg_subset_550 | while read line ; do efetch -db protein -id $line -format docsum | xtract -pattern DocumentSummary -element TaxId >> acn_nirS_bp_fg_subset_550_taxIDs; done

# step2 TaxID_to_Tax Info
source activate ete3
cat acn_nirS_bp_fg_subset_550_taxIDs | while read line ; do ete3 ncbiquery --search $line --info >> acn_nirS_bp_fg_subset_550_taxInfo; done

### formatted
cut -f1 -d ";" blastp_fg_subset_output.txt>seq_ID_fg_subset_550
sed '/#/d' acn_nirS_bp_fg_subset_550_taxInfo | cut -f 4 > tax_blastp_fg_subset_550
paste seq_ID_fg_subset_550 tax_blastp_fg_subset_550>output_nirS_fg_subset_tax_550

scp xiangn@XXX:/home/xiangn/project/nirS/nirS_tax_assignment/output_nirS_fg_subset_tax_550 /Users/nanxiang/Desktop/PhD/WP3_Aiptasia_nirS/nirS_analysis/bash_output/nirS_tax/

## organize output file in R ##
