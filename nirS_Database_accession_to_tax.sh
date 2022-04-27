### ### ### ### ### ### ### ### ### ### ### ### ### ###
### NCBI Accesion number_to_Taxonomic information   ###
### ### ### ### ### ### ### ### ### ### ### ### ### ###


# extract ID
grep -o -E "^>\w+" nirS_prot_ref_25330.fasta | tr -d ">" >nirS_NCBI_accessionNumber


# install Entrez Direct (EDirect)
sh -c "$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"


# activate EDirect for this terminal session,
export PATH=${PATH}:${HOME}/edirect


# transfer one nucl accession number to full taxonomy
esearch -db assembly -query "GCA_000005845" | elink -target taxonomy | efetch -format native -mode xml | grep ScientificName | awk -F ">|<" 'BEGIN{ORS=", ";}{print $3;}'
# Escherichia coli str. K-12 substr. MG1655, cellular organisms, Bacteria, Proteobacteria, Gammaproteobacteria, Enterobacterales, Enterobacteriaceae, Escherichia, Escherichia coli, Escherichia coli K-12


# transfer one prot accession number to full taxonomy
esearch -db protein -query "AXS91626" | elink -target taxonomy | efetch -format native -mode xml | grep ScientificName | awk -F ">|<" 'BEGIN{ORS=", ";}{print $3;}'
# Pseudomonas aeruginosa, cellular organisms, Bacteria, Proteobacteria, Gammaproteobacteria, Pseudomonadales, Pseudomonadaceae, Pseudomonas, Pseudomonas aeruginosa group,


# a list of prot Accession to full_tax_info
# two-step commands

# step1 protAccession_to_TaxID
cat nirS_NCBI_accessionNumber | while read line ; do efetch -db protein -id $line -format docsum | xtract -pattern DocumentSummary -element TaxId >> nirS_Tax_IDs; done
# step2 TaxID_to_Tax Info
source activate ete3
cat nirS_Tax_IDs | while read line ; do ete3 ncbiquery --search $line --info >> nirS_tax_info ; done


# combine accession_taxID_taxInfo together
# delete the # information line
sed '/#/d' nirS_tax_info | cut -f 1,4 > nirS_accession_taxID_info
paste nirS_NCBI_accessionNumber nirS_Tax_IDs nirS_accession_taxID_info >output_nirS_NCBI_accession_tax.txt
