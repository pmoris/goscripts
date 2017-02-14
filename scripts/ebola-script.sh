#!/bin/bash
# Main script to reproduce all steps of the analysis.
# Downloads human-ebola PPI from HPIDB2 and Gene Ontology files.
# Performs basic data preprocessing by calling helper scripts (i.e. filtering ebola interactions from HPIDB2 database, finding human genes with existing bat orthologs, etc.).
# Calls GO enrichment script.

echo "Creating directory hpidb2 in data directory to store HPIDB database"
cd ..
mkdir -p data
# File was downloaded in December 2016, HPIDB2.0 version Last Updated June 28, 2016
# PSI-MITAB(2.5) file format.
mkdir -p data/hpidb2
wget -nc http://www.agbase.msstate.edu/hpi/downloads/hpidb2.mitab.zip -O data/hpidb2/hpidb2.zip
cd data/hpidb2
unzip -u hpidb2.zip
cd ../..
rm data/hpidb2/hpidb2.zip

# Extract all interactions involving ebola
echo "Extracting PPI's involving ebola"
head -1 data/hpidb2/*plus* > data/hpidb2/ppi-human-ebola.csv
grep ebola data/hpidb2/*plus* >> data/hpidb2/ppi-human-ebola.csv

# Write all included taxid's to file (note: all human proteins are listed in column 10)
echo "Writing ebola taxids to taxonomyID.csv"
 { tail -n +2 data/hpidb2/ppi-human-ebola.csv | cut -f10  ; tail -n +2 data/hpidb2/ppi-human-ebola.csv | cut -f11 ; } | sort -u > data/hpidb2/taxid-list.txt

# Print number of PPIs for each taxon to screen
# cut -f1 taxonomyID.csv | tail -n +2 | while read x; do echo "$x count: $(cut -f11 ppi-human-ebola.csv | grep $x | wc -l)"; done
cat taxonomyID.csv | tail -n +2 | while read x; do echo "${x} count: $(cut -f10,11 ppi-human-ebola.csv | grep $(echo $x | awk '{print $1}') | wc -l)"; done

# Search Entrez for list of taxonomy ID's using python script, save as taxonomyID.csv
echo "Retrieving taxonomy ID's from Entrez"
cd data/hpidb2
python3 ../../scripts/entrezTaxIDlookup.py taxid-list.txt
cd ../..

# Download 1-1 genome pairs between human and myotis lucifugus from omabrowser.org 
echo "Creating directory oma_orthologs in data directory to store OMA orthologs"
mkdir -p data/oma_orthologs
wget -nc 'http://omabrowser.org/cgi-bin/gateway.pl?f=PairwiseOrthologs&p1=HUMAN&p2=MYOLU&p3=UniProt' -O data/oma_orthologs/orthologs_human_myotis.csv
sed -i '1s/^/uniprot_AC_1\tuniprot_AC_2\torthology_type\tOMA_group\n/' data/oma_orthologs/orthologs_human_myotis.csv

# Check for each human-ebola PPI if a bat ortholog exists for the human protein
echo "Processing PPI data to include information about the existance of bat orthologs, saved to ppi-human-ebola-bat-ortholog.csv in main data directory"
cd data
python3 ../scripts/OMA-orthology-search-pandas.py hpidb2/ppi-human-ebola.csv	oma_orthologs/orthologs_human_myotis.csv 
# outputs:
#	updated hpidb2-human-ebola-ortholog.csv that contains additional column "batOrthologExists"
#	updated hpidb2-human-ebola-ortholog-dedup.csv with same output but omitting duplicate pairs
#	ppi-human-ebola-present-bat.csv that contains only entries with a bat ortholog
#	ppi-human-ebola-missing-bat.csv	that contains only entries without a bat ortholog
cut -f27 ppi-human-ebola-ortholog-dedup.csv > background.txt
cut -f27 ppi-human-ebola-missing-bat.csv > interest.txt
cd ..

# Download gaf and obo files for human GOA - version 2.1 !Generated: 2016-11-28 09:32
echo "Creating directory go_data in data directory to store GO obo and gaf files"
mkdir -p data/go_data
wget -nc ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz -O data/go_data/goa_human.gaf.gz
gzip -d goa_human.gaf.gz
# wget -nc http://purl.obolibrary.org/obo/go/releases/2016-11-26/go.owl data/go.owl
wget -nc http://purl.obolibrary.org/obo/go/releases/2016-11-26/go.obo -O data/go_data/go.obo


# find orthologs
# https://www.biostars.org/p/7568/
# inparanoid, https://www.ncbi.nlm.nih.gov/guide/howto/find-homolog-gene , orthoMCL , https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4779612/

# Requires SOAPpy package
# sudo pip3 install SOAPpy



# NOTE: KEEP RESTON VIRUS OR NOT?
