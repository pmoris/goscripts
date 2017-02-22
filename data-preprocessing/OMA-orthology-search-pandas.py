import pandas as pd
import sys, os

# Read in human-ebola PPI network from first positional argumnet
print('Reading PPI entries from',sys.argv[1])
ppi = pd.read_csv(sys.argv[1], delimiter='\t')

# Extract uniprot accession numbers and store in new columns
ppi['uniprot_AC_1'] = ppi.iloc[:, 0].str.rsplit('kb:').str[1]    # human
ppi['uniprot_AC_2'] = ppi.iloc[:, 1].str.rsplit('kb:').str[1]    # ebola

# print(ppi.loc[:,'# protein_xref_1'])
# print(ppi.iloc[:,0].str.rsplit('kb:').str[1])
# print(ppi['# protein_xref_1'].str.extract(r'kb:(.*)',expand=True))
# humanAC = ppi['# protein_xref_1']
# print(humanAC[0][-6:])
# print(humanAC[0].split('kb:')[1])
# import re
# p = re.search(r'kb:(.*)',humanAC[0])
# print(p.group(1))

# Read in human - Myotis lucifugus ortholog pairs, from second positional argument
# uniprot_AC_1 = human       AC_2 = Myotis lucifugus
# orthologs = pd.read_csv('oma_orthologs/orthologs_human_myotis.csv', delimiter='\t')
print('Reading orthologs from', sys.argv[2])
orthologs = pd.read_csv(sys.argv[2], delimiter='\t')

# Add binary indicator to human-bat PPI list if a bat ortholog exists for the human protein
ppi['batOrthologExists'] = ppi['uniprot_AC_1'].isin(orthologs['uniprot_AC_1'])
print('Writing new PPI files in',os.getcwd())
ppi.to_csv('ppi-human-ebola-bat-ortholog.csv', sep='\t', index=False)

# Subset PPIs where human partner (AC1) does not appear in human-bat
# ortholog list
ppi_missingInBat = ppi[~ppi['uniprot_AC_1'].isin(orthologs['uniprot_AC_1'])]
ppi_missingInBat.to_csv('ppi-human-ebola-missing-bat.csv', sep='\t', index=False)

# Subset PPIs where human partner (AC1) does appear in human-bat ortholog list
ppi_presentInBat = ppi[~ppi.index.isin(ppi_missingInBat.index)]
ppi_presentInBat.to_csv('ppi-human-ebola-present-bat.csv', sep='\t', index=False)
# ppi[ppi['uniprot_AC_1'].isin(orthologs['uniprot_AC_1'])]

# Deduplicate AC1/AC2 pairs to separate file for network visualisation
ppi_dedup = ppi.drop_duplicates(subset=['uniprot_AC_1', 'uniprot_AC_2'])
ppi_dedup.to_csv('ppi-human-ebola-bat-ortholog-dedup.csv', sep='\t', index=False)


# Check if subsets add up to original ppi list
# print(ppi_missingInBat['uniprot_AC_1'].unique().size + len(ppi_presentInBat['uniprot_AC_1'].unique()) == ppi['uniprot_AC_1'].unique().size)

# Subset human-bat orthologs where human partner is present in human-ebola
# network
orthologs_in_ppi = orthologs[
    orthologs['uniprot_AC_1'].isin(ppi['uniprot_AC_1'])]
# orthologs_in_ppi2 = orthologs[
#     orthologs['uniprot_AC_1'].isin(ppi_presentInBat['uniprot_AC_1'])]
# orthologs_in_ppi3 = orthologs[orthologs['uniprot_AC_1'].isin(
#     ppi_presentInBat['uniprot_AC_1'])]


# print(ppi_missingInBat['uniprot_AC_1'].unique().size)

# print(ppi_presentInBat['uniprot_AC_1'].unique().size)

# print(ppi['uniprot_AC_1'].unique().size)


# print(ppi)

'''
background for enrichment: positief ortholoog met alle interacties
positief ortholoog met alle mogelijke orthologen
geen ortholoog met alle interacties
'''





















# print(orthologs_in_ppi)

# print(len(orthologs_in_ppi['uniprot_AC_1'].unique()))


# print(orthologs[orthologs['uniprot_AC_1'].isin(ppi_presentInBat['uniprot_AC_1']) &
#                 orthologs['orthology_type'].isin(['1:1', 'n:1'])])

# print(len(orthologs[orthologs['uniprot_AC_1'].isin(ppi_presentInBat['uniprot_AC_1']) &
# orthologs['orthology_type'].isin(['1:1',
# 'n:1'])]['uniprot_AC_1'].unique()))

# # print(ppi)
# print(ppi.iloc[:,-2:])

# print(ppi_missingInBat.index)

# print(ppi_missingInBat.equals(ppi_2))

# print(ppi.loc[ppi_missingInBat.index])
# print(ppi[ppi.index.isin(ppi_missingInBat.index)])
# print(ppi[ppi['uniprot_AC_1'].isin(orthologs['uniprot_AC_1'])].equals(ppi[~ppi.index.isin(ppi_missingInBat.index)]))


# print(ppi_missingInBat)

# import numpy as np

# people_all = pd.DataFrame({ 'ID' : np.arange(5)})
# people_usa = pd.DataFrame({ 'ID' : [3,4,6,7,100]})

# print(people_all)
# print(people_usa)

# print(people_usa[~people_usa['ID'].isin(people_all['ID'])])
