#!/usr/bin/env python3

# http://biopython.org/DIST/docs/tutorial/Tutorial.html

import sys,csv, time
from Bio import Entrez
Entrez.email = 'pieter.moris@uantwerpen.be'



with open(sys.argv[1],'r') as f:
    taxIDlist = [line.rstrip() for line in f]

with open("taxonomyID.csv",'w') as o:
    writer = csv.writer(o,lineterminator='\n',delimiter='\t')
    writer.writerow(['Taxonomy ID','Scientific Name'])

    for id in taxIDlist:
        time.sleep(0.5)
        handle = Entrez.efetch(db="Taxonomy", id = id)
        records = Entrez.read(handle)
#        o.write(id+'\t'+records[0]["ScientificName"].strip()+'\n')
#        writer.writerow(id + '\t' + records[0]["ScientificName"].strip())
        writer.writerow([id,records[0]['ScientificName'].rstrip()])
