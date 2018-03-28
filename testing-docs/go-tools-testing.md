# Testing go-tools output after re-writing functions

# TODO:

- Namespace filtering should happen during initial import, otherwise there could be jumps across namespaces (if capable_of and part_of relations are included).
- download files on request

## Commands used

### Converting copy-paste of GO terms from QuickGo to python list.

`pieter@pieter-XPS-15-9560 ~/Desktop $ grep GO temp.txt | tr '\n' ',' | sed 's/,/","/g' | sed 's/.\{2\}$//g'`

### Arguments used for go enrichment test:

```
pieter@pieter-XPS-15-9560 /media/pieter/DATA/Wetenschap/Doctoraat/host-pathogen-project/go-tools/go-tools $ python go_enrichment_script.py -b ../../ebola-project/ebola-network-GO/data/background.txt -s ../../ebola-project/ebola-network-GO/data/interest.txt -o ../../ebola-project/ebola-network-GO/data/go_data/go.obo -g ../../ebola-project/ebola-go/data/go_data/gene_association.gaf -O output-old.txt -n biological_process -m 3 -T 0.05 -t 0.1
```

```
pieter@pieter-XPS-15-9560 /media/pieter/DATA/Wetenschap/Doctoraat/host-pathogen-project/go-tools/go-tools $ python go_enrichment_script.py -b /media/pieter/DATA/Wetenschap/Doctoraat/host-pathogen-project/ebola-project/ebola-network-GO/data/background.txt -s /media/pieter/DATA/Wetenschap/Doctoraat/host-pathogen-project/ebola-project/ebola-network-GO/data/interest.txt -o /media/pieter/DATA/Wetenschap/Doctoraat/host-pathogen-project/ebola-project/ebola-network-GO/data/go_data/go.obo -g /media/pieter/DATA/Wetenschap/Doctoraat/host-pathogen-project/ebola-project/ebola-go/data/go_data/gene_association.gaf -O /media/pieter/DATA/Wetenschap/Doctoraat/host-pathogen-project/go-tools/go-tools/output-1.txt -n biological_process -m 3 -l 0.05 -t 0.1
```

# Differences between Cytoscape BINGO and go-tools
## go-tools

47,GO:0034641,cellular nitrogen compound metabolic process,biological_process,0.2229713478320156,0.97467920494363,12/25 (48.00%),55/140 (39.29%)

## BINGO

GO-ID	p-value	corr p-value	x	n	X	N	Description	Genes in test set
34641	8.0298E-2	8.8290E-1	14	58 25	140		cellular nitrogen compound metabolic process	Q8TAQ2|A8CG34|P35637|P31943|Q32P51|P0DMV9|P49411|Q71UM5|Q5VTE0|P35268|Q14839|P51532|P27708|Q16695

## Difference?

12/25 vs 14/25
55/140 vs 58/140

## Debugging

```
viable = GOterms['GO:0034641'].recursive_children
viable.add('GO:0034641')
len(viable)
2096

enrichment_stats.countGOassociations(viable, gafDict)
55
enrichment_stats.countGOassociations(viable, gafSubset)
12
len(gafSubset)
25
len(gafDict)
140

for gene, GOids in gafSubset.items():
    if not viable.isdisjoint(GOids):
        print(gene)

A8CG34
P0DMV9
P27708
P31943
P35268
P35637
P51532
Q14839
Q16695
Q32P51
Q71UM5
Q8TAQ2
```

**Additional proteins present in BINGO**:
- P49411
- Q5VTE0

Find additional GO terms:

    protein in ['P49411', 'Q5VTE0']:
      for j in gafDict[protein]:
        if j not in viable:
          print(j)

      GO:0045471
      GO:0006414
      GO:0008150
      GO:0006414

```
gafDict['P49411']
{'GO:0006414', 'GO:0045471'}

gafDict['Q5VTE0']
{'GO:0006414', 'GO:0008150'}

'GO:0006414' in viable
False
'GO:0045471' in viable
False
'GO:0045471' in GOterms
True
'GO:0006414' in GOterms
True
```

    viable.update(['GO:0006414', 'GO:0045471'])
    enrichment_stats.countGOassociations(viable, gafSubset)
    14
    enrichment_stats.countGOassociations(viable, gafDict)
    58

**Reason**: 'GO:0006414' is only a child of GO:0034641 if **'part_of'** relations are taken into account.

***Fixed by changing ignore_partof to False!***

# More differences between GO tools and BINGO
    48731	3.7889E-1	8.8290E-1	6	28  25	140		system development	Q8TAQ2|P27708|P35268|P35613|P51532|P12036

    2,GO:0048731,system development,biological_process,0.08221993833504583,0.97467920494363,2/25 (8.00%),3/140 (2.14%)

    P51532
    Q8TAQ2

    for protein in ['P27708', 'P35268', 'P35613', 'P12036']:
      for j in gafDict[protein]:
          if j not in viable:
              print(j)

    GO:0017144
    GO:0031100
    GO:0042594
    GO:0007595
    GO:0046134
    GO:0031000
    GO:0051414
    GO:0006207
    GO:0071364
    GO:0006541
    GO:0035690
    GO:0032868
    GO:0044205
    GO:0006228
    GO:0007507
    GO:0007565
    GO:0046777
    GO:0019240
    GO:0014075
    GO:0018107
    GO:0001889
    GO:0033574
    GO:0006614
    GO:0019083
    GO:0006364
    GO:0006412
    GO:0002181
    GO:0046632
    GO:0006413
    GO:0000184
    GO:0051591
    GO:0030198
    GO:0042475
    GO:0022617
    GO:0050900
    GO:0007566
    GO:0007166
    GO:0046697
    GO:0015718
    GO:0043434
    GO:0072661
    GO:0046689
    GO:0006090
    GO:1902513
    GO:0033693
    GO:0061564
    GO:0048936
    GO:0000226
    GO:0007409
    GO:0030031

GO:0007595, GO:0007507 are children of system development according to *part_of* relations. Further terms were not checked, but likely same reason.


# Even more
    429,GO:0010467,gene expression,biological_process,0.8018958974197338,0.97467920494363,1/25 (4.00%),8/140 (5.71%)

P31943

    10467	9.0520E-2	8.8290E-1	11	43 25	140		gene expression	Q8TAQ2|A8CG34|P35637|P31943|Q32P51|P49411|Q71UM5|Q5VTE0|P35268|Q14839|P51532

# Additional top hits in BINGO not present in either method
***Fixed by changing the minimum number of genes to 2!***

    14074	3.0832E-2	8.8290E-1	2	22  5	140		response to purine-containing compound	P27708|P35613

    51091	3.0832E-2	8.8290E-1	2	22  5	140		positive regulation of sequence-specific DNA binding transcription factor activity	P0DMV9|P51532

**Not reached in either new or old!**

## 14074 first

```
viable = GOterms['GO:0014074'].recursive_children
viable.add('GO:0014074')
len(viable)
19
enrichment_stats.countGOassociations(viable, gafDict)
2
enrichment_stats.countGOassociations(viable, gafSubset)
2
viable
{'GO:1901560', 'GO:0070305', 'GO:1905795', 'GO:0072764', 'GO:0071320', 'GO:1905794', 'GO:1904310', 'GO:0014074', 'GO:0033198', 'GO:0071313', 'GO:0031319', 'GO:0071321', 'GO:0031000', 'GO:0071415', 'GO:0051591', 'GO:0071318', 'GO:0072754', 'GO:1904309', 'GO:1901596'}
for gene, GOids in gafSubset.items():
    if not viable.isdisjoint(GOids):
        print(gene)
P27708
P35613

for gene, GOids in gafDict.items():
    if not viable.isdisjoint(GOids):
        print(gene)
P27708
P35613
```

**Exactly the same proteins are found for subset, but also for background!**

Can't find out which additional terms were found for background in BINGO...

## Comparison between new and old
###new:
    GOterms['GO:0014074'].children
    {'GO:1901596', 'GO:0071415', 'GO:0031000', 'GO:0070305', 'GO:0033198', 'GO:1904309', 'GO:1905794', 'GO:0051591', 'GO:1901560'}
    GOterms['GO:0014074'].parents
    {'GO:0010243', 'GO:0014070'}
    GOterms['GO:0014074'].recursive_parents
    {'GO:0010033', 'GO:0010243', 'GO:0042221', 'GO:0008150', 'GO:0009719', 'GO:1901698', 'GO:0014070', 'GO:0050896'}
    GOterms['GO:0014074'].recursive_children
    {'GO:0071313', 'GO:1901596', 'GO:0071415', 'GO:0031319', 'GO:0071318', 'GO:0031000', 'GO:0070305', 'GO:0071321', 'GO:0033198', 'GO:0072754', 'GO:0072764', 'GO:1904309', 'GO:1904310', 'GO:1905794', 'GO:0051591', 'GO:1901560', 'GO:1905795', 'GO:0071320'}
    len(GOterms['GO:0014074'].recursive_children)
    18
    len(GOterms['GO:0014074'].recursive_parents)
    8
    len(GOterms['GO:0014074'].parents)
    2
    len(GOterms['GO:0014074'].children)
    9

### old:
    len(GOterms['GO:0014074'].parents)
    8
    GOterms['GO:0014074'].parents
    {'GO:0009719', 'GO:0008150', 'GO:0010033', 'GO:0050896', 'GO:0042221', 'GO:0014070', 'GO:1901698', 'GO:0010243'}

## Ancestors

### 0010243
285,GO:0010243,response to organonitrogen compound,biological_process,0.6201797372031462,0.97467920494363,2/25 (8.00%),11/140 (7.86%)

10243	6.2018E-1	9.0240E-1	2	11 25	140		response to organonitrogen compound	P27708|P35613

No difference in association counts.

### 0014070
61,GO:0014070,response to organic cyclic compound,biological_process,0.29110760339774083,0.97467920494363,2/25 (8.00%),6/140 (4.29%)

14070	6.2018E-1	9.0240E-1	2	11 25	140		response to organic cyclic compound	P27708|P35613

Slightly more proteins associated with term according to BINGO.

## None of its children are being tested.

# Mismatches between old and new functionality
***THIS WAS FIXED IN COMMIT a67ac4e0522356ab293ba58f70bfb5b787287cc2 BY REWRITING DEF RECURSIVE_TESTER***

old:
`6,GO:0043044,ATP-dependent chromatin remodeling,biological_process,0.1079838499700676,0.97467920494363,3/25 (12.00%),7/140 (5.00%)`

new:
`0,GO:0043044,ATP-dependent chromatin remodeling,biological_process,0.039773983426274626,1.0,3/25 (12.00%),5/140 (3.57%)`

BINGO:
`43044	1.0798E-1	8.8290E-1	3/25  12.0%	7/140  5.0%	Q8TAQ2 Q14839 P51532`

new:
```
GOterms['GO:0043044'].parents
{'GO:0006338'}
GOterms['GO:0043044'].children
{'GO:0043486'}
GOterms['GO:0043044'].recursive_children
{'GO:0035093', 'GO:0035042', 'GO:0043486', 'GO:0034080'}
GOterms['GO:0043044'].recursive_parents
{'GO:0006338', 'GO:0008150', 'GO:0009987', 'GO:0006325', 'GO:0071840', 'GO:0016043'}
len(GOterms['GO:0043044'].recursive_parents)
6
```

old:
```
GOterms['GO:0043044'].parents
{'GO:0071840', 'GO:0008150', 'GO:0016043', 'GO:0009987', 'GO:0006338', 'GO:0006325'}
6
GOterms['GO:0043044'].children
{'GO:0043486', 'GO:0035042', 'GO:0035093', 'GO:0034080'
4
```

New and old method still result in same GO annotation tree...









mismatch

old:
285,GO:0010243,response to organonitrogen compound,biological_process,0.6201797372031462,0.97467920494363,2/25 (8.00%),11/140 (7.86%)


new:
49,GO:0010243,response to organonitrogen compound,biological_process,0.44835560123329843,1.0,1/25 (4.00%),3/140 (2.14%)

BINGO:
false	10243	6.2018E-1	9.0240E-1	2/25  8.0%	11/140  7.8%	P27708 P35613

new:
len(GOterms['GO:0010243'].recursive_parents)
6
len(GOterms['GO:0010243'].recursive_children)
264
len(GOterms['GO:0010243'].children)
74
len(GOterms['GO:0010243'].parents)
3
GOterms['GO:0010243'].parents
{'GO:1901698', 'GO:0009719', 'GO:0010033'}
GOterms['GO:0010243'].recursive_parents
{'GO:1901698', 'GO:0042221', 'GO:0009719', 'GO:0050896', 'GO:0008150', 'GO:0010033'}

old:
GOterms['GO:0010243'].parents
{'GO:0009719', 'GO:0008150', 'GO:0010033', 'GO:0050896', 'GO:0042221', 'GO:1901698'}
len(GOterms['GO:0010243'].parents)
6
len(GOterms['GO:0010243'].children)
264

## Old mismatch example

GOterms['GO:0008150'].children
{'GO:0032502', 'GO:0001906', 'GO:0032501', 'GO:0051179', 'GO:0098754', 'GO:0050896', 'GO:0023052', 'GO:0000003', 'GO:0065007', 'GO:0040007', 'GO:0044699', 'GO:0071840', 'GO:0044848', 'GO:0051704', 'GO:0022414', 'GO:0022610', 'GO:0007610', 'GO:0099531', 'GO:0040011', 'GO:0008152', 'GO:0098743', 'GO:0002376', 'GO:0048511', 'GO:0009987'}
l = GOterms['GO:0008150'].children
len(l)
24
biological_process_children = ["GO:0008150","GO:0099531","GO:0002376","GO:0009758","GO:0000003","GO:0023052",
                                  "GO:0051179","GO:0098754","GO:0050896","GO:0022414","GO:0065007","GO:0001906",
                                  "GO:0040007","GO:0050789","GO:0071840","GO:0032502","GO:0032501","GO:0044848",
                                  "GO:0022610","GO:0048511","GO:0048518","GO:0008152","GO:0098743","GO:0006791",
                                  "GO:0019740","GO:0006794","GO:0051704","GO:0007610","GO:0015976","GO:0043473",
                                  "GO:0009987","GO:0040011","GO:0048519","GO:0008283"]
len(biological_process_children)
34












<!-- # Pytest

def test_children():
    GOterms = obo_tools.importOBO('/media/pieter/DATA/Wetenschap/Doctoraat/host-pathogen-project/ebola-project/ebola-network-GO/data/go_data/go.obo')
    obo_tools.buildGOtree(GOterms, ['GO:0008150', 'GO:0005575', 'GO:0003674'])
    biological_process_children = ["GO:0099531","GO:0002376","GO:0009758","GO:0000003","GO:0023052",
                                  "GO:0051179","GO:0098754","GO:0050896","GO:0022414","GO:0065007","GO:0001906",
                                  "GO:0040007","GO:0050789","GO:0071840","GO:0032502","GO:0032501","GO:0044848",
                                  "GO:0022610","GO:0048511","GO:0048518","GO:0008152","GO:0098743","GO:0006791",
                                  "GO:0019740","GO:0006794","GO:0051704","GO:0007610","GO:0015976","GO:0043473",
                                  "GO:0009987","GO:0040011","GO:0048519","GO:0008283"]
    assert GOterms['GO:0008150'].children == set(biological_process_children)



# def test_recursive_children():
#
#
#
# def test_parents():
#
#
# def test_recursive_parents():
#
#
# def test_namespace_uniformity():
# -->
