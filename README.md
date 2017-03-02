# ebola-go
Mining ebola-human PPI GO terms for presence in bats.

Collaboration with Ben Verhees.

## GO enrichment analysis of Ebola virus - human protein-protein interaction network with respect to the existance of orthologous genes in Myotis lucifugus.

[The full analysis is described here.](docs/Overview.ipynb)

The human partners in the protein-protein interaction network (PPI) between Ebola virus and Homo sapiens were investigated to uncover GO divergences between two of its subsets: the proteins lacking an orthologue in the Myotis lucifugus genome and those that do have a corresponding orthologue.

Initial results obtained by comparing these two subsets against each other, i.e. by using the larger set as the background reference set and searching for relative enrichment of GO terms in the smaller set), suggested that no such difference exists.

However, additional comparisons between the two subsets and the entire human genome annotation set revealed a small pocket of GO terms, related to chromatin DNA and nucleosome DNA binding, which were absent from the set of proteins with a corresponding orthologue. Furthermore, this enrichment was still present when the complete protein set was compared against the human annotation set.

Comparing multiple (dependent) GO enrichment analyses in this manner could be plagued by severe statistical issues, not least of which the issue of multiple testing and how to correct for it in this setting.

## GO enrichment tool
[GO enrichment analysis tool](go-enrichment-tool/go_enrichment_script.py). Based on code by Pieter Meysman.

Performs one-sided hypergeometric tests, starting from the most specific (child) GO terms associated with the genes in the set of interest. If the p-value of the test does not fall below the specified significance level alpha, the test wll be carried out for all of the term's parent terms, otherwise the process will terminate. This method attempts to limit the total number of tests that need to be carried out, since a term that is enriched will likely also have enriched parent terms. Furthermore, GO terms associated with a small number of genes are skipped. Next, the Benjamini-Hochberg FDR or Bonferroni multiple testing correction are applied to the test results. Finally, a csv file containing all the GO terms that were evaluated and their p-values are returned.

## Data acquisition

This analysis can be reproduced by executing the provided [bash script](data-preprocessing/data_setup.sh) to retrieve the utilised data and preprocess it, followed by the GO enrichment analyses described in the [documentation](docs/Overview.ipynb).

---

Dependencies:
    `numpy, pandas, scipy.stats, statsmodels.sandbox.stats.multicomp`
