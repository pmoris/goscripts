## GO enrichment analysis tools

Pieter Moris 2017
ADReM - biomina - UAntwerpen

Performs one-sided hypergeometric tests, starting from the most specific (child) GO terms associated with the genes in the set of interest. If the p-value of the test does not fall below the specified significance level alpha, the test wll be carried out for all of the term's parent terms, otherwise the process will terminate. This method attempts to limit the total number of tests that need to be carried out, since a term that is enriched will likely also have enriched parent terms. Furthermore, GO terms associated with a small number of genes are skipped. Next, the Benjamini-Hochberg FDR or Bonferroni multiple testing correction are applied to the test results. Finally, a csv file containing all the GO terms that were evaluated and their p-values are returned. More information is available in the docstrings.

---

Dependencies:
    `numpy, pandas, scipy.stats, statsmodels.sandbox.stats.multicomp`
