
def importGAF(path, geneSet):
    """
    Imports a GAF file and generates a dictionary mapping the gene ID to the GO ID.

    Information on the GAF 2.1 format can be found at 
        http://geneontology.org/page/go-annotation-file-gaf-format-21 

    Parameters
    ----------
    path : str
        The path to the file.
    geneSet : set
        A set containing the ID's of the genes under consideration.

    Returns
    -------
    dict of str
        A dictionary mapping gene ID's to GO ID's.

    Possible improvements:
        Check for `is_obsolete` and `replaced_by`, although the replacement term should be in OBO file as an entry.

        Check for inclusion in provided geneset afterwards using:
            gafDict = {key: value for key, value in gafDict.items() if key in geneSet }

        To do: double dictionary? already filter based on gene list? what to do with NOT?
    """
    gafPath = os.path.abspath(path)

    if not geneSet:
        print('Empty gene set provided during creation of GAF dictionary.')
        exit()

    with open(gafPath, 'r') as gafFile:
        gafDict = {}
        for line in gafFile:
            # ignore comments in gaf file
            if not line.startswith('!'):
                splitLine = line.split('\t')                # split column-wise
                uniprotAC = splitLine[1]
                goTerm = splitLine[4]
                goQualifier = splitLine[3]
                if 'NOT' not in goQualifier:                # ignore annotations with "NOT"
                    if uniprotAC in geneSet:                # only keep genes present in input gene set
                        if uniprotAC not in gafDict:        # Create new key if AC does not already appear in dictionary
                            # make dictionary containing uniprot AC as key and
                            # set of GO's
                            gafDict[uniprotAC] = {goTerm}
                        else:
                            gafDict[uniprotAC].add(goTerm)
    return gafDict

