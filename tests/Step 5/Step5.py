def secstructurelib(Library,full16srRNAseq):
    """
    Build a library of secondary structure of 16s rRNA corresponding to a given library of ASD and SD sequence.


    """
    NewLib = {}
    for n in range(0,int(len(Library)/2)):
        full16srRNAseq = full16srRNAseq[0:-12] + Library['ASD{0}'.format(n+1)]
        secstructure = RNAfoldseq(full16srRNAseq)
        NewLib['foldedseq{0}'.format(n+1)] = str(secstructure)
    return NewLib


def ASD_2rystructure_narrow(full16srRNAseq):
    """ Narrow down the library by discarding sequences that forms secondary structure in ASD region
     Import the most narrowed Libdict up to step 5 as well as the secondary structure dictionary.

     :param full16srRNAseq: The rRNA sequence in string form.
     :return: The narrowed dictionary of ASD/ SD that do not have secondary structure in the ASD region.
     """
    WtASD = full16srRNAseq[-12:-1] + full16srRNAseq[-1]
    WtSD = revcomp(WtASD[4:10])
    Library = narrow_binding(WtASD, WtSD)
    print('Library after step 2: ', Library)
    secstructureseqlib = secstructurelib(Library, full16srRNAseq)
    # Locate the index of the first bp of the ASD for all 16S rRNA sequences
    import re
    for pseudoindex in re.finditer('AUCA', full16srRNAseq):
        if full16srRNAseq[pseudoindex.start() + 10:pseudoindex.start() + 12] == 'UA':
            index = int(pseudoindex.start())
    """Iterate through the whole secondary structure dictionary, locating a constant positioned
    ASD and determining whether there is secondary structure in that ASD"""
    for i in range(0, int(len(Library)/2)):
        secname = "foldedseq" + str(i + 1)
        ASDname = "ASD" + str(i + 1)
        SDname = "SD" + str(i + 1)
        secstructureseq = secstructureseqlib[secname]
        ASDrandomregion = secstructureseq[index+4:index + 10]
        if "(" in ASDrandomregion or ")" in ASDrandomregion or "{" in ASDrandomregion or "}" in ASDrandomregion or "[" in ASDrandomregion or "]" in ASDrandomregion:
            del Library[ASDname]
            del Library[SDname]
    Library = updatelibindex(Library)
    print('Library after step 5: ',Library)
    return Library
