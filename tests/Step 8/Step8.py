def getallTIRs():
    """
    *****************************************************************
    Builds a dictionary of strings that contains all the TIRS, accessible via
    the TIR number.

    eg TIRdict['TIR1'] returns the first TIR in the genome. Total of
    3746 TIRs in the E. coli genome.
    :param file: The file used here is a .txt with the sequence 50 bps upstream of all transcription start sites.
    :return: A dictionary containing all translation initiation regions, starting 20 bps in front of the start codon and
    including the start codon itself. Formatted to access the first TIR of the genome by key "TIR1"
    """

    """NOTE: WE WILL HAVE TO REDESIGN HOW THE TIRS ARE OBTAINED FOR USE WITH OTHER CELL LINES """
    with open('C:/users/siddu/Desktop/Python Code/fiftyupstreamecoli.txt', 'r') as genefile:
        i = 1
        TIRdict = {}
        for line in genefile:
            line = line.strip()
            seq = line[-20:]
            TIR = str(convertTtoU(seq))+'AUG'
            TIRname = 'TIR' + str(i)
            i = i + 1
            TIRdict[TIRname] = TIR
        return TIRdict
