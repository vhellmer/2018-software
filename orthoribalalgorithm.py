# orthoribalalgorithm.py
# Determine Shine Dalgarno (SD) and Anti-Shine Dalgarno sequences for strains of bacteria to design orthogonal ribosomes
# Rice University iGEM 2018
# Vu Hoang Anh and Sai Sriram
# Correspondence: ss162@rice.edu; ahv1@rice.edu
# ---------------------------------------Step 1----------------------------------------------------------


def convertTtoU(sequence):
    """
    Convert thymine to uracil

    :param sequence: The DNA sequence that is being converted to mRNA, in string form
    :return: Another string with T's replaced with U's.
    """
    return sequence.replace('T', 'U')


def revcomp(sequence):
    """
    Find reverse complementary sequence

    :param sequence: The RNA sequence in string form
    :return: The reverse complement sequence in string form
    """
    complement = {"A": "U", "U": "A", "C": "G", "G": "C", "N": "N"}
    revcompseq = ""
    sequence_list = list(sequence)
    sequence_list.reverse()
    for letter in sequence_list:
        revcompseq += complement[letter.upper()]
    return revcompseq


def SDseqnum2str(code):
    """
    Translate a 6-digit number to its responding 6-base SD sequence

    :param code: The 6 digit number code (with numbers 1-4) in string form
    :return: The corresponding RNA sequence
    """
    seq = ''
    codex = 'AUCG'
    for i in range(0, len(code)):
        index = int(code[i])
        seq = seq + codex[index-1]
    return seq


def ASDseqnum2str(code):
    """
    Translate a 6-digit number to its responding 12-base ASD sequence

    :param code: The 6 digit code (numbers 1-4) as a string
    :return: The ASD sequence
    """
    seq = revcomp(SDseqnum2str(code))
    seq = 'AUCA' + seq + 'UA'
    return seq


def libbuild():
    """
    Build a library of 4096 pairs

    Uses 6 for loops to iterate through all possible combinations of 6 digit numbers consisting of 1, 2, 3, and 4.
    Additionally, formats the initial dictionary so that the first SD has the key "SD1", and so on.
    :return: The initial dictionary containing all possible ASD sequences
    """
    Lib = []
    for b1 in range(1,5):
        for b2 in range(1,5):
            for b3 in range(1,5):
                for b4 in range(1,5):
                    for b5 in range(1,5):
                        for b6 in range(1,5):
                            seqcode = str(b1) + str(b2) + str(b3) + str(b4) + str(b5) + str(b6)
                            SDseq = SDseqnum2str(seqcode)
                            ASDseq = ASDseqnum2str(seqcode)
                            Lib.append([SDseq, ASDseq])
    Libdict = {}
    for x in range(0,4096):
        pair = Lib[x]
        Libdict['SD{0}'.format(x+1)] = pair[0]
        Libdict['ASD{0}'.format(x+1)] = pair[1]
    print('4096 pairs are: ', Libdict)
    return Libdict

# ------------------------------------------Step 2---------------------------------------------------------------


def RNAduplexval(ASD,SD):
    """
    Call RNAduplex in python shell; change RNAduplex path used below corresponding to position of RNAduplex.exe

    Pipes results from RNAduplex back into the algorithm
    :param ASD: The ASD sequence as a string
    :param SD: The SD sequence as a string
    :return: A value corresponding to the binding energy of the ASD and SD in kcal/mol
    """
    import subprocess
    input = ASD + '\n' + SD + '\n'
    input = input.encode('utf-8')
    from subprocess import Popen,PIPE,STDOUT
    result = Popen(["C:\\Users\\siddu\\Desktop\\Python Code\\Vienna\\RNAduplex.exe"], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    m = result.communicate(input=input)[0]
    val = m.decode('utf-8')
    n = val.find(':')
    val = val[n+11:-3]
    return val

# -----------------------------------------Step 3----------------------------------------------------------------


def updatelibindex(Library):
    """
    Update the index of the sequence in the library to be continuous

    After narrowing down the dictionary through the various steps of the algorithm, certain terms are deleted. Thus, the
    dictionary is no longer continuous, and could be ordered "ASD1" and the next term "ASD9". This function changes the
    "ASD9" into an "ASD2"
    :param Library: The dictionary whose keys need to be updated.
    :return: The dictionary with updated keys.
    """
    liblen = len(Library) / 2
    n = 1
    while (n < liblen + 1):
        testidx = n
        while 'ASD' + str(testidx) not in Library:
            testidx = testidx + 1
        Library['ASD{0}'.format(n)] = Library.pop('ASD{0}'.format(testidx))
        Library['SD{0}'.format(n)] = Library.pop('SD{0}'.format(testidx))
        n = n + 1
    return Library


def narrow_binding(WtASD,WtSD):
    """
    Input the wildtype ASD and SD, output is the narrowed library with pairs that fall within 0.5 kcal range of the
    wildtype binding energy

    :param WtASD: ASD in string form
    :param WtSD: SD in string form
    :return: The narrowed dictionary with updated key indices only consisting of ASD/ SD pairs that have binding
    energies of -0.5 <= x <= 0.5
    """
    Wildval = RNAduplexval(WtASD,WtSD)
    SequenceLib = libbuild()
    for n in range(0,4096):
        ASD = SequenceLib['ASD'+str(n+1)]
        SD = SequenceLib['SD' + str(n+1)]
        val = RNAduplexval(ASD,SD)
        if float(Wildval) - 0.5 >= float(val) or float(val) >= float(Wildval) + 0.5:
            del SequenceLib['ASD'+str(n+1)]
            del SequenceLib['SD'+str(n+1)]
    Libdict = updatelibindex(SequenceLib)
    return Libdict

# ---------------------------Step 4-------------------------------------------------------------------


def RNAfoldseq(fullRNAsequence):
    """
    Call RNAfold program in python, ***change full path of program accordingly

    :param fullRNAsequence: The RNA sequence in string form.
    :return: The free energy of the folded structure.
    """
    import subprocess
    input = fullRNAsequence
    input = input.encode('utf-8')
    from subprocess import Popen,PIPE,STDOUT
    result = Popen(["C:\\Users\\siddu\\Desktop\\Python Code\\Vienna\\RNAfold.exe"], stdout=subprocess.PIPE,stdin=subprocess.PIPE)
    m = result.communicate(input=input)[0]
    val = m.decode('utf-8')
    val = val[len(fullRNAsequence)+1:len(fullRNAsequence)*2+1]
    return val
#-------------------------------Step 5--------------------------


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

# --------------------------Step 6-7-----------------------------------


def narrow_crossbinding(rRNA):
    """
    Narrow down the library more by eliminating pairs with ASD that has < 1 kcal/mol binding energy with
    wild type SD

    :param rRNA: The rRNA sequence in string form.
    :return: Further narrowed dictionary.
    """
    Library = ASD_2rystructure_narrow(rRNA)
    WtSD = revcomp(rRNA[1534:1540])
    for n in range(0,int(len(Library)/2)):
        ASD = 'ASD' + str(n+1)
        SD = 'SD' + str(n+1)
        if float(RNAduplexval(Library[ASD],WtSD)) < -1:
            del Library[ASD]
            del Library[SD]
    Library = updatelibindex(Library)
    return Library


def orthoribalgorithm():
    """
    Read the file with the genome and call other functions to narrow library, basically a calling function.

    :return: The library after step 7.
    """
    with open(r'C:\Users\siddu\Desktop\Python Code\ecoligenome.txt', 'r') as myfile:
        data = myfile.read().replace('\n', '')
    rRNA = data[4166658:4168200]
    rRNA = convertTtoU(rRNA)
    Libraryafterstep7 = narrow_crossbinding(rRNA)
    print('Library after step 7: ', Libraryafterstep7)
    return Libraryafterstep7

# ----------------------- Step 8 -----------------------------------------------------


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

# -----------------------------------Step 9---------------------------------------------


def ASDTIRbinding(ASD, TIR):
    """
    Call RNAduplex in python shell to determine binding energies of ASD and TIRs.

    :param ASD: ASD sequence as a string.
    :param TIR: TIR sequence as a string.
    :return: The binding energy of the two in kcal/mol calculated by RNAduplex.
    """
    import subprocess
    input = ASD + '\n' + TIR + '\n'
    input = input.encode('utf-8')
    from subprocess import Popen, PIPE, STDOUT
    result = Popen(["C:\\Users\\siddu\\Desktop\\Python Code\\Vienna\\RNAduplex.exe"], stdout=subprocess.PIPE,
                   stdin=subprocess.PIPE)
    m = result.communicate(input=input)[0]
    val = m.decode('utf-8')
    n = val.find(':')
    val = val[n + 11:-3]
    return val

# ---------------------------------Step 10-----------------------------------------------------


def allASDTIRpairs():
    """
    Iterate through all possible ASD TIR pairs and find the ones with the highest average binding energies with host
    TIRs

    We wish to find the ASDs that do not bind well with the host translation initiation regions to assure orthogonality
    of the ribosomes, so here we choose the ASDs that have the highest binding energies (i.e. don't bind well with host
    TIRs)
    :return: Prints the list of the top ten ASD candidates.
    """
    currentlib = orthoribalgorithm()
    TIRdict = getallTIRs()
    dictofvals = {}
    print("Number of TIRs: " + str(len(TIRdict)))
    listofaverages = []
    for i in range(0, round(len(currentlib) / 2)): # iterate through all ASDs
        listofvals = []
        ASDname = str('ASD' + str(i + 1))
        for j in range(0, len(TIRdict)): # for each ASD, iterate through all TIRs in the genome
            TIRname = str('TIR' + str(j + 1))
            val = float(ASDTIRbinding(currentlib[ASDname], TIRdict[TIRname]))
            listofvals.append(val)
        average = sum(listofvals) / len(listofvals)
        dictofvals[average] = ASDname  # calculate the average binding energy between the ASD
        # and all TIRs; here we store the key as the average so that the we can call the names of the highest ASDs after
        # the list is sorted
        listofaverages.append(average)

    listofaverages.sort(reverse=True)
    print('Here are the 10 top candidates with highest ASD-host binding values:')
    for i in range(0, 10):
        print(dictofvals[listofaverages[i]])

allASDTIRpairs()