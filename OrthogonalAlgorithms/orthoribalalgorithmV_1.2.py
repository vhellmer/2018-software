# orthoribalalgorithm.py
# Determine Shine Dalgarno (SD) and Anti-Shine Dalgarno sequences for strains of bacteria to design orthogonal ribosomes
# Rice University iGEM 2018
# Vu Hoang Anh and Sai Sriram
# Correspondence: ss162@rice.edu; ahv1@rice.edu


# ---------------------------------------Step 1----------------------------------------------------------
# -------------------------------- Initialize library --------------------------------------------------


def libbuild():
    """
    Read the initial library of 4096 pairs from a InitialLibrary file, including
    all possible combinations from randomizing 6 bases in the SD/ASD region.
    Additionally, formats the initial dictionary so that the first SD has the key "SD1", and so on.
    :return: The initial dictionary containing all possible SD/ASD sequences
    """
    import pickle
    pickle_in = open("InitialLibrary","rb")
    IniLibrary = pickle.load(pickle_in)
    #print('4096 pairs are: ', Libdict)
    return IniLibrary

# --------------------------------- Supporting functions -------------------------------------------------------

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


def RNAduplexval(firstseq, secondseq):
    """
    Call RNAduplex in python shell; change RNAduplex path used below corresponding to position of RNAduplex.exe

    Pipes results from RNAduplex back into the algorithm
    :param firstseq: The ASD sequence as a string
    :param secondseq: The SD sequence as a string
    :return: A value corresponding to the binding energy of the ASD and SD in kcal/mol
    """
    import subprocess
    import os
    input = firstseq + '\n' + secondseq + '\n'
    input = input.encode('utf-8')
    from subprocess import Popen,PIPE,STDOUT
    RNAduplexpath = os.path.join(os.path.dirname(__file__), "Vienna\\RNAduplex.exe")
    result = Popen([RNAduplexpath], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    m = result.communicate(input=input)[0]
    val = m.decode('utf-8')
    n = val.find(':')
    val = val[n+11:-3]
    return val


def updatelibindex(library):
    """
    Update the index of the sequence in the library to be continuous

    After narrowing down the dictionary through the various steps of the algorithm, certain terms are deleted. Thus, the
    dictionary is no longer continuous, and could be ordered "ASD1" and the next term "ASD9". This function changes the
    "ASD9" into an "ASD2"
    :param library: The dictionary whose keys need to be updated.
    :return: The dictionary with updated keys.
    """
    liblen = len(library) / 2
    n = 1
    while (n < liblen + 1):
        testidx = n
        while 'ASD' + str(testidx) not in library:
            testidx = testidx + 1
        library['ASD{0}'.format(n)] = library.pop('ASD{0}'.format(testidx))
        library['SD{0}'.format(n)] = library.pop('SD{0}'.format(testidx))
        n = n + 1
    return library


def RNAfoldseq(fullRNAsequence):
    """
    Call RNAfold program in python using the default option utilizing MFE structure.

    :param fullRNAsequence: The RNA sequence in string form.
    :return: The free energy of the folded structure.
    """
    import subprocess
    import os
    input = fullRNAsequence
    input = input.encode('utf-8')
    from subprocess import Popen,PIPE,STDOUT
    RNAfoldpath = os.path.join(os.path.dirname(__file__), "Vienna\\RNAfold.exe")
    result = Popen([RNAfoldpath], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    m = result.communicate(input=input)[0]
    val = m.decode('utf-8')
    val = val[len(fullRNAsequence)+1:len(fullRNAsequence)*2+1]
    return val



def RNAfoldcentroidseq(fullRNAsequence):
    """
    Call RNAfold program in python using the centroid option.

    :param fullRNAsequence: The RNA sequence in string form.
    :return: The free energy of the folded structure.
    """
    import subprocess
    import os
    input = fullRNAsequence
    input = input.encode('utf-8')
    from subprocess import Popen,PIPE, STDOUT
    RNAfoldpath = os.path.join(os.path.dirname(__file__), "Vienna\\RNAfold.exe")
    result = Popen([RNAfoldpath, "-p"], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    m = result.communicate(input=input)[0]
    val = m.decode('utf-8')
    index = val.find('d')
    val = val[index-10-len(fullRNAsequence):index-10]
    return val


def getallTIRs():
    """
    *****************************************************************
    Builds a dictionary of strings that contains all the TIRS, accessible via
    the TIR number.

    eg TIRdict['TIR1'] returns the first TIR in the genome. TIRs are determined from the CDS file of the species as well
    as a .txt of the genome file, containing only base pairs.
    :param files: The files used here are the species CDS and genome.
    :return: A dictionary containing all translation initiation regions, starting 20 bps in front of the start codon and
    including the start codon itself. Formatted to access the first TIR of the genome by key "TIR1"
    """

    with open('C:/users/siddu/Desktop/Python Code/E coli K-12 MG1655 CDSs.txt', 'r') as cdsfile:  # get all CDSs
        # Note that the CDS file sometimes skips some numbers -- for example goes from 54 to 56.
        cdslist1 = []
        cdslist2 = []
        for line in cdsfile:
            idx1 = line.find('[location=')
            if idx1 != -1:
                idx2 = line.find('..')
                if line[idx1 + 10: idx1 + 13] == 'com':
                    # complement cases -- take the second number then reverse complement
                    if line[idx1 + 21] == 'j':  # join cases -- unaccounted for
                        pass
                    else:
                        idx2_1 = line.find(')]\n')

                        cdslist2.append(int(line[idx2 + 2: idx2_1]) - 1)  # -1 to account for Python indexing
                        pass
                else:
                    cdslist1.append(int(line[idx1 + 10: idx2]) - 1)
            else:
                pass
        print('Number of CDS = ' + str(len(cdslist1)))
        print('Number of complement CDS = ' + str(len(cdslist2)))
    # MUST ACCOUNT FOR THE JOIN CASES
    # ------------------------ Regular cases when TIRs are found on the genome strand ---------------------------------
    with open('C:/users/siddu/Desktop/Python Code/ecoligenome.txt', 'r') as genomefile:
        index = 0
        i = 0
        previousline = ''
        TIRdict = {}
        for line in genomefile:  # cdslist1 deals with the forward sequences
            line = line.rstrip()
            currentandprevious = previousline.rstrip() + line.rstrip()
            # we concatenate the current and previous lines to account for
            # CDS that overlap over two lines
            if i < len(cdslist1):
                idx3 = cdslist1[i] - index

                if idx3 > 66:
                    # this is the case when the end 3 bps overlap into the next line or we haven't yet reached
                    # the line of the CDS yet
                    previousline = line.rstrip()
                    index += 70  # each line is 70 bps long
                    pass

                else:
                    # here we have the case when we already passed the start of this TIR, so we have to go back to the
                    # previous line
                    a = convertTtoU(currentandprevious[idx3 + 49:idx3 + 73])

                    TIRdict['TIR' + str(i + 1)] = a  # 59 to 73
                    i = i + 1
                    previousline = line.rstrip()
                    index += 70  # each line is 70 bps long
                    pass

    # ------------------------------- Complement cases where TIRs are found on complementary strand ----------------------
    with open('C:/users/siddu/Desktop/Python Code/ecoligenome.txt', 'r') as genomefile:
        index = 0
        j = i
        i = 0
        previousline = ''
        for line in genomefile:  # This loop handles cdslist2, which contains the reverse complement sequences.
            line = line.rstrip()
            currentandprevious = previousline.rstrip() + line.rstrip()
            # we concatenate the current and previous lines to account for
            # CDS that overlap over two lines
            if i < len(cdslist2):
                idx3 = cdslist2[i] - index
                if idx3 > 49:
                    # this is the case when the end 3 bps overlap into the next line or we haven't yet reached
                    # the line of the CDS yet
                    previousline = line.rstrip()
                    index += 70  # each line is 70 bps long
                    pass
                else: 
                    # here we have the case when we already passed the start of this TIR, so we have to go back to the
                    # previous line
                    a = convertTtoU(revcompDNA(currentandprevious[idx3 + 68:idx3 + 92]))
                    TIRdict['TIR' + str(i + 1 + j)] = a
                    i = i + 1
                    previousline = line.rstrip()
                    index += 70  # each line is 70 bps long
                    pass
    return TIRdict


def secstructurelib(Library):
    """
    Build a library of secondary structure of 16s rRNA corresponding to a given library of ASD and SD sequence.
    :param Library: The library of ASD SD sequence
    :return: A library containing secondary structure of rRNA with each ASD
    """
    full16srRNAseq = get16srRNAseq()
    NewLib = {}
    for n in range(0,int(len(Library)/2)):
        full16srRNAseq = full16srRNAseq[0:-12] + Library['ASD{0}'.format(n+1)]
        secstructure = RNAfoldcentroidseq(full16srRNAseq)
        NewLib['foldedseq{0}'.format(n+1)] = str(secstructure)
    return NewLib


def get16srRNAseq():
    """
    Get the full 16s rRNA sequence from the genome
    :return: the full 16s rRNA sequence from the genome
    """
    import os
    filepath = os.path.join(os.path.dirname(__file__), "ecoligenome.txt")   
    with open(filepath,'r') as myfile:  # determine rRNA sequence
        data = myfile.read().replace('\n', '')
    full16srRNAseq = convertTtoU(data[4166658:4168200])
    return full16srRNAseq

# ---------------------- Define a class to model the library and apply each step to that object -----------------


class Library:
    """A model for the library of ASD/ SD sequences that will be narrowed as more steps are applied"""

    def __init__(self, library):
        """

        :param genome: The string corresponding to the file path with the strain's genome file
        :param library: The library itself
        """
      # ******************* Still have to make self.genome a parameter by moving get16srRNAseq into the class
        self.library = library

    # -----------------------------------------Steps 2-3----------------------------------------------------------------

    def narrow_binding(self):
        """
        Input the wildtype ASD and SD, output is the narrowed library with pairs that fall within 0.5 kcal range of the
        wildtype binding energy

        :param WtASD: Wild-type ASD in string form
        :param WtSD: Wild-type SD in string form
        :return: The narrowed dictionary with updated key indices only consisting of ASD/ SD pairs that have binding
        energies of -0.5 <= x <= 0.5
        """
        full16srRNAseq = get16srRNAseq()

        WtASD = full16srRNAseq[-12:]
        WtSD = revcomp(WtASD[4:10])
        Wildval = RNAduplexval(WtASD, WtSD)  # Calculate the binding energies of Wild-type ASD/SD pair
        # self.library is the initial library after step 1
        for n in range(0, 4096):
            ASD = self.library['ASD'+str(n+1)]
            SD = self.library['SD' + str(n+1)]
            val = RNAduplexval(ASD,SD)
            if float(Wildval) - 0.5 >= float(val) or float(val) >= float(Wildval) + 0.5:  # Compare each ASD/SD binding
                # energy to wild-type binding energy
                del self.library['ASD'+str(n+1)]
                del self.library['SD'+str(n+1)]
        self.library = updatelibindex(self.library)
        print('Library after step 3: ', self.library)
        return self.library

    # -------------------------------Steps 4-5--------------------------
    
    
    def narrow_crossbinding(self):
        """
        Narrow down the library more by eliminating pairs with ASD that has < 1 kcal/mol binding energy with
        wild type SD

        :param rRNA: The rRNA sequence in string form.
        :return: Further narrowed dictionary.
        """
        rRNA = get16srRNAseq()
        WtSD = revcomp(rRNA[1534:1540])
        for n in range(0,int(len(self.library)/2)):
            ASD = 'ASD' + str(n+1)
            SD = 'SD' + str(n+1)
            if float(RNAduplexval(self.library[ASD],WtSD)) < -1:
                del self.library[ASD]
                del self.library[SD]
        self.library = updatelibindex(self.library)
        print('Library after step 5:', self.library)
        return self.library

    # --------------------------Steps 6-7-----------------------------------
    
    
    def ASD_2rystructure_narrow(self):
        """ Narrow down the library by discarding sequences that forms secondary structure in ASD region
         Import the most narrowed Libdict up to step 5 as well as the secondary structure dictionary.

         :param full16srRNAseq: The rRNA sequence in string form.
         :return: The narrowed dictionary of ASD/ SD that do not have secondary structure in the ASD region.
         """

        full16srRNAseq = get16srRNAseq()
        secstructureseqlib = secstructurelib(self.library)
        # Locate the index of the first bp of the ASD for all 16S rRNA sequences
        import re
        for pseudoindex in re.finditer('AUCA', full16srRNAseq):
            if full16srRNAseq[pseudoindex.start() + 10:pseudoindex.start() + 12] == 'UA':
                index = int(pseudoindex.start())
        """Iterate through the whole secondary structure dictionary, locating a constant positioned
        ASD and determining whether there is secondary structure in that ASD"""
        for i in range(0, int(len(self.library)/2)):
            secname = "foldedseq" + str(i + 1)
            ASDname = "ASD" + str(i + 1)
            SDname = "SD" + str(i + 1)
            secstructureseq = secstructureseqlib[secname]
            ASDrandomregion = secstructureseq[index+4:index + 10]
            if "(" in ASDrandomregion or ")" in ASDrandomregion or "{" in ASDrandomregion or "}" in ASDrandomregion or "[" in ASDrandomregion or "]" in ASDrandomregion:
                del self.library[ASDname]
                del self.library[SDname]
        self.library = updatelibindex(self.library)
        print('Library after step 7: ', self.library)
        return self.library


# ---------------------------------Steps 8-10-----------------------------------------------------

    def allASDTIRpairs(self):
        """
        Iterate through all possible ASD TIR pairs and find the ones with the highest average binding energies with host
        TIRs

        We wish to find the ASDs that do not bind well with the host translation initiation regions to assure orthogonality
        of the ribosomes, so here we choose the ASDs that have the highest binding energies (i.e. don't bind well with host
        TIRs)
        :return: Prints the list of the top ten ASD candidates.
        """
        TIRdict = getallTIRs()
        dictofvals = {}
        print("Number of TIRs: " + str(len(TIRdict)))
        listofaverages = []
        for i in range(0, round(len(self.library) / 2)):  # iterate through all ASDs
            listofvals = []
            ASDname = str('ASD' + str(i + 1))
            for j in range(0, len(TIRdict)):  # for each ASD, iterate through all TIRs in the genome
                TIRname = str('TIR' + str(j + 1))
                val = float(RNAduplexval(self.library[ASDname], TIRdict[TIRname]))
                listofvals.append(val)
            average = sum(listofvals) / len(listofvals)
            dictofvals[average] = ASDname  # calculate the average binding energy between the ASD
            # and all TIRs; here we store the key as the average so that the we can call the names of the highest ASDs after
            # the list is sorted
            listofaverages.append(average)

        listofaverages.sort(reverse=True)
        print('Here are the 10 top candidates with highest ASD-host binding values:')
        for i in range(0, 10):
            print(self.library[str(dictofvals[listofaverages[i]])])


# ------------------------ Create the E. coli instance and run all the steps on it----------------------
initiallib = libbuild() # Create the initial library, this will be replaced by a function that just reads the inital
# library from a file
ecolistrain = Library(initiallib) # Create the E. coli instance.
ecolistrain.narrow_binding() # Step 3
ecolistrain.narrow_crossbinding() # Step 5
ecolistrain.ASD_2rystructure_narrow() # Step 7
ecolistrain.allASDTIRpairs() # Last step
