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
