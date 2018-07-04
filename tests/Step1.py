def convertTtoU(sequence):
    """
    Convert thymine to uracil

    :param sequence: The DNA sequence that is being converted to mRNA, in string form
    :return: Another string with T's replaced with U's.
    """
    return sequence.replace('T', 'U')
