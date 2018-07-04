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
