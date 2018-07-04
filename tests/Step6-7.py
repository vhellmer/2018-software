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
