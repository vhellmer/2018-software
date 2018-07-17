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

