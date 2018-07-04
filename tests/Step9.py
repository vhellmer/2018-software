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
