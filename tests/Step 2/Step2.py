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
