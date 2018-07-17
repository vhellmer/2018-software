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
