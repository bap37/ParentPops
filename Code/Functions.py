def NAndR(filename):
    """Takes in a FITRES file and outputs the variable names and startline for the data

    Outputs are [Names, Startrow]. List and number respectively.
    """
    with open(filename) as fp:
        for i, line in enumerate(fp):
            if line.startswith('VARNAMES:'):
                line = line.replace(',',' ')
                line = line.replace('\n','')
                Names = line.split()
            elif line.startswith('SN'):
                Startrow = i
                break
    return Names, Startrow

    