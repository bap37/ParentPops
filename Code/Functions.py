import numpy as np

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

def agauss(x,a,x0,sigma_l,sigma_r,n): #We need to propose an actual input distribution. Here I've chosen an asymmetric Gaussian. We can talk about whether or not this is the right choice, but it provides four different parameters. If that sounds overfitted, I have opinions!
    if x < x0:
        return a*np.exp(-np.abs(x-x0)**n/(n*sigma_l**n))
    else:
        return a*np.exp(-np.abs(x-x0)**n/(n*sigma_r**n))