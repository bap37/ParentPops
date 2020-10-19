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

def ggn(x,a,x0,sigma_l,sigma_r,n): #We need to propose an actual input distribution. Here I've chosen an asymmetric Gaussian. We can talk about whether or not this is the right choice, but it provides four different parameters. If that sounds overfitted, I have opinions!
    if x < x0:
        return a*np.exp(-np.abs(x-x0)**n/(n*sigma_l**n))
    else:
        return a*np.exp(-np.abs(x-x0)**n/(n*sigma_r**n))

def gauss(x,a,x0,sigma): 
    return a*np.exp(-np.abs(x-x0)**2/(2*sigma**2))

#def dgauss(x,a1,x0_1,sigma_1,a2,x0_2, sigma_2): 
#    if x < 0:
#        return a1*np.exp(-np.abs(x-x0_1)**2/(2*sigma_1**2))
#    else:
#        return a2*np.exp(-np.abs(x-x0_2)**2/(2*sigma_2**2))

def dgauss(x,a1,x0_1,sigma_1,a2,x0_2, sigma_2): 
    return a1*np.exp(-np.abs(x-x0_1)**2/(2*sigma_1**2)) + a2*np.exp(-np.abs(x-x0_2)**2/(2*sigma_2**2))

def agauss(x,a,x0,sigma_l,sigma_r): #We need to propose an actual input distribution. Here I've chosen an asymmetric Gaussian. We can talk about whether or not this is the right choice, but it provides four different parameters. If that sounds overfitted, I have opinions!
    if x < x0:
        return a*np.exp(-np.abs(x-x0)**2/(2*sigma_l**2))
    else:
        return a*np.exp(-np.abs(x-x0)**2/(2*sigma_r**2))

def poisson_interval(k, alpha=0.32):
    """                                                                                                                                                                                             
    uses chisquared info to get the poisson interval. Uses scipy.stats                                                                                                                              
    (imports in function).                                                                                                                                                                          
    (http://stackoverflow.com/questions/14813530/poisson-confidence-interval-with-numpy)                                                                                                            
    """
    from scipy.stats import chi2
    a = alpha
    low, high = (chi2.ppf(a/2, 2*k) / 2, chi2.ppf(1-a/2, 2*k + 2) / 2)
    low[k == 0] = 0.0
    #if k == 0:                                                                                                                                                                                     
    #    low = 0.0                                                                                                                                                                                  
    return low, high