#Contains the functions needed to find the chi2. 
EXTRA_FUNC_FILEPATH = 'replace me'

import numpy as np
import pandas as pd
import sys
#sys.path.insert(0, EXTRA_FUNC_FILEPATH) 
import Functions #I will need to include this

def Matrix_c(params, dfk, cI_m, binsize, SHAPE2, debug=False):
    bins = np.arange(-.4,.41,binsize)
    cdatI = np.histogram(dfk.c.values, bins=bins)[0] #A histogram of the data
    input_cI = []
    if SHAPE2 == 'GGN':    
        cI_mean, cI_l, cI_r = params
        cI = [1, cI_mean, cI_l, cI_r]
        
        for x in ((bins[1:] + bins[:-1])/2):
            input_cI.append(Functions.agauss(x, *cI))
        input_cI = np.array(input_cI)
    elif SHAPE2 == 'Gaussian':
        cI_mean, cI_std = params
        cI = [1, cI_mean, cI_std]
        for x in ((bins[1:] + bins[:-1])/2):
            input_cI.append(Functions.gauss(x, *cI))
        input_cI = np.array(input_cI)
    elif SHAPE2 == 'DGaussian':
        a1, mean1, std1, a2, mean2, std2 = params
        cI = [a1, mean1, std1, a2, mean2, std2]
        for x in ((bins[1:] + bins[:-1])/2):
            input_cI.append(Functions.dgauss(x, *cI))
        input_cI = np.array(input_cI)


    
    MP = np.matmul(input_cI, cI_m) #We take our true populaiton and process it through the probability matrix. This tells us what it should look like after being observed in a real telescope. Roughly.
    MP = MP*((np.sum(cdatI))/np.sum(MP)) #normalising it.
 
    Delta_c = cdatI - MP #This is detailed better in Scolnic and Kessler 2016, the paper that's attached in the email you got. But this is the difference between the data distribution and our guess.
    Delta_cerr = np.copy(cdatI) #we use Poisson statistics, so the error is the square root of the occupation of the bins. The data isn't normalised. 
    Delta_cerr[Delta_cerr == 0 ] = 1 #setting the unpopulated bins to 1 to avoid nans
    Delta_cerr = np.sqrt(np.abs(Delta_cerr)) #taking the square root.
    chi2 = []
    for i in range(len(Delta_c)):
        temp = ((Delta_c[i]))/((Delta_cerr[i])) #This is leftover from an error flagging protocol I had for a while. Certain bins were causing some huge statistcal biases, so this will ignore them in the chi2 calculation. flag always starts as 1.
        temp = temp**2
        chi2.append(temp) #It's a chi2.
    chi2 = np.array(chi2) 
    LL = -np.sum(chi2)/2. #Converted to log-likelihood. The closer it is to 0 the better. Negative values likely.
    #these if statements are designed to discourage the MCMC from certain values. 
    if LL != LL:
        LL = -np.inf #no nans.
    if (cI_l < 0) or (cI_r < 0): #Don't want negative standard deviations.
        LL = -np.inf
    if (cI_l > .3) or (cI_r > .3): # |0.3| is the max allowed value for c in cosmological analysis. if the std is 0.3, that's flat. No good. 
        LL = -np.inf
    if (np.abs(cI_mean) > .3): #same as standard deviation.
        LL = -np.inf
    if debug == True:
        return LL, (Delta_c/Delta_cerr)**2 #I did a lot of debugging.
    else:
        return LL



def Matrix_x(params, dfk, xI_m, binsize, SHAPE2, debug=False): #same as Matrix_c, but for stretch.
    bins = np.arange(-3,3.1,binsize)
    xdatI = np.histogram(dfk.x1.values, bins=bins)[0]
    input_xI = []
    if SHAPE2 == 'GGN':    
        xI_mean, xI_l, xI_r, xI_n = params
        xI = [1, xI_mean, xI_l, xI_r, xI_n]
        
        for x in ((bins[1:] + bins[:-1])/2):
            input_xI.append(Functions.ggn(x, *xI))
        input_xI = np.array(input_xI)
    elif SHAPE2 == 'Gaussian':
        xI_mean, xI_std = params
        xI = [1, xI_mean, xI_std]
        for x in ((bins[1:] + bins[:-1])/2):
            input_xI.append(Functions.gauss(x, *xI))
        input_xI = np.array(input_xI)
    elif SHAPE2 == 'DGaussian':
        a1, mean1, std1, a2, mean2, std2 = params
        xI = [a1, mean1, std1, a2, mean2, std2]
        for x in ((bins[1:] + bins[:-1])/2):
            input_xI.append(Functions.dgauss(x, *xI))
        input_xI = np.array(input_xI)
    
    MP = np.matmul(input_xI, xI_m)
    MP = MP*((np.sum(xdatI))/np.sum(MP))
    
    Delta_x = xdatI - MP
    #print(Delta_c)
    Delta_xerr = np.copy(xdatI)
    Delta_xerr[Delta_xerr == 0 ] = 1
    #Delta_xerr[xdatI < 5] = 20
    Delta_xerr = np.sqrt(np.abs(Delta_xerr))
    chi2 = []
    for i in range(len(Delta_x)):
        temp = ((Delta_x[i]))/((Delta_xerr[i]))
        temp = temp**2
        chi2.append(temp)
    chi2 = np.array(chi2)
    LL = -np.sum(chi2)/2.
    if LL != LL: 
        LL = -np.inf
    if SHAPE2 == 'GGN': #if statements sorted by thematic consistency 
        if (xI_l < 0) or (xI_r < 0):
            LL = -np.inf
        if (xI_l > 3) or (xI_r > 3):
            LL = -np.inf
        if (np.abs(xI_mean) > 3):
            LL = -np.inf
        if (np.abs(xI_n) > 4):
            LL = -np.inf
    elif SHAPE2 == 'DGaussian': #a1, mean1, std1, a2, mean2, std2 = params
        if (std1 < 0) or (std2 < 0):
            LL = -np.inf
        if (std1 > 3) or (std2 > 3):
            LL = -np.inf
        if (np.abs(mean1) > 3) or (np.abs(mean2) > 3):
            LL = -np.inf
        if (a1 < 0 ) or (a2 < 0) or (a1 > 5) or (a2 > 5):
            LL = -np.inf
    if debug == True:
        return LL, (Delta_x/Delta_xerr)**2
    else:
        return LL
