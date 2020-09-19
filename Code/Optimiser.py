EXTRA_FUNC_FILEPATH = 'replace me'

import pandas as pd
import numpy as np
import sys
sys.path.insert(0, EXTRA_FUNC_FILEPATH) 
import Functions #I will need to include this
import MI
import Matrix
import argparse
parser=argparse.ArgumentParser()
import emcee
from scipy.stats import binned_statistic

def Optimiser(Param, dfdata, SHAPE2, dfpre, dfpost, binsize, SHAPE): #better explained in Optimiser_Mass
    dfk = dfdata
    DIM = len(SHAPE2)
    if Param == 'c':
        cI_m = MI.Matrix_c_init(dfpre, dfpost, binsize)
        if DIM == 4:
            nwalkers = 2*(DIM + 1) - 1
            ndim = DIM - 1 
        else: 
            nwalkers = 2*(DIM + 1) 
            ndim = DIM 
        p0 = np.random.rand(nwalkers, ndim)
        p0 = p0/100
        p0 = np.abs(p0)    
    else:
        cI_m = MI.Matrix_x_init(dfpre, dfpost, binsize)
        nwalkers = 2*(DIM + 1)
        ndim = DIM
        p0 = np.random.rand(nwalkers, ndim)
        p0 = p0/100
        p0 = np.abs(p0)
    if Param == 'c':
        sampler = emcee.EnsembleSampler(nwalkers, ndim, Matrix.Matrix_c, args=[dfk, cI_m, binsize, SHAPE])
    else:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, Matrix.Matrix_x, args=[dfk, cI_m, binsize, SHAPE])
    state = sampler.run_mcmc(p0, 1000, progress=True)
    sampler.reset()
    sampler.run_mcmc(state, 10000, progress=True)
    samples = sampler.get_chain()
    print(Param)
    for h in range((ndim)):
        print([np.mean(samples[:,:,h]), np.std(samples[:,:,h])], SHAPE2[h])

    
#Now we need to cycle through masses
#masses = np.arange(8, 13.6, 0.2)
#With window of 0.6 right now - eg 0.2 +/- 0.6 

#If we want redshift, it'll need to be in log scale - get back to you on that



""" 

if __name__ == "__main__": 
    Optimiser_Mass()

 """
