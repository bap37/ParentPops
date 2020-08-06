EXTRA_FUNC_FILEPATH = 'replaceme'

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

def Optimiser(Param, dfdata, DIM, dfpre, dfpost, binsize): #better explained in Optimiser_Mass
    dfk = dfdata
    if Param == 'c':
        cI_m = MI.Matrix_c_init(dfpre, dfpost, binsize)
    else:
        cI_m = MI.Matrix_x_init(dfpre, dfpost, binsize)
    nwalkers = 2*(DIM + 1)
    ndim = DIM
    p0 = np.random.rand(nwalkers, ndim)
    p0 = p0/100
    p0 = np.abs(p0)
    if Param == 'c':
        sampler = emcee.EnsembleSampler(nwalkers, ndim, Matrix.Matrix_c, args=[dfk, cI_m, binsize])
    else:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, Matrix.Matrix_x, args=[dfk, cI_m, binsize])
    state = sampler.run_mcmc(p0, 1000, progress=True)
    sampler.reset()
    sampler.run_mcmc(state, 10000, progress=True)
    samples = sampler.get_chain()
    print(Param)
    print([np.mean(samples[:,:,0]), np.std(samples[:,:,0])])
    print([np.mean(samples[:,:,1]), np.std(samples[:,:,1])])
    print([np.mean(samples[:,:,2]), np.std(samples[:,:,2])])
    print([np.mean(samples[:,:,3]), np.std(samples[:,:,3])])
    
#Now we need to cycle through masses



""" 

if __name__ == "__main__": 
    Optimiser_Mass()

 """
