#EXTRA_FUNC_FILEPATH = 'replaceme'

import pandas as pd
import numpy as np
import sys
#sys.path.insert(0, EXTRA_FUNC_FILEPATH) 
import Functions #I will need to include this
import MI
import Matrix
import argparse
parser=argparse.ArgumentParser()
import emcee
from scipy.stats import binned_statistic

def Optimiser(Param, dfdata, DIM, dfpre, dfpost, binsize, MASS): #better explained in Optimiser_Mass
    dfk = dfdata
    if MASS != None:
        window = 0.6
        dfk = dfk.loc[(dfk.HOST_LOGMASS < MASS +window) & (dfk.HOST_LOGMASS >= MASS - window)]
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
    print([np.mean(samples[:,:,0]), np.std(samples[:,:,0])]) # 0 == mean
    print([np.mean(samples[:,:,1]), np.std(samples[:,:,1])]) # stdl
    print([np.mean(samples[:,:,2]), np.std(samples[:,:,2])]) # stdr
    print([np.mean(samples[:,:,3]), np.std(samples[:,:,3])]) # n
    # TODO: Tab deliniation csv? If do arange == true, append to the tsv output of each iteration (m is default HOST_LOGMASS from dfdata, otherwise the arrange)
    # TODO: mean, stdl, stdr, n, m
    
#Now we need to cycle through masses

class Optimizer_Calculation:
    
    def optimize(self, Param, dfdata, DIM, dfpre, dfpost, binsize, MASS):
        dfk = dfdata
        if MASS != None:
            window = 0.6
            dfk = dfk.loc[(dfk.HOST_LOGMASS < MASS +window) & (dfk.HOST_LOGMASS >= MASS - window)]
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
        print([np.mean(samples[:,:,0]), np.std(samples[:,:,0])]) # 0 == mean
        print([np.mean(samples[:,:,1]), np.std(samples[:,:,1])]) # stdl
        print([np.mean(samples[:,:,2]), np.std(samples[:,:,2])]) # stdr
        print([np.mean(samples[:,:,3]), np.std(samples[:,:,3])]) # n
        # TODO: Tab deliniation csv? If do arange == true, append to the tsv output of each iteration (m is default HOST_LOGMASS from dfdata, otherwise the arrange)
        # TODO: mean, stdl, stdr, n, m
        return [[],[],[],[],[]]
    
    def write_to_file(self):
        print("Written To File")

    def optimize_in_range(self, Param, dfdata, DIM, dfpre, dfpost, binsize, MASS):
        mean = []
        stdl= []
        stdr = []
        n = []
        m = []
        for m in np.arange(6,14,.2):
            calculation = self.optimize(Param, dfdata, DIM, dfpre, dfpost, binsize, m) # RETURNS: An array of arrays, for mean, stdl, stdr, n, m
            mean = mean.append(calculation[0])
            stdl = stdl.append(calculation[1])
            stdr = stdr.append(calculation[2])
            n = n.append(calculation[3])
            m = m.append(calculation[4])
        return [mean, stdl, stdr, n, m] #Allows us potentially to print from here



        
        


""" 

if __name__ == "__main__": 
    Optimiser_Mass()

 """
