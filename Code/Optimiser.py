EXTRA_FUNC_FILEPATH = 'replace me'

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
