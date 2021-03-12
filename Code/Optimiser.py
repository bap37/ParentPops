EXTRA_FUNC_FILEPATH = 'replace me'

import pandas as pd
import numpy as np
import sys
import Functions #I will need to include this
import MI
import Matrix
import argparse
parser=argparse.ArgumentParser()
import emcee
from scipy.stats import binned_statistic
import os
import csv     
#If we want redshift, it'll need to be in log scale - get back to you on that

class Optimizer_Calculation:

    def optimize(self, Param, dfdata, SHAPE2, Migration_Matrix, binsize, SHAPE, STEP, CORR, BIN_THINGS):
        dfk = dfdata
        return_dictionary = {}

        if STEP != None: #Renamed MASS to STEP, since it's no longer just Mass
            if (CORR == 'zHD') or (CORR == 'zCMB'):
                lowbound, highbound = STEP
                dfk = dfk.loc[(dfk[CORR] >= lowbound) & (dfk[CORR] <=highbound)]
                STEP = np.mean(lowbound+highbound)

            else:
                window = BIN_THINGS[3] #Final entry in BIN_THINGS should be window value.
                dfk = dfk.loc[(dfk[CORR] < STEP + window) & (dfk[CORR] >= STEP - window)]
            
            #Minimum length: There's no point calculating if the number of supernovae is below 50
            if len(dfk) < 50:
                print(len(dfk), "Supernovae in this bin")
                print('Not enough supernova in this bin:', np.around(STEP,3))
                print("Skipping for now...")
                return return_dictionary

            return_dictionary[CORR] = str(round(STEP,3))
        else:
            return_dictionary[CORR] = "ALL"

        print("We have this many supernova in the current bin:", len(dfk))
        DIM = len(SHAPE2)
        if Param == 'c':
            cI_m = Migration_Matrix
            if DIM == 4: #transfer this into optimizer method. 
                nwalkers = 2*(DIM + 1) - 1
                ndim = DIM - 1 
            else: 
                nwalkers = 2*(DIM + 1) 
                ndim = DIM 
            p0 = np.random.rand(nwalkers, ndim)
            p0 = p0/100
            p0 = np.abs(p0)    
        else:
            cI_m = Migration_Matrix
            if SHAPE == 'GGN': #transfer this into optimizer method.  
                nwalkers = 2*(DIM + 1) - 1
                ndim = DIM - 1 
            else:
                nwalkers = 2*(DIM + 1)
                ndim = DIM 
            if SHAPE == 'DGaussian':
                p0 = np.random.rand(nwalkers, ndim)
                p0 = np.random.rand(nwalkers, ndim)      
                #p0 = p0/100
                p0 = np.abs(p0) 
                for q in range(len(p0[:,0])):
                    p0[q,1] = np.random.normal(-1,.1,1)
                    p0[q,4] = np.random.normal(1,.1,1)
                    p0[q,0] = np.random.normal(1,.1,1)
                    p0[q,3] = np.random.normal(1,.1,1)
            else:
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
            return_dictionary[SHAPE2[h]] = [np.mean(samples[:,:,h]), np.std(samples[:,:,h])]
        
        return return_dictionary
    
    def optimize_in_range(self, Param, dfdata, SHAPE2, Migration_Matrix, binsize, SHAPE, CORR, BIN_THINGS):
        collected_result = []
        
        if (CORR == 'zHD') or (CORR == 'zCMB'):
            dfdata = dfdata.sort_values(by=['zHD'])       
            ranger = dfdata.zHD.values[::BIN_THINGS[2]]    
#            for q in range(len(ranger)):       
#                try:     
#                    print(ranger[q], ranger[q+1])   
#                except IndexError:       
#                    print(ranger[q-1], dfdata.zHD.values[-1]) 
            TOTAL_BINS = ranger

        else:
            TOTAL_BINS = np.arange(BIN_THINGS[0], BIN_THINGS[1], BIN_THINGS[2]) 

        for m in range(len(TOTAL_BINS)):
            if (CORR == 'zHD') or (CORR == 'zCMB'):
                try:
                    _ = TOTAL_BINS[m+1]
                    calculation = self.optimize(Param, dfdata, SHAPE2, Migration_Matrix, binsize, SHAPE, (TOTAL_BINS[m], TOTAL_BINS[m+1]), CORR, BIN_THINGS)
                except IndexError:
                    pass
            else:
                calculation = self.optimize(Param, dfdata, SHAPE2, Migration_Matrix, binsize, SHAPE, TOTAL_BINS[m], CORR, BIN_THINGS)
            #def optimize(self, Param, dfdata, SHAPE2, Migration_Matrix, binsize, SHAPE, STEP, CORR, BIN_THINGS):
            #A return of an empty dictionary means that specific mass doesn't have enough supernovae
            #to have been worth calculating on. So it's not part of the output.
            if len(calculation) == 0:
                pass
            else:
                collected_result.append(calculation)
            
        
        return collected_result
    
    def write_to_file(self, result_dictionary, SHAPE2, SURVEY, TYPE, SHAPE, MODEL, is_series, Param, REF_FP, CORR):
        if not os.path.exists(REF_FP+"output"):
            os.mkdir(REF_FP+"output")

        file_name = REF_FP + "output/" + TYPE + "_" + SHAPE + "_" +  MODEL + "_" + Param + "_" + CORR
        for survey in SURVEY:
            file_name += "_" + survey
        file_name += ".tsv" 

        with open(file_name, mode='w', newline='\n') as outFile:
            outFileWriter = csv.writer(outFile, delimiter = '\t')

            header = [CORR]
            for column in SHAPE2:
                header.append(column)
                header.append(column + '_error')
            outFileWriter.writerow(header)

            if(is_series):
                # We are dealing with multiple iterations
                for result in result_dictionary:
                    outFileWriter.writerow(self.row_from_dictionary(result, SHAPE2, Param, CORR))

            else:
                outFileWriter.writerow(self.row_from_dictionary(result_dictionary, SHAPE2, Param, CORR))

        print("Written To File")
    
    def row_from_dictionary(self, dictionary, SHAPE2, Param, CORR):
        row = [dictionary[CORR]]
        for shape in SHAPE2:
            try:
                for result in dictionary[shape]:
                    row.append(str(result))
            except KeyError:
                    for rest in range(2): #This absolutely needs to be fixed.
                        if rest == 0:
                            if Param == 'c':
                                row.append(str(2))
                            else:
                                row.append(str(3))
                        else:
                            row.append(str(0))

            #for result in dictionary[shape]:
             #   row.append(str(result))
        return row


    def write_MM(self, MM, SHAPE2, SURVEY, TYPE, SHAPE, MODEL, Param, REF_FP, CORR):         
        if not os.path.exists(REF_FP+"output"):
            os.mkdir(REF_FP+"output")   
        file_name = REF_FP + "output/" + TYPE + "_" + SHAPE + "_" +  MODEL + "_" + Param + "_" + CORR    
        for survey in SURVEY:                                   
            file_name += "_" + survey                  
        file_name += "_MM"  
        np.save(file_name, MM)
        print("Saved Migration Matrix")
