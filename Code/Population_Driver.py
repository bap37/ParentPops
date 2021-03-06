#!/usr/bin/env python
#Define config variables up here
#Right now this only works when run in command line with the relevant arguments.
import os
if os.getenv("REF_FILEPATH") is None:
    print("Please define $REF_FILEPATH as an environmental variable. It needs to point to where your input files are.")
    print("REF_FILEPATH=", os.getenv("REF_FILEPATH"))
    quit()
else:
    REF_FP = os.getenv("REF_FILEPATH") + "/"

#Ideally we'd like to be able to load a couple different options here.
#Select arbitrary combination of SURVEY/TYPE/MODEL
#SURVEY = DES, SDSS, PS1, SNLS, FOUND, LOWZ
#TYPE = PHOT, SPEC
#MODEL = G10, C11
#keep in mind that SPEC exists for every survey, but PHOT does not - need error trapping there
#Also, LOWZ needs to be handled separately always. 

import pandas as pd
import numpy as np
import sys
#sys.path.insert(0, EXTRA_FUNC_FILEPATH) 
import Functions 
import MI
import Matrix
import matplotlib.pyplot as plt

import argparse
parser=argparse.ArgumentParser()
import emcee
from scipy.stats import binned_statistic



#Commented options don't work right now.

#Most important parameters are SURVEY, Param, MODEL, TYPE

parser.add_argument('--SURVEY', help='Choose the Survey', nargs='+') #It would be nice to get this to work. To switch from the cumulative to a specific survey on command.
parser.add_argument('--PARAM', help='x1 or c')
parser.add_argument('--I_c', help='Initial guess for c', default=[0,.1,.1,2], nargs = '+')
parser.add_argument('--I_x1', help='Initial guess for x1', default=[0,1,1,2], nargs = '+')
parser.add_argument('--SHAPE', help='How do you want to model these parent populations? GGN, Gaussian, DGaussian?', default='GGN')
#parser.add_argument('--Choice', help='Which bin you wish to investigate') LOW priority, but would be nice to select a specific range to investigate
#parser.add_argument('--Redshift', help='If you want to investigate redshift bins rather than mass', default=False) MEDIUM priority. Not the point of the paper, but super easy to implement. Just change masses to redshifts.
#parser.add_argument('--Test', help='Tests an iteration, more print statements', default=False) Junk.
#parser.add_argument('--Total', help='Cycle through all the surveys', default=False) MEDIUM Priority. Right now, it only does one survey per run. A flag to cycle through all surveys at once would be useful, but not necessary.
parser.add_argument('--MODEL', help='Toggle to True if using the C11 scatter model')  
#parser.add_argument('--Iteration', help='Are you using a Gaussian distribution? If so, this should be true.', default=False) #junk.
parser.add_argument('--TYPE', help='Is your sample spectroscopic? default= No.') #HIGH Priority. Photometry is the future of SNIa cosmology
parser.add_argument('--ITERATIVE', help='Make this option False if you want to marginalise over MASS', default=True)
parser.add_argument('--REF', help='Set a custom REF directory.', default=os.getenv("REF_FILEPATH") + "/")
#A choice of distribution shapes is desired but difficult. Can discuss more. 
parser.add_argument('--BINS', help="""For use with custom external properties.

First bin is the minimum value.
Second bin is the maximum value.
Third bin is the step size (for np.arange), UNLESS you choose zHD or zCMB. If using redshift, we select a certain number of events per bin, so this will need to be the number of SNe per bin that you want. Needs to be over 50.
Fourth bin is window size. 

""", default=[6,14,0.2,0.6], nargs= '+')
parser.add_argument('--CORR', help='Name of the parameter you want to investigate correlations with. Default is HOST_LOGMASS.', default="HOST_LOGMASS")

parser.add_argument('--MM', help="Will save the associated Migration Matrix, but not run the program. Default is False.", default=False, type=bool)

args = parser.parse_args()
Param = args.PARAM
I_c = args.I_c
I_c = [float(i) for i in I_c]
I_x1 = args.I_x1
I_x1 = [float(i) for i in I_x1]
#choice = args.choice
#redshift = args.redshift
SURVEY= args.SURVEY
MODEL = args.MODEL
TYPE = args.TYPE
SHAPE = args.SHAPE
IT = args.ITERATIVE
REF_FP = args.REF
BIN_THINGS = args.BINS 
BIN_THINGS = [float(i) for i in BIN_THINGS]
BIN_THINGS = [np.around(i,4) for i in BIN_THINGS]
CORR = args.CORR
MM = args.MM

import distutils.util
IT = distutils.util.strtobool(IT)

xbinsize=.2
cbinsize=0.02

if CORR != 'HOST_LOGMASS': 
    if BIN_THINGS == [6.0, 14.0, 0.2, 0.6]:
        print("The bins you gave me are the default for HOST_LOGMASS. But you did not specify to use HOST_LOGMASS. Please change your bins; check the help string for --BINS for more info.")
        print("quitting")
        quit()


if SURVEY == ['HIZ']: #if we want to do all targeted surveys at once.
    SURVEY = ['DES', 'SNLS', 'SDSS', 'PS1']


IDSURVEY_Dictionary = {1:'SDSS', 4:'SNLS', 10:'DES', 15:'PS1', 150:'FOUND'}

SH_DIC = {'GGNN': ['means', 'stdl', 'stdr', 'n'], 'Gaussian':['means', 'std'], 
'DGaussian':['a1', 'mean1', 'std1', 'a2', 'mean2', 'std2'], 'AGaussian':['means', 'stdl', 'stdr'], 'GGN': ['means', 'stdl', 'stdr', 'nfake']}
SHAPE2 = SH_DIC[SHAPE]
import distutils.util

if SURVEY is None:
    print("You haven't specified a Survey!")
    quit()


 
print('You are looking at '+str(SURVEY)+' with the '+str(MODEL)+' model and confirmed as '+str(TYPE))


# ['DES', 'SDSS', 'FOUND', 'SNLS', 'PS1'] The five high-redshift surveys. 
matrixdic = {} #store the probability matrix for each survey

lensdic = {} #Store relative contributions of each survey. SDSS has twice as many SNIa as DES, for instance. Need to account for that. 

#Param = 'c'
count = 0
for s in SURVEY:
    #Everything is in FITOPTS. From there the structure is SURVEY_TYPE_MODEL.
    #SURVEY is what's in slist, line 106
    #TYPE is SPEC or PHOT, depending. Not all surveys have both, so some error trapping is necessary.
    #MODEL is G10 or C11.
    #Right now this is hardcoded for G10, which is unlabeled, but something that will be fixed when this moves off of my computer and into a more public space. 
    filename1 = REF_FP+s+'_'+TYPE+'_'+MODEL+'/MATRIX_PRE.FITRES'
    filename2 = REF_FP+s+'_'+TYPE+'_'+MODEL+'/MASS_0_DATA.FITRES'
    filename3 = REF_FP+s+'_'+TYPE+'_'+MODEL+'/MATRIX_POST.FITRES'
    try:
        Names1, StartRow1 = Functions.NAndR(filename1)
        Names2, StartRow2 = Functions.NAndR(filename2)
        Names3, StartRow3 = Functions.NAndR(filename3)
    except FileNotFoundError or NameError:
        print('Probably forgot a slash somewhere, check your filepath')
        print("Your current filepath is: ", REF_FP)
        sys.stdout.flush()
        quit()
        
    #This block needs to be run for every SURVEY+TYPE+MODEL. 
    print('Loading True Simulated Population...')
    sys.stdout.flush()
    dfpre = pd.read_csv(filename1, header=None, skiprows=StartRow1,names=Names1, delim_whitespace=True, skip_blank_lines=True, error_bad_lines=False, dtype={'S2c':float, 'S2x1':float, 'CID': int}, comment='#')
    print('Loading Observed Simulated Population...')
    sys.stdout.flush()
    dfdata = pd.read_csv(filename2, header=None, skiprows=StartRow2,names=Names2, delim_whitespace=True, skip_blank_lines=True, error_bad_lines=False, comment='#')
    print('Loading Real Data...')
    sys.stdout.flush()
    dfpost = pd.read_csv(filename3, header=None, skiprows=StartRow3,names=Names3, delim_whitespace=True, skip_blank_lines=True, error_bad_lines=False, comment='#')
    #dfpre.CID = pd.to_numeric(dfpre.CID, errors='coerce') #sometimes the CIDs are stored as strings. Don't want that. 
    dfpre = dfpre.loc[dfpre.CID == dfpre.CID] #Getting rid of nans
    #dfpre.CID = dfpre.CID.astype(int) #setting type. 
    #dfpre.S2c = pd.to_numeric(dfpre.S2c, errors='coerce') #S2c is the "true" colour. The simulation chooses a "true" value at the start of the simulation. The process of observing the supernova, measurement noise, ambient noise, etc, can make the "true" colour not equal to the observed one. 
    #dfpre.S2x1 = pd.to_numeric(dfpre.S2x1, errors='coerce') #Likewise, this is the "true" stretch.
    #I think I've upgraded the loading procedures enough that the forcing numeric can be commented out. Leaving it just in case.

    if CORR not in list(dfdata): #Then we make sure everything's a valid parameter.    
        print('You need to give a valid variable to look at! You gave me', CORR) 
        print('Please check your FITRES file. Quitting...')    
        quit() 

    print(str(100*len(dfdata.loc[dfdata.HOST_LOGMASS > 3])/len(dfdata))+"% of the sample has valid mass...") #Just reading out the number of supernovae in each survey
    sys.stdout.flush()
    try: 
        print('That was '+str(IDSURVEY_Dictionary[np.unique(dfdata.IDSURVEY.values)[0]])+' that you just loaded!') #printing the SURVEY ID. 
        sys.stdout.flush()
    except KeyError:
        if s == 'LOWZ':
            print('LOWZ')
            sys.stdout.flush()
        else:
            print("Woops! This doesn't correspond to a valid IDSURVEY! Are you loading files correctly?")
            sys.stdout.flush()
            quit()

    except TypeError:
        print("uhhh, looks like something's up with the IDSURVEY value in your FITRES file! I can't look up the value in the name dictionary. This doesn't mean anything's gone terribly wrong, but it's worth checking out!")
        sys.stdout.flush()

    if Param == 'c':
        matrixdic[s+Param] = MI.Matrix_c_init(dfpre, dfpost, cbinsize)
        binsize  = cbinsize
    else:
        matrixdic[s+Param] = MI.Matrix_x_init(dfpre, dfpost, xbinsize)
        binsize = xbinsize
    if count == 0:
        DATOT = dfdata[['CID','x1', 'c', CORR]] #The important parameters.
    else:
        DATOT = DATOT.append(dfdata[['CID','x1', 'c', CORR]])
    if CORR == 'HOST_LOGMASS':
        lensdic[s+Param] = len(dfdata.loc[dfdata.HOST_LOGMASS > 3])
    else:
        lensdic[s+Param] = len(dfdata)
    count += 1 #I think this is fairly obvious?


q = 0
for s in SURVEY:
    q += lensdic[s+Param]

for num,s in enumerate(SURVEY):
    if num == 0:
        newmatrix = matrixdic[s+Param]*(lensdic[s+Param] / q) #Here we're making sure each survey gives a proportional contribution
    else:
        newmatrix += matrixdic[s+Param]*(lensdic[s+Param] / q)

import Optimiser

if SURVEY == ['DES', 'SNLS', 'SDSS', 'PS1']:
    SURVEY = ['HIZ']


optimizer = Optimiser.Optimizer_Calculation()

if MM == True:
    optimizer.write_MM(newmatrix, SHAPE2, SURVEY, TYPE, SHAPE, MODEL, Param, REF_FP, CORR)
    print("Quitting... This is an expected closure.")
    quit()

if IT == True:
    result = optimizer.optimize_in_range(Param, DATOT, SHAPE2, newmatrix, binsize, SHAPE, CORR, BIN_THINGS)
    optimizer.write_to_file(result, SHAPE2, SURVEY, TYPE, SHAPE, MODEL, True, Param, REF_FP, CORR)
elif IT == False:
    result = optimizer.optimize(Param,DATOT, SHAPE2, newmatrix, binsize, SHAPE, None, CORR, BIN_THINGS)
    paramslist = []
    for vals in result:
        paramslist.append(result[vals][0])
    paramslist = paramslist[1:]
    [float(i) for i in paramslist]
    paramslist = np.array(paramslist)
    if Param == 'x1':
        LL, plotPredicted, plotData, plotbins = Matrix.Matrix_x(paramslist, DATOT, newmatrix, xbinsize, SHAPE, debug=True)
    else:
        LL, plotPredicted, plotData, plotbins = Matrix.Matrix_c(paramslist, DATOT, newmatrix, cbinsize, SHAPE, debug=True)
    errl,erru = Functions.poisson_interval(plotData)
    plt.figure()
    plt.errorbar(plotbins, plotData, yerr=[plotData-errl, erru-plotData], label='Data', c='k', fmt='o')
    plt.plot(plotbins, plotPredicted, label='Predicted', drawstyle='steps-mid')
    plt.xlabel(Param)
    plt.legend()
    plt.ylabel('Count')
    plt.title(TYPE + "_" + SHAPE + "_" +  MODEL + "_" + Param)
    plt.savefig(REF_FP + "output/" + TYPE + "_" + SHAPE + "_" +  MODEL + "_" + Param + ".pdf", format='pdf')
    print(np.sum(((plotData-plotPredicted)**2)/(erru-plotData)**2))
    optimizer.write_to_file(result, SHAPE2, SURVEY, TYPE, SHAPE, MODEL, False, Param, REF_FP, CORR)
else:
    print('oops, you hecked up!')
