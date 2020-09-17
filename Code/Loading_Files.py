#Define config variables up here
#Right now this only works when run in command line with the relevant arguments.

#EXTRA_FUNC_FILEPATH = 'replaceme' #'I need to point to ParentPop/Code/'
REF_FP = 'replaceme' #I need to point to the unzipped version of REF.zip!'

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


import argparse
parser=argparse.ArgumentParser()
import emcee
from scipy.stats import binned_statistic

#Commented options don't work right now.

parser.add_argument('--SURVEY', help='Choose the Survey', nargs='+') #It would be nice to get this to work. To switch from the cumulative to a specific survey on command.
parser.add_argument('--Param', help='x1 or c')
parser.add_argument('--I_c', help='Initial guess for c', default=[0,.1,.1,2], nargs = '+')
parser.add_argument('--I_x1', help='Initial guess for x1', default=[0,1,1,2], nargs = '+')
#parser.add_argument('--Choice', help='Which bin you wish to investigate') LOW priority, but would be nice to select a specific range to investigate
#parser.add_argument('--Redshift', help='If you want to investigate redshift bins rather than mass', default=False) MEDIUM priority. Not the point of the paper, but super easy to implement. Just change masses to redshifts.
#parser.add_argument('--Test', help='Tests an iteration, more print statements', default=False) Junk.
#parser.add_argument('--Total', help='Cycle through all the surveys', default=False) MEDIUM Priority. Right now, it only does one survey per run. A flag to cycle through all surveys at once would be useful, but not necessary.
parser.add_argument('--MODEL', help='Toggle to True if using the C11 scatter model')  
#parser.add_argument('--Iteration', help='Are you using a Gaussian distribution? If so, this should be true.', default=False) #junk.
parser.add_argument('--TYPE', help='Is your sample spectroscopic? default= No.') #HIGH Priority. Photometry is the future of SNIa cosmology
#A choice of distribution shapes is desired but difficult. Can discuss more. 

args = parser.parse_args()
Param = args.Param
I_c = args.I_c
I_c = [float(i) for i in I_c]
I_x1 = args.I_x1
I_x1 = [float(i) for i in I_x1]
#choice = args.choice
#redshift = args.redshift
SURVEY= args.SURVEY
MODEL = args.MODEL
TYPE = args.TYPE

xbinsize=.2
cbinsize=0.02


if SURVEY == 'HZ': #if we want to do all targeted surveys at once.
    SURVEY = ['DES', 'SNLS', 'SDSS', 'PS1', 'FOUND']



IDSURVEY_Dictionary = {1:'SDSS', 4:'SNLS', 10:'DES', 15:'PS1', 150:'FOUND'}


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
    try:
        filename1 = REF_FP+s+'_'+TYPE+'_'+MODEL+'/MATRIX_PRE.FITRES'
        filename2 = REF_FP+s+'_'+TYPE+'_'+MODEL+'/MASS_0_DATA.FITRES'
        filename3 = REF_FP+s+'_'+TYPE+'_'+MODEL+'/MATRIX_POST.FITRES'
        Names1, StartRow1 = Functions.NAndR(filename1)
        Names2, StartRow2 = Functions.NAndR(filename2)
        Names3, StartRow3 = Functions.NAndR(filename3)
    except FileNotFoundError:
        if TYPE == 'PHOT':
            print(REF_FP+SURVEY+TYPE+MODEL)
            print("That file does not exist! Are you sure there's a photometric sample available!")
        
    #This block needs to be run for every SURVEY+TYPE+MODEL. 
    print('Loading True Population...')
    dfpre = pd.read_csv(filename1, header=None, skiprows=StartRow1,names=Names1, delim_whitespace=True, skip_blank_lines=True, error_bad_lines=False, dtype={'S2c':float, 'S2x1':float, 'CID': int}, comment='#')
    print('Loading Observed Population...')
    dfdata = pd.read_csv(filename2, header=None, skiprows=StartRow2,names=Names2, delim_whitespace=True)
    print('Loading Real Data...')
    dfpost = pd.read_csv(filename3, header=None, skiprows=StartRow3,names=Names3, delim_whitespace=True)
    #dfpre.CID = pd.to_numeric(dfpre.CID, errors='coerce') #sometimes the CIDs are stored as strings. Don't want that. 
    dfpre = dfpre.loc[dfpre.CID == dfpre.CID] #Getting rid of nans
    #dfpre.CID = dfpre.CID.astype(int) #setting type. 
    #dfpre.S2c = pd.to_numeric(dfpre.S2c, errors='coerce') #S2c is the "true" colour. The simulation chooses a "true" value at the start of the simulation. The process of observing the supernova, measurement noise, ambient noise, etc, can make the "true" colour not equal to the observed one. 
    #dfpre.S2x1 = pd.to_numeric(dfpre.S2x1, errors='coerce') #Likewise, this is the "true" stretch.
    #I think I've upgraded the loading procedures enough that the forcing numeric can be commented out. Leaving it just in case.

    print(str(len(dfdata.loc[dfdata.HOST_LOGMASS > 3])/len(dfdata))+"% of the sample can be used...") #Just reading out the number of supernovae in each survey

    try: 
        print('That was '+str(IDSURVEY_Dictionary[np.unique(dfdata.IDSURVEY.values)[0]])+' that you just loaded!') #printing the SURVEY ID. 
    except TypeError:
        if SURVEY == 'LOWZ':
            print('LOWZ')
        else:
            print("Woops! This doesn't correspond to a valid IDSURVEY! Are you loading files correctly?")
            quit()

 
    if Param == 'c':
        matrixdic[s+Param] = MI.Matrix_c_init(dfpre, dfpost, cbinsize)
    else:
        matrixdic[s+Param] = MI.Matrix_x_init(dfpre, dfpost, xbinsize)
    if count == 0:
        DATOT = dfdata[['CID','x1', 'c', 'HOST_LOGMASS', 'zHD']] #The important parameters.
    else:
        DATOT = DATOT.append(dfdata[['CID','x1', 'c', 'HOST_LOGMASS', 'zHD']])
    lensdic[s+Param] = len(dfdata.loc[dfdata.HOST_LOGMASS > 3])
    count += 1 #I think this is fairly obvious?

q = 0
for s in SURVEY:
    q += lensdic[s+Param]

for num,s in enumerate(SURVEY):
    if num == 0:
        newmatrix = matrixdic[s+Param]*(lensdic[s+Param] / q) #Here we're making sure each survey gives a proportional contribution
    else:
        newmatrix += matrixdic[s+Param]*(lensdic[s+Param] / q)

print(len(DATOT))

import Optimiser #TODO: Have loading_files function to be main, then, grab arange from optimizer_mass and print the outputs to a file, instead of stdout
Optimiser.Optimiser(Param,DATOT, 4, dfpre, dfpost, .2, None)


