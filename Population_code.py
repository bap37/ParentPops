import pandas as pd
import numpy as np
import sys
sys.path.insert(0, '/home/bap37/Documents/Cosmology/SDSS/SDSS_allCandidates+BOSS/Functions') 
import BF #I will need to include this
import argparse
parser=argparse.ArgumentParser()
import emcee
from scipy.stats import binned_statistic

#Commented options don't work right now.

#parser.add_argument('--SURVEY', help='Choose the Survey') It would be nice to get this to work. To switch from the cumulative to a specific survey on command.
parser.add_argument('--Param', help='x1 or c')
parser.add_argument('--I_c', help='Initial guess for c', default=[0,.1,.1,2], nargs = '+')
parser.add_argument('--I_x1', help='Initial guess for x1', default=[0,1,1,2], nargs = '+')
#parser.add_argument('--Choice', help='Which bin you wish to investigate') LOW priority, but would be nice to select a specific range to investigate
#parser.add_argument('--Redshift', help='If you want to investigate redshift bins rather than mass', default=False) MEDIUM priority. Not the point of the paper, but super easy to implement. Just change masses to redshifts.
#parser.add_argument('--Test', help='Tests an iteration, more print statements', default=False) Junk.
#parser.add_argument('--Total', help='Cycle through all the surveys', default=False) MEDIUM Priority. Right now, it only does one survey per run. A flag to cycle through all surveys at once would be useful, but not necessary.
parser.add_argument('--C11', help='Toggle to True if using the C11 scatter model', default=False) #This does not work at present. It would be part of the loading data process down below. 
#parser.add_argument('--Iteration', help='Are you using a Gaussian distribution? If so, this should be true.', default=False) #junk.
#parser.add_argument('--Spec', help='Is your sample spectroscopic? default= No.', default=False) HIGH Priority. Photometry is the future of SNIa cosmology
parser.add_argument('--AB', help='Ugh', default=False) #Junk. We fixed that problem.
#A choice of distribution shapes is desired but difficult. Can discuss more. 
args = parser.parse_args()
Param = args.Param
I_c = args.I_c
I_c = [float(i) for i in I_c]
I_x1 = args.I_x1
I_x1 = [float(i) for i in I_x1]
#choice = args.choice
#redshift = args.redshift
Chotard= args.C11
#SPEC = args.Spec
AB = args.AB



SPEC = True
AS = 'AS'
import distutils.util
Chotard = distutils.util.strtobool(Chotard)
#SPEC = distutils.util.strtobool(SPEC)
AB = distutils.util.strtobool(AB)


#This function creates the colour (c) probability matrix. It's invoked once at the start of the "Optimiser_Mass" process
def Matrix_c_init():
    #cI_mean, cI_l, cI_r = params                                                                                                    
    send = {}
    length = {}
    bins = np.arange(-.4,.41,.02) #The range of colour values we consider.
    for i in bins:
        if i < .4:
            #0.4 is the upper bound for permissible colour values. 
            located = dfpre.loc[(dfpre.S2c > i ) & (dfpre.S2c <= i+0.02)] #Look in bins of colour ranging from i to i+0.02
            #print(located)                                                                                                          
            passed = dfpost.loc[dfpost.CID.isin(np.unique(located.CID.values))] #Not all supernova in "pre" are detected. This checks the unique ID, the CID, for those that are.
            dfo = np.histogram(passed.c.values, bins=bins)[0] #Convert values to a histogram. Establishing a 1D probability map.
            send[str(np.round(i+.005,3))] = np.array(dfo) #Starting to build a dictionary for each c input bin.
            length[str(np.round(i+.005,3))] = len(located) #Some colours are more likely than others. Here, we're going to use this to normalise each row of the matrix to make sure we aren't weighting anything improperly.

    for num,i in enumerate(send):
        if num == 0: #If this is the first iteration, we need to establish a normal row before doing any operations. 
            if np.sum(length[i]) != 0:
                 combined = (send[i])/length[i] #it's entirely possible the row is empty. Don't want to divide by 0. 
            else:
                 combined = (send[i])

        else:
            if np.sum(length[i]) != 0:
                 send[i] = (send[i])/length[i]
            combined = np.vstack((combined, send[i])) #It stacks rows. 
    cI_m = np.copy(combined)
    return cI_m #and the end result - an NxN probability matrix that we'll use to process our data. 

def Matrix_x_init(): #same as the colour probability matrix, but for stretch (x1).
    send = {}
    length = {}
    bins = np.arange(-3,3.1,.2)
    for i in bins:
        if i < 3:
            located = dfpre.loc[(dfpre.S2x1 > i ) & (dfpre.S2x1 <= i+0.2)]
            #print(located)                                                                                                         \
                                                                                                                                     
            passed = dfpost.loc[dfpost.CID.isin(np.unique(located.CID.values))]
            dfo = np.histogram(passed.x1.values, bins=bins)[0]
            send[str(np.round(i+.05,3))] = np.array(dfo)
            length[str(np.round(i+.05,3))] = len(located)

    for num,i in enumerate(send):
        if num == 0:
            if np.sum(length[i]) != 0:
                 combined = (send[i])/length[i]
            else:
                 combined = (send[i])

        else:
            if np.sum(length[i]) != 0:
                 send[i] = (send[i])/length[i]
            combined = np.vstack((combined, send[i]))
    xI_m = np.copy(combined)
    return xI_m

slist = ['DES', 'SDSS', 'FOUND', 'SNLS', 'PS1'] #The four high-redshift surveys. 
matrixdic = {} #store the probability matrix for each survey

lensdic = {} #Store relative contributions of each survey. SDSS has twice as many SNIa as DES, for instance. Need to account for that. 

#Param = 'c'
count = 0
for s in slist:
    #Everything is in FITOPTS. From there the structure is SURVEY_TYPE_MODEL.
    #SURVEY is what's in slist, line 106
    #TYPE is SPEC or PHOT, depending. Not all surveys have both, so some error trapping is necessary.
    #MODEL is G10 or C11.
    #Right now this is hardcoded for G10, which is unlabeled, but something that will be fixed when this moves off of my computer and into a more public space. 
    if s == 'DES':
        filename1 = '/home/bap37/Documents/Cosmology/MASS/FITOPTS/DES_SPEC/MATRIX_PRE.FITRES'
        filename2 = '/home/bap37/Documents/Cosmology/MASS/FITOPTS/DES_SPEC/MASS_0_DATA.FITRES'
        filename3 = '/home/bap37/Documents/Cosmology/MASS/FITOPTS/DES_SPEC/MATRIX_POST.FITRES'
        Names1, StartRow1 = BF.NAndR(filename1)
        Names2, StartRow2 = BF.NAndR(filename2)
        Names3, StartRow3 = BF.NAndR(filename3)

    elif s == 'SDSS':
        filename1 = '/home/bap37/Documents/Cosmology/MASS/FITOPTS/SDSS_SPEC/MATRIX_PRE.FITRES'
        filename2 = '/home/bap37/Documents/Cosmology/MASS/FITOPTS/SDSS_SPEC/MASS_0_DATA.FITRES'
        filename3 = '/home/bap37/Documents/Cosmology/MASS/FITOPTS/SDSS_SPEC/MATRIX_POST.FITRES'
        Names1, StartRow1 = BF.NAndR(filename1)
        Names2, StartRow2 = BF.NAndR(filename2)
        Names3, StartRow3 = BF.NAndR(filename3)
        
    elif s == 'PS1':
        filename1 = '/home/bap37/Documents/Cosmology/MASS/FITOPTS/PS1_SPEC/MATRIX_PRE.FITRES'
        filename2 = '/home/bap37/Documents/Cosmology/MASS/FITOPTS/PS1_SPEC/MASS_0_DATA.FITRES'
        filename3 = '/home/bap37/Documents/Cosmology/MASS/FITOPTS/PS1_SPEC/MATRIX_POST.FITRES'
        Names1, StartRow1 = BF.NAndR(filename1)
        Names2, StartRow2 = BF.NAndR(filename2)
        Names3, StartRow3 = BF.NAndR(filename3)
        

    else:
        filename1 = '/home/bap37/Documents/Cosmology/MASS/FITOPTS/'+s+'/MATRIX_PRE.FITRES'
        filename2 = '/home/bap37/Documents/Cosmology/MASS/FITOPTS/'+s+'/MASS_0_DATA.FITRES'
        filename3 = '/home/bap37/Documents/Cosmology/MASS/FITOPTS/'+s+'/MATRIX_POST.FITRES'
        Names1, StartRow1 = BF.NAndR(filename1)
        Names2, StartRow2 = BF.NAndR(filename2)
        Names3, StartRow3 = BF.NAndR(filename3)

    dfpre = pd.read_csv(filename1, header=None, skiprows=StartRow1,names=Names1, delim_whitespace=True)
    dfdata = pd.read_csv(filename2, header=None, skiprows=StartRow2,names=Names2, delim_whitespace=True)
    dfpost = pd.read_csv(filename3, header=None, skiprows=StartRow3,names=Names3, delim_whitespace=True)
    #Up until this comment, the above forloop is a fairly boilerplate process. For each survey, there are three "fitres" files. One for the actual data, one for the raw simulation output (pre), and one for simulations that have been processed as if they were data (POST)
    dfpre.CID = pd.to_numeric(dfpre.CID, errors='coerce') #sometimes the CIDs are stored as strings. Don't want that. 
    dfpre = dfpre.loc[dfpre.CID == dfpre.CID] #Getting rid of nans
    dfpre.CID = dfpre.CID.astype(int) #setting type. 
    dfpre.S2c = pd.to_numeric(dfpre.S2c, errors='coerce') #S2c is the "true" colour. The simulation chooses a "true" value at the start of the simulation. The process of observing the supernova, measurement noise, ambient noise, etc, can make the "true" colour not equal to the observed one. 
    dfpre.S2x1 = pd.to_numeric(dfpre.S2x1, errors='coerce') #Likewise, this is the "true" stretch.
    
    print(str(len(dfdata.loc[dfdata.HOST_LOGMASS > 3]))+' of '+str(len(dfdata))) #Just reading out the number of supernovae in each survey
    print(np.unique(dfdata.IDSURVEY.values)) #printing the SURVEY ID. 
    
    if Param == 'c':
        matrixdic[s+Param] = Matrix_c_init()
    else:
        matrixdic[s+Param] = Matrix_x_init()
    if count == 0:
        DATOT = dfdata[['CID','x1', 'c', 'HOST_LOGMASS']] #The important parameters.
    else:
        DATOT = DATOT.append(dfdata[['CID','x1', 'c', 'HOST_LOGMASS']])
    lensdic[s+Param] = len(dfdata.loc[dfdata.HOST_LOGMASS > 3])
    count += 1 #I think this is fairly obvious?

q = 0
for s in slist:
    q += lensdic[s+Param]

for num,s in enumerate(slist):
    if num == 0:
        newmatrix = matrixdic[s+Param]*(lensdic[s+Param] / q) #Here we're making sure each survey gives a proportional contribution
    else:
        newmatrix += matrixdic[s+Param]*(lensdic[s+Param] / q)

#newmatrix = 3*newmatrix
    
print(filename1)

def agauss(x,a,x0,sigma_l,sigma_r,n): #We need to propose an actual input distribution. Here I've chosen an asymmetric Gaussian. We can talk about whether or not this is the right choice, but it provides four different parameters. If that sounds overfitted, I have opinions!
    if x < x0:
        return a*np.exp(-np.abs(x-x0)**n/(n*sigma_l**n))
    else:
        return a*np.exp(-np.abs(x-x0)**n/(n*sigma_r**n))


def Matrix_c(params, dfk, cI_m, flag, debug=False):

    cI_mean, cI_l, cI_r, cI_n = params
    cI_n = 2
    #send = {}
    #length = {}
    bins = np.arange(-.4,.41,.02)
    
    
    cdatI = np.histogram(dfk.c.values, bins=bins)[0] #A histogram of the data

    cI = [1, cI_mean, cI_l, cI_r, cI_n]
    input_cI = []
    for x in ((bins[1:] + bins[:-1])/2):
        input_cI.append(agauss(x, *cI))
    input_cI = np.array(input_cI) #this creates a probability distribution for our estimate at the true population
    
    MP = np.matmul(input_cI, cI_m) #We take our true populaiton and process it through the probability matrix. This tells us what it should look like after being observed in a real telescope. Roughly.
    MP = MP*((np.sum(cdatI))/np.sum(MP)) #normalising it.
 
    Delta_c = cdatI - MP #This is detailed better in Scolnic and Kessler 2016, the paper that's attached in the email you got. But this is the difference between the data distribution and our guess.
    Delta_cerr = np.copy(cdatI) #we use Poisson statistics, so the error is the square root of the occupation of the bins. The data isn't normalised. 
    Delta_cerr[Delta_cerr == 0 ] = 1 #setting the unpopulated bins to 1 to avoid nans
    Delta_cerr = np.sqrt(np.abs(Delta_cerr)) #taking the square root.
    chi2 = []
    for i in range(len(Delta_c)):
        temp = ((Delta_c[i]*flag[i]))/((Delta_cerr[i])) #This is leftover from an error flagging protocol I had for a while. Certain bins were causing some huge statistcal biases, so this will ignore them in the chi2 calculation. flag always starts as 1.
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



def Matrix_x(params, dfk, xI_m, flag, debug=False): #same as Matrix_c, but for stretch.
    xI_mean, xI_l, xI_r, xI_n = params
    bins = np.arange(-3,3.1,.2)

    xdatI = np.histogram(dfk.x1.values, bins=bins)[0]
    
    
    xI = [1, xI_mean, xI_l, xI_r, xI_n]
    input_xI = []
    for x in ((bins[1:] + bins[:-1])/2):
        input_xI.append(agauss(x, *xI))
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
        temp = ((Delta_x[i]*flag[i]))/((Delta_xerr[i]))
        temp = temp**2
        chi2.append(temp)
    chi2 = np.array(chi2)
    LL = -np.sum(chi2)/2.
    if LL != LL:
        LL = -np.inf
    if (xI_l < 0) or (xI_r < 0):
        LL = -np.inf
    if (xI_l > 3) or (xI_r > 3):
        LL = -np.inf
    if (np.abs(xI_mean) > 3):
        LL = -np.inf
    if (np.abs(xI_n) > 4):
        LL = -np.inf
    if debug == True:
        return LL, (Delta_x/Delta_xerr)**2
    else:
        return LL

def Optimiser(Param): #better explained in Optimiser_Mass
    dfk = dfdata
    if Param == 'c':
        cI_m = Matrix_c_init()
    else:
        cI_m = Matrix_x_init()
    nwalkers = 7
    ndim = 3
    p0 = np.random.rand(nwalkers, ndim)
    p0 = p0/100
    p0 = np.abs(p0)
    if Param == 'c':
        sampler = emcee.EnsembleSampler(nwalkers, ndim, Matrix_c, args=[dfk, cI_m])
    else:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, Matrix_x, args=[dfk, cI_m])
    state = sampler.run_mcmc(p0, 1000, progress=True)
    sampler.reset()
    sampler.run_mcmc(state, 10000, progress=True)
    samples = sampler.get_chain()
    print(Param)
    print([np.mean(samples[:,:,0]), np.std(samples[:,:,0])])
    print([np.mean(samples[:,:,1]), np.std(samples[:,:,1])])
    print([np.mean(samples[:,:,2]), np.std(samples[:,:,2])])
    
#Now we need to cycle through masses
def Optimiser_Mass():
#The empty lists will be used to store the result for each parameter.
    stdls = []
    track = []
    stdrs = []
    means = []
    mass = []
    ns = []
    c = .6 #This is window size. Each mass bin will be from (mass - c) to (mass + c). 
    massstepper = np.arange(6,14,.2) #the masses. 
    for m in massstepper:
        fbf=False #This is vestigial. 
        dfk = DATOT.loc[(DATOT.HOST_LOGMASS <= m + c) & (DATOT.HOST_LOGMASS > m - c)] #The data, but cut down to the specified mass window. 
        print(len(dfk)) #How many supernova?
        if SPEC == True: #Spectroscopic data and photometric data is different. Spectroscopic data is more difficult to obtain, but is cleaner. So we can sacrifice some statistics for the spectroscopic dat.
            min_count = 50
        else:
            min_count = 150
        if (len(dfk) > min_count) & fbf==False:
            track.append(m)
            fbf = True
        if len(dfk) < min_count:
            pass #if the data doesn't have at least 50 events (or 150) for non-spectroscopic samples, we don't consider it. too noisy. 
        else:
            if Param == 'x1':
                flag = np.zeros_like(np.arange(-3,3.1,.2))[:-1]+1 #The flag array. Will disregard specific bins if they're big outliers.
                xI_m = newmatrix #The probability matrix
                nwalkers = 9 #This is specific to emcee (emcee.readthedocs.io), the python MC I'm using. But it sets up how many walkers in dimension space will be used.
                ndim = 4 #number of variables.
                p0 = np.random.rand(nwalkers, ndim)
                p0 = p0
                p0 = np.abs(p0) #Our initial guess for the variables.
                sampler = emcee.EnsembleSampler(nwalkers, ndim, Matrix_x, args=[dfk, xI_m, flag]) #Walkers, variables, function to run on. The args are what you would feed Matrix_x()
                state = sampler.run_mcmc(p0, 2000, progress=True) #Burn in. This would originally prime us to discard certain bins, then we would redo with the bins ignored. That's commented out now. 
                #Here we do a burn to find which bins to cut
                #samples = sampler.get_chain()
                #tempmean, tempstdl, tempstdr, tempn = np.mean(samples[:,:,0]), np.mean(samples[:,:,1]), np.mean(samples[:,:,2]), np.mean(samples[:,:,3])
                #chi2, OneD = Matrix_x([tempmean, tempstdl, tempstdr, tempn], dfk, xI_m,flag, debug=True)
                #flag[OneD > np.mean(OneD)*10] = 0
                #print(flag)
                #print('Number of discarded bins', str(len(flag[flag != 1])), str(np.arange(-3,3.1,.2)[:-1][flag == 0]), 'roughly')
                #print('originalish chi2' ,str(chi2))
                #
                sampler = emcee.EnsembleSampler(nwalkers, ndim, Matrix_x, args=[dfk, xI_m, flag]) #if the flag bit is commented out, this is a direct repeat
                state = sampler.run_mcmc(p0, 2000, progress=True)                
                sampler.reset() #We burned in our guess, now for the actual MCMC
                sampler.run_mcmc(state, 10000, progress=True)  #This walks around parameter space to minimise chi2.          
                samples = sampler.get_chain()
                tempmean, tempstdl, tempstdr, tempn = np.mean(samples[:,:,0]), np.mean(samples[:,:,1]), np.mean(samples[:,:,2]), np.mean(samples[:,:,3]) #our results.
                chi2, OneD = Matrix_x([tempmean, tempstdl, tempstdr, tempn], dfk, xI_m,flag, debug=True)   #I wanted to see the chi2. That's why debug is present.
                print(chi2) #Make sure it's not garbage.
                mass.append(m)
                means.append([np.mean(samples[:,:,0]), np.std(samples[:,:,0])])
                stdls.append([np.mean(samples[:,:,1]), np.std(samples[:,:,1])])
                stdrs.append([np.mean(samples[:,:,2]), np.std(samples[:,:,2])])
                ns.append([np.mean(samples[:,:,3]), np.std(samples[:,:,3])])
                
            elif Param == 'c': #same as x1. 
                flag = np.zeros_like(np.arange(-.4,.41,.02))[:-1]+1
                cI_m = newmatrix
                nwalkers = 9
                ndim = 4
                p0 = np.random.rand(nwalkers, ndim)
                p0 = p0/1000
                p0 = np.abs(p0)
                sampler = emcee.EnsembleSampler(nwalkers, ndim, Matrix_c, args=[dfk, cI_m, flag])
                state = sampler.run_mcmc(p0, 2000, progress=True)
                samples = sampler.get_chain()
                tempmean, tempstdl, tempstdr = np.mean(samples[:,:,0]), np.mean(samples[:,:,1]), np.mean(samples[:,:,2])
                chi2, OneD = Matrix_c([tempmean, tempstdl, tempstdr, 2], dfk, cI_m,flag, debug=True)
                #flag[OneD > np.mean(OneD)*6] = 0
                #print(OneD)
                print('Number of discarded bins', str(len(flag[flag != 1])), str(np.arange(-.4,.41,.02)[:-1][flag == 0]), 'roughly')
                #      
                try: #unfortunately, colour is a lot more likely to break than stretch. I'm not sure why, exactly. But the try/except thing works. Slows it down, but works. 
                    state = sampler.run_mcmc(p0, 2000, progress=True)
                except ValueError:
                    state = sampler.run_mcmc(p0, 2000, progress=True)
                print(np.mean(sampler.acceptance_fraction))
                sampler.reset()
                
                try:
                    sampler.run_mcmc(state, 10000, progress=True)
                except ValueError:
                    sampler.run_mcmc(state, 10000, progress=True)
                samples = sampler.get_chain()
                mass.append(m)
                means.append([np.mean(samples[:,:,0]), np.std(samples[:,:,0])])
                stdls.append([np.mean(samples[:,:,1]), np.std(samples[:,:,1])])
                stdrs.append([np.mean(samples[:,:,2]), np.std(samples[:,:,2])])
                ns.append(1)
    print(Param)
    if SPEC == False:
        #a horrible series of if statements designed to write out the results to a newly created file. The idea is that we have SURVEY_SCATTERMODEL_VALUE_PARAMETER.
        #SCATTER = G10 or C11
        #VALUE = stretch (x1) or colour (c)
        #PARAMETER = masses, means, stdls, stdrs, ns.
        #SURVEY = AS stands for "All surveys", which is the combination of DES, SNLS, PS1, and SDSS. 
                
        if Chotard == True:
            file = open('/home/bap37/Documents/Cosmology/MASS/PARENTS/'+AS+'_C11_'+Param+'_masses.txt', "w+")
        else:
            file = open('/home/bap37/Documents/Cosmology/MASS/PARENTS/'+AS+'_'+Param+'_masses.txt', "w+")
        for m in mass:
            file.write(str(m)+'\n')
        file.close()
    
        if Chotard == True:
            file = open('/home/bap37/Documents/Cosmology/MASS/PARENTS/'+AS+'_C11_'+Param+'_means.txt', "w+")
        else:
            file = open('/home/bap37/Documents/Cosmology/MASS/PARENTS/'+AS+'_'+Param+'_means.txt', "w+")
        for f in means:
            file.write(str(f)+'\n')
        file.close()

        if Chotard == True:
            file = open('/home/bap37/Documents/Cosmology/MASS/PARENTS/'+AS+'_C11_'+Param+'_stdls.txt', "w+")
        else:
            file = open('/home/bap37/Documents/Cosmology/MASS/PARENTS/'+AS+'_'+Param+'_stdls.txt', "w+")
        for f in stdls:
            file.write(str(f)+'\n')
        file.close()

        if Chotard == True:
            file = open('/home/bap37/Documents/Cosmology/MASS/PARENTS/'+AS+'_C11_'+Param+'_stdrs.txt', "w+")
        else:
            file = open('/home/bap37/Documents/Cosmology/MASS/PARENTS/'+AS+'_'+Param+'_stdrs.txt', "w+")
        for f in stdrs:
            file.write(str(f)+'\n')
        file.close()
    else:
        if Chotard == True:
            file = open('/home/bap37/Documents/Cosmology/MASS/PARENTS/'+AS+'_C11_'+Param+'_masses.txt', "w+")
        else:
            file = open('/home/bap37/Documents/Cosmology/MASS/PARENTS/'+AS+'_'+Param+'_masses.txt', "w+")
        for m in mass:
            file.write(str(m)+'\n')
        file.close()

        if Chotard == True:
            file = open('/home/bap37/Documents/Cosmology/MASS/PARENTS/'+AS+'_C11_'+Param+'_means.txt', "w+")
        else:
            file = open('/home/bap37/Documents/Cosmology/MASS/PARENTS/'+AS+'_'+Param+'_means.txt', "w+")
        for f in means:
            file.write(str(f)+'\n')
        file.close()

        if Chotard == True:
            file = open('/home/bap37/Documents/Cosmology/MASS/PARENTS/'+AS+'_C11_'+Param+'_stdls.txt', "w+")
        else:
            file = open('/home/bap37/Documents/Cosmology/MASS/PARENTS/'+AS+'_'+Param+'_stdls.txt', "w+")
        for f in stdls:
            file.write(str(f)+'\n')
        file.close()

        if Chotard == True:
            file = open('/home/bap37/Documents/Cosmology/MASS/PARENTS/'+AS+'_C11_'+Param+'_stdrs.txt', "w+")
        else:
            file = open('/home/bap37/Documents/Cosmology/MASS/PARENTS/'+AS+'_'+Param+'_stdrs.txt', "w+")
        for f in stdrs:
            file.write(str(f)+'\n')
        file.close()

        if Chotard == True:
            file = open('/home/bap37/Documents/Cosmology/MASS/PARENTS/'+AS+'_C11_'+Param+'_ns.txt', "w+")
        else:
            file = open('/home/bap37/Documents/Cosmology/MASS/PARENTS/'+AS+'_'+Param+'_ns.txt', "w+")
        for f in ns:
            file.write(str(f)+'\n')
        file.close()
    return 





if __name__ == "__main__": 
    Optimiser_Mass()

