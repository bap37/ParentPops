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
                cI_m = newmatrix
                nwalkers = 9
                ndim = 4
                p0 = np.random.rand(nwalkers, ndim)
                p0 = p0/1000
                p0 = np.abs(p0)
                sampler = emcee.EnsembleSampler(nwalkers, ndim, Matrix_c, args=[dfk, cI_m, flag])
                state = sampler.run_mcmc(p0, 2000, progress=True)
                samples = sampler.get_chain()
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