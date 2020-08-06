#File contains the two matrix initialiser functions. 
import numpy as np
import pandas as pd

def Matrix_c_init(dfpre, dfpost, binsize):
    #cI_mean, cI_l, cI_r = params                                                                                                    
    send = {}
    length = {}
    bins = np.arange(-.4,.41,binsize) #The range of colour values we consider.
    for i in bins:
        if i < .4:
            #0.4 is the upper bound for permissible colour values. 
            located = dfpre.loc[(dfpre.S2c > i ) & (dfpre.S2c <= i+binsize)] #Look in bins of colour ranging from i to i+binsize
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

def Matrix_x_init(dfpre, dfpost, binsize): #same as the colour probability matrix, but for stretch (x1).
    send = {}
    length = {}
    bins = np.arange(-3,3.1,binsize)
    for i in bins:
        if i < 3:
            located = dfpre.loc[(dfpre.S2x1 > i ) & (dfpre.S2x1 <= i+binsize)]
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
