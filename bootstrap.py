#basic bootstrap functions for general use with gridded spatio-temporal data. Assumption is that data is in form long-lat-time,
#with bootstrapping generally performed along time axis.

#elementary bootstrap - given data, returns standard deviation of mean by bootstrapping (as well as mean)
import numpy as np
import numpy.random as rnd

#V is the sample space, nsamp is the number of samples that we're drawing
def bs_resample(V,nsamp):
    nn = len(V)
    indices = np.floor(nn * rnd.random_sample((nsamp,))).astype(int)
    Vnew = np.array(V[indices])
    return Vnew
    
#using the preceding bs_resample function, finds std dev of mean by bootstrap redrawing
#returns actual mean and its standard deviation based on bootstrapping
def bs_mean(V,niter, axis=2):
    nn = len(V)
    means = np.empty((niter,))
    
    for i in np.range(niter):
        means[i] = np.mean(bs_resample(V,nn)))
        
    means_std = np.std(means)
    return np.mean(V), means_std

def bs_mean_blocks(V,blklen,niter):
    
    return

def bs_diff(V1,V2,blklen,niter):

    return
    
def bs_diff_blocks(V1,V2,blklen,niter):

    return