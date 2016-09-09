#basic bootstrap functions for general use with gridded spatio-temporal data. Assumption is that data is in form long-lat-time,
#with bootstrapping generally performed along time axis.

#elementary bootstrap - given data, returns standard deviation of mean by bootstrapping (as well as mean)
import numpy as np
import numpy.random as rnd
import time


#Bootstrap resampling of sample data
#V is the sample space, nsamp is the number of samples that we're drawing from V
def bs_resample(V,nsamp):
    
    nn = len(V)
    indices = np.floor(nn * rnd.random_sample((nsamp,))).astype(int)
    Vnew = np.array(V[indices])
    return Vnew

#a resampling method that takes into account temporal autocorrelation by drawing data in consecutive blocks:
#the "Moving Blocks Bootstrap" - Künsch (1988)
def bs_resample_block(V,nsamp,blklen):
    
    nn = len(V)
    nblks = np.ceil(nsamp/blklen)
    
    #the number of possible different blocks is nn-blklen+1
    indices = np.floor((nn-blklen+1) * rnd.random_sample((nblks,))).astype(int)    
    Vnew = np.zeros(nsamp)
    
    for i in np.arange(nblks):
        
        #final block may not be of full length - must account for it
        myblklen = np.minimum(blklen,nsamp-blklen*i)
        Vnew[blklen*i : (blklen*i+myblklen)] = V[indices[i] : indices[i]+myblklen]

    return Vnew

#using the preceding bs_resample function, finds std dev of mean by bootstrap redrawing
#returns actual mean and its standard deviation based on bootstrapping
def bs_stdofmean(V, niter):
    
    nn = len(V)
    means = np.zeros((niter,))
    
    for i in np.arange(niter):
        means[i] = np.mean(bs_resample(V,nn))
        
    means_std = np.std(means)
    return np.mean(V), means_std

#equivalent to bs_stdofmean, but again, with samples drawn in consecutive blocks of length blklen
def bs_stdofmean_block(V, niter, blklen):
    
    nn = len(V)
    means = np.zeros((niter,))
    
    for i in np.arange(niter):
        means[i] = np.mean(bs_resample_block(V, nn, blklen))
        
    means_std = np.std(means)
    return np.mean(V), means_std


#code assesses the statistical significance of the difference in means between two data sets.
#two different (largely equivalent) techniques are implemented: bootstrapping with replacement and a permutation method (bootstrapping without replacement)
#Default: Bootstrap with replacement, can be specified as keyword.

def bs_means_diff(V1,V2,niter,method='bootstrap'):

    n1 = len(V1)
    n2 = len(V2)
    V = np.append(V1,V2)
    nn = len(V)
    
    return V
    
def bs_diff_blocks(V1,V2,blklen,niter):

    
    return