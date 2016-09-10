#basic bootstrap functions for general use with gridded spatio-temporal data. Assumption is that data is in form long-lat-time,
#with bootstrapping generally performed along time axis.

#elementary bootstrap - given data, returns standard deviation of mean by bootstrapping (as well as mean)
import matplotlib.pyplot as plt
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
#three different (largely equivalent) techniques are implemented: 
#bootstrapping with replacement, bootstrapping with replacement and mixing of samples, and a permutation method (bootstrapping without replacement)
#Default: Bootstrap with replacement without mixing, can be specified as keyword.
#Keywords: bs_nomix, bs_mix and perm
#Importantly, all the methods implemented here should lead to equivalent answers

def bs_means_diff(V1, V2, niter, method='bs_nomix', debug='n'):

    n1 = len(V1)
    n2 = len(V2)
    V = np.append(V1,V2)
    nn = len(V)
    diffs = np.zeros((niter,))
    
    #bs_mix - draw samples from joint pool of V1 and V2
    #p-value - percentage of samples whose difference is less than actual diff
    if method == 'bs_mix':
        
        for i in np.arange(niter):
            
            #debugging module to show result of each iteration.
            if debug == 'y':   
                
                s1 = bs_resample(V, n1)
                s2 = bs_resample(V, n2)
                print(s1)
                print(np.mean(s1))
                print(s2)
                print(np.mean(s2))
                diffs[i] = np.mean(s1)-np.mean(s2)
                print(diffs[i])
                time.sleep(3)
            
            else: #actual iteration is a 1-liner
                diffs[i] = np.mean(bs_resample(V, n1))-np.mean(bs_resample(V, n2))
                
        #calculate effective p-value
        actualdiff = np.mean(V1)-np.mean(V2)
        pval = (sum(actualdiff > diffs)+1)/(niter+1) 
        #normalized parameter to account for bias
            
    elif method == 'perm':
        
        for i in np.arange(niter):
            
            #debugging module to show result of each iteration.
            if debug == 'y':      
                
                #print(V1)
                #print(V2)
                #print(V)
                pm = rnd.permutation(nn)
                #print(pm)
                Vnew = np.array(V[pm])
                #print(Vnew)
                #time.sleep(10)
                s1 = Vnew[0:n1]
                s2 = Vnew[n1:nn]
                print(s1)
                print(s2)
                print(np.mean(s1))
                print(np.mean(s2))
                print(np.mean(s1)-np.mean(s2))
                time.sleep(2)
            
            else:
                
                Vnew = np.array(V[rnd.permutation(nn)])
                s1 = Vnew[0:n1]
                s2 = Vnew[n1:nn]
                
            diffs[i] = np.mean(s1)-np.mean(s2)
    
        #calculate effective p-value
        actualdiff = np.mean(V1)-np.mean(V2)
        pval = (sum(actualdiff > diffs)+1)/(niter+1) 
        #normalized parameter
            
    else: 
    #bs_nomix - no mixing between samples V1 and V2
    #p-value is the probability of obtaining a difference of 0    
        
        for i in np.arange(niter):
            
            #debugging module to show result of each iteration.
            if debug == 'y':      
                
                s1 = bs_resample(V1, n1)
                s2 = bs_resample(V2, n2)
                print(s1)
                print(np.mean(s1))
                print(s2)
                print(np.mean(s2))
                diffs[i] = np.mean(s1)-np.mean(s2)
                print(diffs[i])
                time.sleep(3)
            
            else:
                
                diffs[i] = np.mean(bs_resample(V1, n1))-np.mean(bs_resample(V2, n2))
    
        #calculate effective p-value
        actualdiff = np.mean(V1)-np.mean(V2)
        pval = (sum(diffs > 0)+1)/(niter+1) 
        #normalized parameter to account for bias
        
    
    #optional toggle produces a histogram to see the bootstrapped distribution     
    if debug == 'hist':
            
            plt.hist(diffs,20)
            plt.show()
            ax = plt.axes()
            xmin, xmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()
            print(ymin)
            print(.75*ymax)
            print(.7*ymax)
            
            #have to put the arrow in different spot depending on method
            if (method == 'bs_mix') | (method == 'perm'):
                xar=actualdiff
            else:
                xar=0
            
            ax.arrow(xar, .98*ymax, 0, -.05*ymax, head_width=.015*(xmax-xmin), head_length=.03*ymax)
            
            #draw lines to demarcate 95%/99% confidence intervals.
            diff_hist, bin_edges = np.histogram(diffs,500)
            cdf = np.cumsum(diff_hist)/niter
            l_99 = bin_edges[np.where(cdf > .005)[0][0]]
            l_95 = bin_edges[np.where(cdf > .025)[0][0]]
            u_95 = bin_edges[np.where(cdf > .975)[0][0]] 
            u_99 = bin_edges[np.where(cdf > .995)[0][0]] 
            
            #draw lines to show confidence intervals on histogram
            plt.plot([l_99,l_99], [ymin,ymax], 'r--', lw=2)
            plt.plot([l_95,l_95], [ymin,ymax], 'k--', lw=2)
            plt.plot([u_95,u_95], [ymin,ymax], 'k--', lw=2)
            plt.plot([u_99,u_99], [ymin,ymax], 'r--', lw=2)
            
    return actualdiff, pval

    
#same as the above, but with a moving blocks bootstrap - as in, data chosen in consecutive blocks. Only one method makes sence - bootstrap with replacement from within each sample.
def bs_diff_blocks(V1,V2,blklen,niter):

    
    return