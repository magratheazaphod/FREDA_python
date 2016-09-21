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
#debug parameter can be set to 'y', which gives verbose output, or also to 'hist'
#hist produces nice histogram of bootstrap simulation results and graphic view of 
#where p-value falls.
def bs_means_diff_block(V1,V2,niter,blklen,debug='n'):

    n1 = len(V1)
    n2 = len(V2)
    V = np.append(V1,V2)
    nn = len(V)
    diffs = np.zeros((niter,))
    
    for i in np.arange(niter):
            
        #debugging module to show result of each iteration.
        if debug == 'y':      
                
            s1 = bs_resample_block(V1, n1, blklen)
            s2 = bs_resample_block(V2, n2, blklen)
            print(s1)
            print(np.mean(s1))
            print(s2)
            print(np.mean(s2))
            diffs[i] = np.mean(s1)-np.mean(s2)
            print(diffs[i])
            time.sleep(3)
            
        else:
            
            #difference from bs_diff is a one-liner thanks to implementation
            diffs[i] = np.mean(bs_resample_block(V1, n1, blklen)) \
            -np.mean(bs_resample_block(V2, n2, blklen))
    
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
                        
            ax.arrow(0, .98*ymax, 0, -.05*ymax, head_width=.015*(xmax-xmin), head_length=.03*ymax)
            
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

    
    return

#same as bs_resample_block except that the data now comes in two dimensions, one of which is temporally correlated
#the other dimension is an ensemble dimension where members are NOT correlated (e.g. different years of a phenomenon).

#by convention, FIRST axis is temporally autocorrelated, and second is independent.

#unlike bs_resample_block, the second provided variable is now the TUPLE sampshape, which is the shape of the output (not necessarily same shape as input, but preserving same temporal dependence).
def bs_resample_block_ensemble_old(V,sampshape,blklen):
    
    #SIZE OF INPUT DATA
    Vlen = V.shape[0]
    Vmem = V.shape[1]
    
    #INTENDED SIZE OF OUTPUT DATA
    nn = sampshape[0]
    nblks = np.ceil(nn/blklen).astype(int)
    wdth = sampshape[1]
    
    #the number of possible different blocks is nn-blklen+1
    x_indices = np.floor((Vlen-blklen+1) * rnd.random_sample((nblks,wdth))).astype(int)    
    y_indices = np.floor(Vmem * rnd.random_sample((nblks,wdth))).astype(int)
    
    #print(x_indices)
    #print(y_indices)
    #time.sleep(10)
    
    Vnew = np.zeros(sampshape)
    
    for i in np.arange(nblks):
               
        for j in np.arange(wdth):
        
            #final block may not be of full length - must account for it
            myblklen = np.minimum(blklen,nn-blklen*i)
            Vnew[blklen*i : (blklen*i+myblklen), j] = V[x_indices[i,j] : x_indices[i,j]+myblklen, y_indices[i,j]]

    return Vnew

## SEPTEMBER 20th - REWRITING BS_RESAMPLE_BLOCK_ENSEMBLE to try to speed up bottlenecks.
#result according to timeit - improved performance by 1/3
def bs_resample_block_ensemble(V,sampshape,blklen):
    
    #SIZE OF INPUT DATA
    Vlen = V.shape[0]
    Vmem = V.shape[1]
    
    #INTENDED SIZE OF OUTPUT DATA
    nn = sampshape[0]
    nblks = np.ceil(nn/blklen).astype(int)
    wdth = sampshape[1]
    
    #the number of possible different blocks is nn-blklen+1
    x_indices = np.floor((Vlen-blklen+1) * rnd.random_sample((nblks,wdth))).astype(int)    
    y_indices = np.floor(Vmem * rnd.random_sample((nblks,wdth))).astype(int)
    
    #print(x_indices)
    #print(y_indices)
    #time.sleep(10)
    
    Vnew = np.zeros(sampshape)
    
    #CALCULATE length of last block - whole block may not fit in.
    lastblklen = nn % blklen 
    #print(lastblklen)
    #time.sleep(5)
    
    for j in np.arange(wdth):

        for i in np.arange(nblks-1):       
            x_index = x_indices[i,j]
            Vnew[blklen*i : blklen*(i+1), j] = V[x_index : x_index+blklen, y_indices[i,j]]
            
        #LAST BLOCK may be of different length, in which case we draw whole block but just put in whatever fits.
        x_index = x_indices[nblks-1,j]
        Vnew[nn-lastblklen : nn, j] = V[x_index : x_index+lastblklen, y_indices[nblks-1,j]]

    return Vnew


## SEPTEMBER 21st - continuing to try to speed up bs_resample
# New try - instead of having to look up blocks in 2-D array, create a dictionary of all blocks
#at the beginning. should make lookup much faster?
def bs_resample_block_ensemble_dict(V,sampshape,blklen):
    
    #SIZE OF INPUT DATA
    Vlen = V.shape[0]
    Vmem = V.shape[1]
    
    #INTENDED SIZE OF OUTPUT DATA
    nn = sampshape[0]
    nblks = np.ceil(nn/blklen).astype(int)
    total_rows = sampshape[1]
    shape_out = (nblks,total_rows)
    
    #MAKE DICTIONARY OF ALL POSSIBLE DOMINOS - more memory-intensive than previous method
    dominos=[]
    dominos_per_row = nn-blklen+1 #total number of dominos in each row
    total_dominos = dominos_per_row*total_rows
    
    for j in np.arange(total_rows):
        for i in np.arange(dominos_per_row):
            dominos.append(V[i : i+blklen,j])
    
    #DEFINE DICTIONARY
    domino_dict = { i:dominos[i] for i in range(total_dominos) }
    
    #CALCULATE length of last block - whole block may not fit in.
    lastblklen = nn % blklen 
    Vnew = np.zeros(sampshape)
    
    #fill in the output vector from original data, one domino at a time
    for j in np.arange(total_rows):

        for i in np.arange(nblks-1):       
            
            Vnew[blklen*i : blklen*(i+1), j] = domino_dict[rnd.randint(0,total_dominos)]
            
        #LAST BLOCK may be of different length - put in whatever fits.
        Vnew[nn-lastblklen : nn, j] = domino_dict[rnd.randint(0,total_dominos)][0:lastblklen]

    return Vnew



#same as bs_means_diff_block, with data drawn in consecutive blocks, except the data now come in multiple spatial dimensions.
#one dimension is an ensemble dimension (not correlated in time), and the other dimension is a time dimension (correlated in time).

#convention - FIRST dimension is correlated in time, second dimension cycles through ensemble members.

#debug parameter can be set to 'y' (verbose) or to 'hist' (histogram)

def bs_means_diff_block_ensemble(V1, V2, niter, blklen, debug='n'):

    len1 = V1.shape[0] #CORRELATED dimension in time
    len2 = V2.shape[0]
    
    mem1 = V1.shape[1] #UNCORRELATED ensemble dimension
    mem2 = V2.shape[1]
    
    diffs = np.zeros((niter,))
    
    for i in np.arange(niter):
            
        #debugging module to show result of each iteration.
        if debug == 'y':      
                
            s1 = bs_resample_block_ensemble(V1, (len1,mem1), blklen)
            s2 = bs_resample_block_ensemble(V2, (len2,mem2), blklen)
            print(s1)
            print(np.mean(s1))
            print(s2)
            print(np.mean(s2))
            diffs[i] = np.mean(s1)-np.mean(s2)
            print(diffs[i])
            time.sleep(3)
            
        else:
            
            #difference from bs_diff is a one-liner thanks to implementation
            diffs[i] = np.mean(bs_resample_block_ensemble(V1, (len1,mem1), blklen)) \
            -np.mean(bs_resample_block_ensemble(V2, (len2,mem2), blklen))
    
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
                        
            ax.arrow(0, .98*ymax, 0, -.05*ymax, head_width=.015*(xmax-xmin), head_length=.03*ymax)
            
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


#NEW FUNCTION: given daily rainfall at each latitude, the length of smoothing in days, number of iterations
#and block length blklen, returns a matrix of p-values corresponding to the SIGNIFICANCE OF CHANGE in that data set
#between the two time periods.
import numpy as np

#yrs1 and yrs2 should both be tuples with two elements - beginning year and ending year
#shape of P should be days (1st axis), latitude (2nd axis), ensemble members/years (3rd axis)

def bs_diff_rain(P,yrs1,yrs2,daysmth,niter,blklen):
    
    print(P.shape)
    pval = np.zeros((P.shape[0],P.shape[1]))
    ll = int((daysmth-1)/2) #15-day smoothing -> 7 days before and after
    
    P1 = P[:,:,yrs1[0]-1 : yrs1[1]] 
    P2 = P[:,:,yrs2[0]-1 : yrs2[1]]
    #print(P1.shape)
    #print(P2.shape)

    for lat in np.arange(P1.shape[1]):
        
        samp1 = P1[:,lat,:]
        samp2 = P2[:,lat,:]
        
        #following is important trick that I have to do to get day axis to wrap around...pad axis on each side
        samp1_pad = np.pad(samp1, [(ll, ll), (0, 0)], 'wrap')
        samp2_pad = np.pad(samp2, [(ll, ll), (0, 0)], 'wrap')
                
        for dd in np.arange(P1.shape[0]):
            
            start = time.time()

            s1 = samp1_pad[dd : dd+daysmth, :]
            s2 = samp2_pad[dd : dd+daysmth, :]
                            
            pval[dd,lat] = bs_means_diff_block_ensemble(s1, s2, niter, blklen)[1]
            end = time.time()
            print(end - start)


            
    return pval