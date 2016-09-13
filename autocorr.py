#Homemade autocorrelation function, since evidently none exists as a built-in in Python.

#imported libraries
import matplotlib.pyplot as plt
import numpy as np
import time
import scipy as sp
from scipy.ndimage.interpolation import shift

#Bootstrap resampling of sample data
#V is the sample space, nsamp is the number of samples that we're drawing from V
def autocorr(yy,lags=10):
        
    nn = len(yy)

    if lags >= nn:
        lags = nn-1
            
    mn = np.mean(yy)
    stdv = np.std(yy)
    yy_norm = (yy - mn) / stdv
    yy_ac = np.zeros((lags+1,))
    
    for i in np.arange(lags+1):
        yy_ac[i] = (np.sum(yy_norm * shift(yy_norm, i, cval=0))) / (nn-i)

    #calculate Tau, the "decorrelation time scale." Must be careful because the
    #formula for tau depends on what we're calculating it for (mean, variance, etc.)
    
    #for the mean, it ends up being 1 + sum(2*rho) where rho is autocorrelation
    tau = 1 + 2*np.sum(yy_ac[1:])   
        
    return yy_ac, tau