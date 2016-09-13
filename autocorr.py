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

    tau = np.sum(yy_ac*yy_ac)   
        
    return yy_ac, tau