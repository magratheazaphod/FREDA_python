#Rewriting some of my old MATLAB code for the RDA project into Python form,
#since MATLAB code requires continuing license to run...

#This particular code uses the output of the RDA algorithm and finds the mean of 
#rainband statistics such as frequency, latitude and intensity for time periods of choice.

from bootstrap import bs_means_diff, bs_stdofmean
import datetime
from itertools import compress
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import os
import pandas as pd
import scipy.stats
import time

#Path to RDA results
RDA_path_1 = "/Users/Siwen/Desktop/ferret/bin/meiyu_clean.nc"
RDA_path_2 = "/Users/Siwen/Desktop/ferret/bin/meiyu_2_clean.nc"


#1-liner for calculating standard deviation of Bernoulli process
def std_bernoulli(n,p):
    return (p*(1-p)/n)**.5

## dives into the netCDF file and returns all data matching date criteria
## note: freq_out returns a tuple [number_of_days, band_pct].
def collect_data(years, period, primary_only = False, latrange = [20,40]):
    RDA_1 = nc.Dataset(RDA_path_1, 'r') #all primary events
    RDA_2 = nc.Dataset(RDA_path_2, 'r') #all secondary events

    #load data from NetCDF files to notebook
    lat_1 = RDA_1.variables['lat_115'][:]
    lat_2 = RDA_2.variables['lat_115'][:]
    
    int_1 = RDA_1.variables['intensity'][:] 
    int_2 = RDA_2.variables['intensity'][:]
        
    countit_1 = RDA_1.variables['countit_1'][:]
    
    #Assign a calendar date to each time point
    startday = datetime.datetime(1951,1,1)
    date_list = [datetime.timedelta(days=x) + startday for x in range(0, 20819)]
    
    #the .timetuple().tm_yday command turns a datetime object into a day of the year.
    filter_days = [(dd.timetuple().tm_yday >= period[0]) & (dd.timetuple().tm_yday <= period[1]) & (dd.year >= years[0]) & (dd.year <= years[1]) for dd in date_list]
    days_list = list(compress(range(20819), filter_days))
    
    #filter by latitude. also effectively filters out NaN
    lat_1_pass = (lat_1 > latrange[0]) & (lat_1 < latrange[1])
    lat_1_days = list(compress(range(20819), filter_days&lat_1_pass))
    lat_2_pass = (lat_2 > latrange[0]) & (lat_2 < latrange[1])
    lat_2_days = list(compress(range(20819), filter_days&lat_2_pass))

    ##apply filters
    lat_1_out = lat_1[lat_1_days]
    lat_2_out = lat_2[lat_2_days]
    int_1_out = int_1[lat_1_days]
    int_2_out = int_2[lat_2_days]
   
    freq_out = {'n':sum(filter_days),'p':len(lat_1_out)/sum(filter_days)}
    
    #default is primary and secondary
    if primary_only == False:
        lats_out = np.append(lat_1_out, lat_2_out)
        ints_out = np.append(int_1_out, int_2_out)
    else:
        lats_out = lat_1_out
        ints_out = int_1_out
    
    RDA_1.close()
    RDA_2.close()
    
    return [freq_out, lats_out, ints_out, ]

## expected input is a 2-item list of tuples. Each tuple should contain
## a beginning and end year.
## likewise, period just be a single tuple. an external script runs through all the sets of years.

## calculations of standard deviation of mean and p-value are analytic for frequency
## for intensity and latitude, we use bootstrapping.
def compare_periods(years, period, primary_only = False, latrange = [0,99], tau=1):
    metrics = ['frequency','latitude','intensity']
    stats = ['mean_p1','mean_p2','std_p1','std_p2','diff_p2p1','pval']
    results = pd.DataFrame(columns=metrics,index=stats)
    
    [freq_p1, lats_p1, ints_p1] = collect_data(years[0], period,\
                                               primary_only=primary_only, latrange = latrange)
    [freq_p2, lats_p2, ints_p2] = collect_data(years[1], period,\
                                               primary_only=primary_only, latrange = latrange)
    data = {'latitude':(lats_p1,lats_p2),'intensity':(ints_p1,ints_p2)}
    
    #default number of iterations for bootstrap (min 2000 recommended)
    niter=10000
    
    #note that the autocorrelation time scale Tau reduces the number of observations
    #to N=n/Tau
    results['frequency']['mean_p1'] = freq_p1['p']
    results['frequency']['mean_p2'] = freq_p2['p']
    results['frequency']['std_p1'] = std_bernoulli(freq_p1['n']/tau,freq_p1['p'])
    results['frequency']['std_p2'] = std_bernoulli(freq_p2['n']/tau,freq_p2['p'])
    results['frequency']['diff_p2p1'] = freq_p2['p']-freq_p1['p']
    Z = results['frequency']['diff_p2p1']/\
    (results['frequency']['std_p1']**2+results['frequency']['std_p2']**2)**.5
    results['frequency']['pval'] = scipy.special.ndtr(Z)

    
    ## handles latitude and intensity
    for var,values in data.items():

        results[var]['mean_p1'] = data[var][0].mean()
        results[var]['mean_p2'] = data[var][1].mean()
        results[var]['std_p1'] = bs_stdofmean(data[var][0],niter)[1]
        results[var]['std_p2'] = bs_stdofmean(data[var][1],niter)[1]
        results[var]['diff_p2p1'] = results[var]['mean_p2']-results[var]['mean_p1']
        results[var]['pval'] = bs_means_diff(data[var][1],data[var][0],\
                                             niter, method='perm')[1]
        
    data['frequency']=(freq_p1,freq_p2)
    
    return results