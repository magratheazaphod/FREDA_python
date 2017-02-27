#Rewriting some of my old MATLAB code for the RDA project into Python form,
#since MATLAB code requires continuing license to run...

#This particular code uses the output of the RDA algorithm and finds the mean of 
#rainband statistics such as frequency, latitude and intensity for time periods of choice.

import datetime
from itertools import compress
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import os
import time

#Path to RDA results
RDA_path_1 = "/Users/Siwen/Desktop/ferret/bin/meiyu_clean.nc"
RDA_path_2 = "/Users/Siwen/Desktop/ferret/bin/meiyu_2_clean.nc"

## dives into the netCDF file and spits out all of the data of particular attribu
## dives into the netCDF file and spits out all of the data of particular attribu
def collect_data(years, period):
    
    RDA_1 = nc.Dataset(RDA_path_1, 'r') #all primary events
    RDA_2 = nc.Dataset(RDA_path_2, 'r') #all secondary events

    #load data from NetCDF files to notebook
    lat_1 = RDA_1.variables['lat_115'][:]
    lat_2 = RDA_2.variables['lat_115'][:]

    int_1 = RDA_1.variables['intensity'][:] 
    int_2 = RDA_2.variables['intensity'][:]
    
    #countit_1_all = 0
    #countit_2_all = 0
    
    #Assign a calendar date to each time point
    startday = datetime.datetime(1951,1,1)
    date_list = [datetime.timedelta(days=x) + startday for x in range(0, 20819)]
    
    #the .timetuple().tm_yday command turns a datetime object into a day of the year.
    filter_days = [(dd.timetuple().tm_yday >= period[0]) & \
                   (dd.timetuple().tm_yday <= period[1]) & \
                   (dd.year >= years[0]) & (dd.year <= years[1]) \
                   for dd in date_list]
    days = list(compress(range(20819),filter_days))
    
    #filter full list of latitudes and intensities to only include time period of interest.
    lat_1_out = lat_1[days]
    lat_2_out = lat_2[days]
    int_1_out = int_1[days]
    int_2_out = int_2[days]
    
    lats_out = np.append(lat_1_out, lat_2_out)
    ints_out = np.append(int_1_out, int_2_out)  
    lats_out = lats_out[~np.isnan(lats_out)]
    ints_out = ints_out[~np.isnan(ints_out)] 

    RDA_1.close()
    RDA_2.close()
    
    return [lats_out, ints_out]

## years below should be a tuple with a beginning year and end year.

## this time, the expected input is a 2-item list of tuples. Each tuple should contain
## a beginning and end year.
def compare_periods(years, periods, tau=1):
    return 2