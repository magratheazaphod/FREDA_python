#Rewriting some of my old MATLAB code for the RDA project into Python form,
#since MATLAB code requires continuing license to run...

#This particular code uses the output of the RDA algorithm and finds the mean of 
#rainband statistics such as frequency, latitude and intensity for time periods of choice.

import numpy as np
import time
import os
import netCDF4 as nc
import datetime
import matplotlib.pyplot as plt
%matplotlib notebook

#Path to RDA results
RDA_path_1 = "/Users/Siwen/Desktop/ferret/bin/meiyu_clean.nc"
RDA_path_2 = "/Users/Siwen/Desktop/ferret/bin/meiyu_2_clean.nc"

## dives into the netCDF file and spits out all of the data of particular attribu
def collect_data(years, period):
    
    RDA_1 = nc.Dataset(RDA_path_1, 'r') #all primary events
    RDA_2 = nc.Dataset(RDA_path_2, 'r') #all secondary events

    #load data from NetCDF files to notebook
    lat_1_all =  RDA_1.variables['lat_115'][:]
    lat_2_all =  RDA_2.variables['lat_115'][:]
    intensity_1_all = RDA_1.variables['intensity'][:] 
    intensity_2_all = RDA_2.variables['intensity'][:]
    #close out netcdf files
    
    
    RDA_1.close()
    RDA_2.close()

    
    return 3

## years below should be a tuple with a beginning year and end year.

## this time, the expected input is a 2-item list of tuples. Each tuple should contain
## a beginning and end year.
def compare_periods(years, periods, tau=1):
    return 2