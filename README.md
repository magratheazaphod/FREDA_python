# RDA_python
## Python codes written to finish manuscript on the Rainband Detection Algorithm.

### List of completed projects:

asia_rain_diff.ipynb: Produces figure SFND.eps showing changes in full year rainfall between time periods of interest.

autocorr.py: function to calculate the autocorrelation rho of a time series with itself at arbitrary time lags. also returns the autocorrelation time scale of the mean tau, as described in function header.

bootstrap.py: collection of bootstrap codes to create resamples using different methods, and also to calculate the significance of changes in mean between time periods. Also implements a moving blocks bootstrap (each sample includes multiple consecutive days) to account for the autocorrelation in some of the time series of interest in this project.

china_rain_diff.ipynb: Crucial workhorse notebook that calculates the significance of decadal changes in zonally averaged rainfall using a bootstrap with 2,000 iterations. Roughly 10 hours per run! Need to enter script and edit time period of interest and type of rainfall. Produces output in the form chinarain_diff_pval_notaiwan_{years2}_{years1}_{type}_{block_length}_{num_of_iterations} in NetCDF format.

china_rain_diff_fig.ipynb: Shows results from china_rain_diff in figure format. Produces figures chinarain_diff_notaiwan_8007_5179.pdf and chinarain_diff_notaiwan_9407_8093.pdf.

china_rain_diff_taiwan_test.ipynb: A quick proof that Taiwan is not the exclusive cause of apparent huge rainfall changes in southern China between 1994-2007 and 1980-1993. Nonetheless, went and repeated all of the china_rain_diff analysis in paper leaving out Taiwan.

RDA_bars_with_whiskers.ipynb: Show the mean and std dev of band frequency and intensity during different time periods, and also overlays the p-value of significance of changes between time periods (51-79 v 80-07 and 80-93 v 94-07). Replaces table from thesis. Produces figure RDA_bar_final.pdf

RDA_effective_rainfall.ipynb: Given rainband frequency and intensity, spits out the mean rainfall generated as RDA_effective_rainfall.nc. Decadal changes saved as RDA_effective_rainfall_diff.nc

RDA_freq_diff and RDA_freq_diff_tau: Calculates the statistical significance of decadal changes in rainband frequency. The latter also accounts for the autocorrelation of rainbands (i.e. they tend to persist for multiple days, so effective sample size of observations is smaller).

RDA_frequency_autocorrelation.ipynb - calculates the autocorrelation timescale tau used by RDA_freq_diff_tau.ipynb.

RDA_hov_frequency_intensity_plots.ipynb: Creates Hovmoller plots of climatological rainband frequency, intensity and net rainfall produced, and also figures showing decadal changes. Produces figures hov_freq_int_climo.pdf, hov_freq_int_8007_5179.pdf and hov_freq_int_9407_8093.pdf.

RDA_hovmoller_frequency.ipynb: Uses output of RDA algorithm to create a hovmoller plot of frequency for each latitude and day of year for an arbitrary set of years of interest. Output saved as RDA_freq_hov.nc.

RDA_hovmoller_intensity_construct.ipynb: Arranges all observations of rainband intensity in useful 3-dimensional format - arranged by day of year, latitude and a record direction in case there are multiple occurrences. Saves output as RDA_int_hov.nc

RDA_intensity_bootstrap.ipynb: bootstrapping script that figures out p-value associated with changes in intensity between decades. Value for each day of the year and latitude.

RDA_intensity_climo.ipynb: Sweet script that figures out mean intensity of rainbands at a given latitude and day of the year by aggregating every event within an n-day and m-degree of latitude range.

RDA_intensity_diff.ipynb: Basically an abandoned first try at finding the significance of changes in intensity between decades.

RDA_precip_china_type: Takes 57 years of China rainfall and finds latitude-mean rainfall for each day, latitude and year. Most recent version masks out Taiwan since physics of rainfall there distinct from what's going on on mainland. Saved as Pchina_type_notaiwan.nc

RDA_type_seasonal: Used to eventually produce a 12-panel figure of changes in different types of rainfall during particular stages of the year. Comparison is between two sets of years. Output is RDA_type_season_yrly.nc.

seaborn_sandbox.ipynb: Figuring out tsplot command from seaborn.




### Projects still in progress:

-Rerun significance of intensity changes with permutation method, instead of bootstrapping without mixing between samples.

-A companion to the 12-panel RDA_type_changes figure, except comparing 1994-2007 and 1980-1993.

-A figure showing the seasonal fraction of banded rainfall at different points (as opposed to just the yearly fraction)
 
-P-values of decadal changes in 'effective rainfall'

-Rewrite the original 12-panel rainfall change figure to use units of mm, not mm/day
