# RDA_python
Python codes written to finish manuscript on the Rainband Detection Algorithm.

List of completed sub-projects:
-RDA_bars_with_whiskers.ipynb: Show the mean and std dev of band frequency and intensity during different time periods, and also overlays the p-value of significance of changes between time periods (51-79 v 80-07 and 80-93 v 94-07). Replaces table from thesis. Produces figure RDA_bar_final.pdf

-RDA_effective_rainfall.ipynb: Given rainband frequency and intensity, spits out the mean rainfall generated as RDA_effective_rainfall.nc. Decadal changes saved as RDA_effective_rainfall_diff.nc

-RDA_freq_diff and RDA_freq_diff_tau: Calculates the statistical significance of decadal changes in rainband frequency. The latter also accounts for the autocorrelation of rainbands (i.e. they tend to persist for multiple days, so effective sample size of observations is smaller).

-RDA_frequency_autocorrelation.ipynb - calculates the autocorrelation timescale tau used by RDA_freq_diff_tau.ipynb.

-RDA_hov_frequency_intensity_plots.ipynb: Creates Hovmoller plots of climatological rainband frequency, intensity and net rainfall produced, and also figures showing decadal changes. Produces figures hov_freq_int_climo.pdf, hov_freq_int_8007_5179.pdf and hov_freq_int_9407_8093.pdf.

-RDA_intensity.ipynb: Sweet script that figures out mean intensity of rainbands at a given latitude and day of the year by aggregating every event within an n-day and m-degree of latitude range.

Projects still in pipeline:
-Rerun significance of intensity changes with permutation method, instead of bootstrapping without mixing between samples.

-A companion to the 12-panel RDA_type_changes figure, except comparing 1994-2007 and 1980-1993.

-A figure showing the seasonal fraction of banded rainfall at different points (as opposed to just the yearly fraction)
 
-P-values of decadal changes in 'effective rainfall'
