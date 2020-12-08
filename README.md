# VectorDeviation
Functions for calculating vector deviation metrics in R.
All vecor deviation calculaations are based on the work of Robert G. Pontius. 


Functions for calculating vector deviations between two vector datasets:

magnitude_md() - calculates mean deviation of magnitude 
magnitude_qc() - calculates quantity component of magnitude
magnitude_ac() - calculates allocation component of magnitude
magnitude_mad() - calculates mean absolute deviation of magnitude

direction_md() - calculates mean deviation of direction
direction_qc() - alculates quantity component of direction
direction_ac() - calculates allocation component of direction
direction_mad() - calculates mean absolute deviation of direction

mad_plot() - creates a double stacked bar plot of magnitude and direction mean absolute deviations with components

Other included functions:

raster_input_vector() - converts 4 rasters of magnitude and direction of 2 vector data sets into a table of the values to be used in deviation functions
raster_input_uv() - converts 4 rasters of meteorological u and v component of 2 wind data sets into a table of vector component values to be used in deviation functions

input_uv() - converts a table of meteorological u and v component of 2 wind data sets into a table of vector component values to be used in deviation functions
