# VectorDeviation
Readme

The Vector Deviation functions quantify the differences between two sets of vectors in terms of magnitude and direction. The program’s functions expect input as a table of magnitude and direction of 2 vector datasets. It also accommodates input in the following formats, which is converted to a table of values to be used in the deviation functions:

a.	4 rasters of magnitude and direction of 2 vector datasets

b.	4 rasters of meteorological u and v components of 2 wind data sets

c.	Table of meteorological u and v components of 2 wind data sets

The functions compute the Mean Deviation (MD) and the Mean Absolute Deviation (MAD) for each property. The MAD is further broken down into the Quantity Component (QC) and the Allocation Component (AC), which are displayed in a stacked bar chart. The calculations of these values are based on equations from Metrics that Make a Difference by Dr. Robert Gilmore Pontius.


Functions

Functions for calculating vector deviations between two vector datasets:

•	magnitude_md() - calculates mean deviation of magnitude

•	magnitude_qc() - calculates quantity component of magnitude

•	magnitude_ac() - calculates allocation component of magnitude

•	magnitude_mad() - calculates mean absolute deviation of magnitude

•	direction_md() - calculates mean deviation of direction

•	direction_qc() - calculates quantity component of direction

•	direction_ac() - calculates allocation component of direction

•	direction_mad() - calculates mean absolute deviation of direction

•	mad_plot() - creates a double stacked bar plot of magnitude and direction mean absolute deviations with components

Other included functions:

•	raster_input_vector() - converts 4 rasters of magnitude and direction of 2 vector data sets into a table of the values to be used in deviation functions

•	raster_input_uv() - converts 4 rasters of meteorological u and v component of 2 wind data sets into a table of vector component values to be used in deviation functions

•	input_uv() - converts a table of meteorological u and v component of 2 wind data sets into a table of vector component values to be used in deviation functions

