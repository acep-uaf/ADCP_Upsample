# ADCP_Upsample

This R script takes reads the csv file of stream-wise velocity data and generates a new time series with comparable but distinct noise using the Langevin method and library. 
It simulates a time series measured at 1/10 the sampling rate by claculating a moving average with a window of 10, though other sampling rates can be calculated by adjusting the moving average window.

A stationary signal is calculated by taking the difference of the original signal and the moving average, then normalizing by the standard deviation of the moving average.
The stationary signal is used to calculate drift and diffusion coefficients.
Those coefficients are used to generate a new stationary signal.

The new stationary signal is un-normalized and added to the moving average (or new equivalent signal) to create a new up-sampled velocity time-series.
