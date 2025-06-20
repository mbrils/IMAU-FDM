# IMAU-FDM v1.2G
This git repository contains the Fortran code of version 1.2G of the IMAU Firn Densification Model (IMAU-FDM). The development of IMAU-FDM has moved to the following repository: https://github.com/IMAU-ice-and-climate/IMAU-FDM. For the latest version of the code, please look there.

## Code description

* mainprogram_firnmodel.f90: loads the user supplied input settings (resolution, output frequency, implicit or explicit etc.), loads RACMO forcing and calls subprogram
* subprogram.f90: calls spin-up, contains time loop and calls physics functions and output functions.
* timeloop.f90: contains function definitions that are called during the spin-up and time loop concerning the model physics
* ini_model.f90: contains function definitions for loading input settings and initialising firn and snow properties
* openNetCDF.f90: contains function definitions for loading RACMO forcing
* output.f90: contains function definitions for outputting results (height integrated properties and depth profiles over time)

## Usage

Compile the code using the provided Makefile. In order to run the program the following data needs to be provided:

* NetCDF files containing the mean values of the input variable during the spin-up period. One file for each variable.
* A single NetCDF file containing the time series of all input variables during the simulation period, this file is created by cutting the (RACMO) forcing into strips of 6 or 7 grid cells wide and then concatenating the files containing different years.
* text file containing user provided settings (see ini_model.f90)
* text file containing dimensions of input forcing netcdf files

The model requires the following variables to be provided as input:
* Total precipitation
* Rain
* Snow drift
* Melt
* Skin temperature
* Sublimation
* absolute 10 m wind speed

Add the path to the directories containing the input to ini_model.f90 and openNetCDF.f90 then compile using the provided Makefile and run the executable.
