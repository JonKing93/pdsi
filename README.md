# pdsi
 A Matlab function to calculate Palmer Drought Severity Index (PDSI) using 
monthly temperature and precipitation data. Can determine values for 
multiple sites simultaneously.

### Setup
Matlab 2020b or higher is recommended, however this function should work on
versions as low as 2016b.

To use, download this repository and add it to your MATLAB path.

### Syntax

Please see the demo for detailed instructions on using this function. Basic
syntax is:
 
[X, Xm, Z, PE] = pdsi(T, P, years, lats, awcs, awcu, cafecYears, timeDim, progressbar)

---- Inputs -----

T: A numeric array holding monthly temperature data at various data sites. 
   Data should begin in a January, and end in a December. Should be in
   units of Celsius.

P: Monthly precipitation data. Should correspond same data sites and time
   steps as the temperature data. Should be in units of inches/month.

years: A two element vector indicating the first and last year of the
   temperature data.

lats: The latitudes of the data sites in degrees. A numeric array. Must be
   the same size as the temperature data, except for the monthly time
   dimension, which should have a length of 1.

awcs: Available water capacity of the surface layer for each site. A
   numeric array the same size as lats. A common default value is 1.

awcu: Available water capacity of the underlying layer for each site. A 
   numeric array the same size as lats. A common default value is 5.

cafecYears: A two element vector indicating the first and last years of the
   period to use for CAFEC normalizations.

timeDim: A positive integer indicating which dimension of the temperature
   and precipitation data is the monthly time step. If not provided, pdsi
   assumes time is along the first dimension.

progressbar: An optional input that can be used to display a progress bar
   for long computations. A scalar logical. Set it to true to display progress.

----- Outputs -----

X: Palmer Drought Severity Index at each site in each month.

Xm: Modified PDSI at each site in each month

Z: Z indices at each site in each month.

PE: Computed potential evapotranspiration at each site in each month.

### Tips for slow computations

It is usually fastest to run pdsi on all sites at once. However, the pdsi
function scales poorly for very large datasets that require MATLAB to swap
data to disk.

If you are experiencing slow run times, try running pdsi on chunks of
data sites. Running on 1000 to 10000 sites at a time could be a good
starting point.

### Credits

This PDSI function is based on a script written by Dave Meko. It was
updated by Jonathan King in 2020 to efficiently process multiple sites.