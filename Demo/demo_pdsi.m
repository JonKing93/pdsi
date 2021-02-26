%% This is a demo for the pdsi function. It has two parts
%
% Demo 1 shows how to calculate pdsi for gridded climate model data in the
% Western US.
%
% Demo 2 shows how to calculate pdsi for weather statiion data in the
% Southern Hemisphere.

% Note that you can also access documentation for pdsi using
% >> help pdsi
% or
% >> doc pdsi
% from the console.
clear;
clc;

%% Demo 1: Calculate PDSI gridded climate model data in the Western US

% Load some model climate data. This demo data is temperature and precipitation
% over much of the Western US from 1900-2005 CE.
data = load('demo-data-Western-US.mat');
T = data.temperature;
P = data.precipitation;

% This climate data is a gridded longitude x latitude x time dataset, and
% we have some accompanying metadata.
lon = data.lon;
lat = data.lat;
time = data.time;

% To run the pdsi function, we must provide
% 1. Monthly Temperature in Celsius (beginning in January)
% 2. Monthly precipitation in mm/month (covering the same time period as temperature)
% 3. The first and last year of the data
% 4. The latitudes of the data sites
% 5. Available water capacity in the surface layer (in mm) at the sites, and
% 6. Available water capacity in the underground layer (in mm) at the sites.
% 7. The first and last years in which to apply CAFEC normalization.
% 8. Which dimension is the time dimension (If not the first dimension)
% 9. (Optional) Request a progress bar for lengthy computations.

% We'll start by converting the climate model data into the correct units.
% Currently temperature is in Kelvin and precipitation is in mm/second
disp(data.temperature_units);
disp(data.precipitation_units);

T = T - 273.15;   % From Kelvin to Celsius
P = P * 2.592E06;   % From mm/second to mm/month

% We have inputs 1 and 2, so let's work on the others. Looking at the time
% metadata, we can see the dataset spans the years 1900 to 2005. We decide
% to compute CAFEC normalizations over 1930 to 1970. We also have latitude
% metadata, which we can propagate over each longitude to get the latitude
% at every data site.
years = [1900 2005];
cafecYears = [1930 1970];
lats = repmat(lat, [numel(lon), 1]);

% Reasonable default values for available water capacities are 25.4 mm (in the
% surface layer), and 127 mm (in the underlying layer). We'll use those here.
awcs = 25.4 * ones(size(lats));
awcu = 127 * ones(size(lats));

% Our data is longitude x latitude x time. Since time is NOT the first
% dimension, we should note that time is along the third dimension using
% the eighth input.
timeDim = 3;

% We can now run the script
[X, Xm] = pdsi(T, P, years, lats, awcs, awcu, cafecYears, timeDim);
% Here, the outputs are X (the Palmer Drought Severity Index), and Xm
% (modified PDSI) at each of the data sites in each month.

% PDSI and modified PDSI are the most commonly used outputs, but we could
% also do:
[X, Xm, Z, PE] = pdsi(T, P, years, lats, awcs, awcu, cafecYears, timeDim);
% to also return Z indices (Z) and calculated potential 
% evapotranspiration (PE, in mm) at each site in each month.

% If the calculations are taking a while, we can use the ninth input to
% request a progress bar for the calculations.
progressbar = true;
X = pdsi(T, P, years, lats, awcs, awcu, cafecYears, timeDim, progressbar);


%% Demo 2: Calculate PDSI from weather station data in the Southern Hemisphere
clear;
clc;

% Now let's load some weather station data. This demo data has dimensions
% of (time x recording site). It includes some metadata for the time steps
% and latitude at each site.
data = load('demo-data-Southern-Hemisphere.mat');

T = data.temperature;
P = data.precipitation;

lat = data.lat;
time = data.time;
nTime = numel(time);

% Conveniently, the temperature data is in the correct units. However,
% we'll need to convert precipitation from inches to mm.
disp(data.temperature_units);
disp(data.precipitation_units);

P = P * 25.4;

% Now, the first month is July (instead of January), and the last month
% is June (instead of December). We'll need to remove these points that
% don't fall within a calendar year.
remove = [1:6, ...             % July through December of the first (incomplete) year
          nTime-5:nTime];      % January through June of the last (incomplete) year
T(remove,:) = [];
P(remove,:) = [];
time(remove) = [];

% Let's build the other inputs. Looking at the time metadata, we can see
% the first year is 1850, and the last year is 2004. We decide to do the
% CAFEC normalization from 1901 to 1950.
years = [1850 2004];
cafecYears = [1901 1950];

% We decide to use the default soil water capacities again.
awcs = 25.4 * ones(size(lat));
awcu = 127 * ones(size(lat));

% The time dimension is the first dimension for this dataset, so providing
% it is optional.
timeDim = 1;

% Now we can run the script
[X, Xm] = pdsi(T, P, years, lat, awcs, awcu, cafecYears, timeDim);
% The outputs are PDSI and modified PDSI at each station in each month.