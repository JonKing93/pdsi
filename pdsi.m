function[X, Xm, Z, PE] = pdsi(T, P, years, lats, awcs, awcu, calibrationYears, dim, showprogress)
%% Calculates PDSI
%
% [X, Xm] = pdsi(T, P, years, lats, awcs, awcu, calibrationYears)
% Calculates the Palmer Drought Severity Index and modified Palmer Drought
% Severity Index given monthly temperature and precipitation data from a
% collection of sites.
%
% [X, Xm, Z, PE] = pdsi(...)
% Also return the Z indices and potential evaportranspiration.
%
% [...] = pdsi(..., dim)
% Specify which dimension of the temperature and precipitation data is the
% monthly time step. By default, pdsi treats the first dimension as time.
%
% [...] = pdsi(..., dim, showprogress)
% Specify whether to show a progress bar. By default, does not display a
% progress bar.
%
% ----- Inputs -----
%
% T: Monthly temperature data (in Celsius) for a collection of sites. A
%    numeric array that cannot contain NaN, Inf, or complex values. Assumed
%    to progress from January to December in each included year.
%
% P: Monthly precipitation data (in inches) for a collection of sites. A 
%    numeric array that cannot contain NaN, Inf, or complex values. Must 
%    progress from January to December in each included year. Must have the
%    same size as the temperature data.
%
% years: The first and last year of the T and P data. A two element vector.
%
% lats: The latitudes of the sites. A numeric array. Must have the same size
%    as the temperature data, except for the monthly time dimension, which
%    must have a length of 1.
%
% awcs: The available water capacity of the surface layer for each site. A
%    numeric array. Must have the same size as the temperature data, except
%    for the monthly time dimension, which must have a length of 1.
%
% awcu: The available water capacity of the underlying layer for each site.
%    A numeric array. Must have the same size as the temperature data,
%    except for the monthly time dimension, which must have a length of 1.
%
% calibrationYears: The first and last year of the calibration period. A
%    two element vector.
%
% showprogress: A scalar logical indicating whether to display a progress
%    bar (true) or not (false).
%
% ----- Outputs -----
%
% X: Palmer Drought Severity Index for each site over time.
%
% Xm: Modified Palmer Drought Severity Index for each site over time.
%
% Z: The Z indices used to calculate PDSI for each site over time.
%
% PE: The computed potential evapotranspiration for each site over time.

% ----- Written By -----
% Original function by Dave Meko
% Updated and vectorized by Jonathan King, 2020

%% Setup

% Error check and determine progress bar
if ~exist('showprogress','var') || isempty(showprogress)
    showprogress = false;
end
assert(isscalar(showprogress) && islogical(showprogress), 'showprogress must be a scalar logical');

% Error check dim and determine time dimension
if ~exist('dim','var') || isempty(dim)
    dim = 1;
end
assert(isscalar(dim), 'dim must be scalar');
assert(mod(dim,1)==0 & dim>0, 'dim must be a positive integer');

% Error check numeric data type
args = {T, P, years, lats, awcs, awcu, calibrationYears, dim};
names = {'T','P','years','lats','awcs','awcu','calibrationYears','dim'};
for k = 1:numel(args)
    assert(isnumeric(args{k}) & ~any(isnan(args{k}),'all') & ~any(isinf(args{k}),'all') & all(isreal(args{k}),'all'), ...
        sprintf('%s must be numeric and cannot contain NaN, Inf, or complex values', names{k}));
end

% Error check sizes for data and years
fullSize = size(T);
assert( mod(fullSize(dim),12)==0, ['The monthly time dimension (dimension %.f), ',...
    'must have a length that is divisible by 12 (Each year must have 12 months)']);
assert(isequal(fullSize, size(P)), 'T and P must have the same size');
assert(isvector(years) & numel(years)==2, 'years must be a vector with two elements');
assert(isvector(calibrationYears) & numel(calibrationYears)==2, ...
    'calibrationYears must be a vector with two elements');

% Error check sizes for inputs that require a singleton time dimension
siz = fullSize;
siz(dim) = 1;
nDims = numel(siz);
for k = 4:6
    argSize = size(args{k});
    argSize(end+1:nDims) = 1;
    assert(isequal(siz, argSize), sprintf(['The size of %s must match the size ',...
        'of T, except for the monthly time dimension (dimension %.f), which ',...
        'must have a length of 1.'], names{k}, dim));
end

% Error check bounded quantities
assert(all(P>=0,'all'), 'P cannot have values less than 0');
assert(all(lats>=-90 & lats<=90, 'all'), 'Values in lats must be between -90 and 90 (inclusive).');
assert(all(awcs>=0, 'all'), 'awcs cannot have values less than 0');
assert(all(awcu>=0, 'all'), 'awcu cannot have values less than 0');

% Error check year args
assert(years(2)>=years(1), 'The second element of years cannot be smaller than the first element');
assert(calibrationYears(2)>=calibrationYears(1), 'The second element of calibrationYears cannot be smaller than the first element');
assert(all(mod(years,1)==0), 'The elements of years must be integers');
assert(all(mod(calibrationYears,1)==0), 'The elements of calibrationYears must be integers');

% Get the years and error check size
years = (years(1):years(2))';
nYears = numel(years);
assert(fullSize(dim)/12==nYears, sprintf(['The number of listed years (%.f), ',...
    'does not match the number of years in T (%.f)'], nYears, fullSize(dim)/12));

% Get the calibration period and error check
assert(all(ismember(calibrationYears, years)), 'The calibration years must be within the years of the data');
calib = years>=calibrationYears(1) & years<=calibrationYears(2);
nCalib = sum(calib);

% Permute inputs so that time is along the first dimension
order = 1:nDims;
order(1) = dim;
order(dim) = 1;

T = permute(T, order);
P = permute(P, order);
awcs = permute(awcs, order);
awcu = permute(awcu, order);

% Reshape N-dimensional arrays to matrix form
fullSize = fullSize(order);
fullMatrix = [fullSize(1), prod(fullSize(2:end))];
T = reshape(T, fullMatrix);
P = reshape(P, fullMatrix);

siz = siz(order);
sizMatrix = [siz(1), prod(siz(2:end))];
lats = reshape(lats, sizMatrix);
awcs = reshape(awcs, sizMatrix);
awcu = reshape(awcu, sizMatrix);

%% Calculations

% Separate months from years
[nTime, nSite] = size(T);
T = reshape(T, [12, nYears, nSite]);
P = reshape(P, [12, nYears, nSite]);

% Compute potential evapotranspiration via the Thornthwaite method
PE = pethorn(T, years, lats, calib);

% Get monthly means from the calibration period.
Pcalib = P(:,calib,:);
PEcalib = PE(:,calib,:);
Pmean = mean(Pcalib, 2);
PEmean = mean(PEcalib, 2);

% Initialize the soild moisture model by running it 10 years on the monthly
% means. Assume saturated starting conditions.
nInitial = 10;
Pi = repmat(squeeze(Pmean), [nInitial, 1]);
PEi = repmat(squeeze(PEmean), [nInitial, 1]);
s = soilMoisture(Pi, PEi, awcs, awcu, awcs, awcu, 1);

% Use the soil moisture in the last year of the initialization as the
% starting condition for additional runs of the soil moisture model
ssi = s.ssi(end-11,:);
sui = s.sui(end-11,:);

% Run the soil moisture model over the calibration period.
Pcalib = reshape(Pcalib, [12*nCalib, nSite]);
PEcalib = reshape(PEcalib, [12*nCalib, nSite]);
s = soilMoisture(Pcalib, PEcalib, awcs, awcu, ssi, sui, 2);

% Get the monthly means for recharge, runoff and loss from the soil
% moisture calibration period
et = mean(reshape(s.et, [12, nCalib, nSite]), 2);
r = mean(reshape(s.r, [12, nCalib, nSite]), 2);
pr = mean(reshape(s.pr, [12, nCalib, nSite]), 2);
ro = mean(reshape(s.ro, [12, nCalib, nSite]), 2);
pro = mean(reshape(s.pro, [12, nCalib, nSite]), 2);
loss = mean(reshape(s.loss, [12, nCalib, nSite]), 2);
ploss = mean(reshape(s.ploss, [12, nCalib, nSite]), 2);

% Alpha
alpha = et ./ PEmean;
alpha(PEmean==0 & et==0) = 1;
alpha(PEmean==0 & et~=0) = 0;

% Beta
beta = r ./ pr;
beta(pr==0 & r==0) = 1;
beta(pr==0 & r~=0) = 0;

% Gamma
gamma = ro ./ pro;
gamma(pro==0 & ro==0) = 1;
gamma(pro==0 & ro~=0) = 0;

% Delta
delta = loss ./ ploss;
delta(ploss==0) = 0;

% Run the soil moisture model over the entire time period
P = reshape(P, [nTime, nSite]);
PE = reshape(PE, [nTime, nSite]);
s = soilMoisture(P, PE, awcs, awcu, ssi, sui, 3, showprogress);

pr = reshape(s.pr, [12, nTime/12, nSite]);
pro = reshape(s.pro, [12, nTime/12, nSite]);
ploss = reshape(s.ploss, [12, nTime/12, nSite]);

% Reshape to monthly
P = reshape(P, [12, nYears, nSite]);
PE = reshape(PE, [12, nYears, nSite]);

% Calculate CAFEC precipitation, monthly departures, and monthly K' factors
Pcafec = (alpha.*PE) + (beta.*pr) + (gamma.*pro) - (delta.*ploss);
D = P - Pcafec;
Dmean = mean( abs(D(:,calib,:)), 2);
w = NaN(1,1,nSite);
for k = 1:nSite
    num = PEmean(:,:,k) + r(:,:,k) + ro(:,:,k);
    denom = Pmean(:,:,k) + loss(:,:,k);
    w(:,:,k) = (num' / denom') + 2.8;
end
Kprime = 1.5 * log10(w ./ Dmean) + 0.5;

% Adjust the K' values and use to calculate Z indices
denom = sum(Dmean .* Kprime, 1);
K = (17.67/denom) .* Kprime;
Z = K .* D;

% Compute PDSI from Z indices
Z = reshape(Z, [nTime, nSite]);
[X, Xm] = zPDSI(Z, showprogress);

% Reshape to match size of input grids
X = reshape(X, fullSize);
X = permute(X, order);
if nargout > 1
    Xm = reshape(Xm, fullSize);
    Xm = permute(Xm, order);
end
if nargout > 2
    Z = reshape(Z, fullSize);
    Z = permute(Z, order);
end
if nargout > 3
    PE = reshape(PE, fullSize);
    PE = permute(PE, order);
end

end