function[PE] = pethorn(T, years, lats, calib)
%% Calculates potential evapotranspiration using the Thornthwaite equation
%
% PE = pethorn(T, years, lats, calib)
%
% ----- Inputs -----
%
% T: Monthly temperatures in Celsius. (12 x years x sites)
%
% years: The years of the temperature data. A vector the length of the
%    second dimension of T.
%
% lats: Latitudes of the sites. A vector the length of the third dimension
%    of T.
%
% calib: Indicates which years to use for the calibration period. A logical
%    vector the length of years.
%
% ----- Outputs -----
%
% PE: Potential evapotranspiration in each month. (12 x years x sites)

% Get the monthly means for the calibration period
Tmean = mean(T(:,calib,:), 2);

% Apply Sellers equation to get heat index and exponent a
Tmean(Tmean<0) = 0;
Tmean(Tmean>26.5) = 26.5;
I = sum((Tmean / 5) .^ 1.514, 1);
a = 1E-6 * (0.675*I.^3 - 77.1*I.^2 + 17920*I + 492390);

% Compute PE. Adjust for extrememly hot and below freezing temperatures
T(T>38) = 38;
PE = 16 * (T .* (10./I)) .^ a;
PE(T<=0) = 0;
PE(:,:,I==0) = 0;

% Use a table lookup to get unadjusted PE for temperatures above 26.5
[PEhot, Thot] = unadjustedPE;
hot = T > 26.5;
S = interp1(Thot, PEhot, T(hot));

% Convert unadjusted PE from mm/day to mm/month
[month, y, ~] = ind2sub(size(T), find(hot));
daysPerMonth =[31 28 31 30 31 30 31 31 30 31 30 31]';
S = S .* daysPerMonth(month);

leap = month==2 & mod(years(y), 4)==0;
S(leap) = S(leap) * (29/28);
PE(hot) = S;

% Load table of mean monthly possible duration of sunlight
dayz = daylength;

% Adjust PE for sunlight duration. Poleward of 50 degrees, use the values
% at 50 degrees.
lats(lats>50) = 50;
lats(lats<-50) = -50;
sunFactor = interp1(dayz.lats, dayz.dayz, lats) / 30;
sunFactor = permute(sunFactor, [2 3 1]);
PE = PE .* sunFactor;

leap = mod(years,4)==0;
PE(2,leap,:) = PE(2,leap,:) * (29/28);

% Convert to inches and return as monthly time array
PE = PE / 25.4;

end