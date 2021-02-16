function[x, xm] = zPDSI(z)
%% Calculates PDSI and modified PDSI from Z indices
%
% [x, xm] = zPDSI(Z)
%
% ----- Inputs -----
%
% Z: A set of monthly Z indices. (nMonths x nSites)
%
% ----- Outputs -----
%
% x: Palmer drought severity index for each month. (nMonths x nSites)
%
% xm: Modified Palmer drought severity index in each month. (nMonths x nSites)

% Preallocate
[nTime, nSite] = size(z);
x1 = zeros(nTime, nSite);
x2 = zeros(nTime, nSite);
x3 = zeros(nTime, nSite);
probEnd = NaN(nTime, nSite);

useX1 = false(nTime, nSite);
useX2 = false(nTime, nSite);
useX12 = false(nTime, nSite);

nump = zeros(1, nSite);
newPeriod = false(1, nSite);
inSubperiod = false(1, nSite);

% Divide by three
z3 = z/3;

% Initialize values for the first month
x1(1,:) = max(0, z3(1,:));
x2(1,:) = min(0, z3(1,:));
x3(1,:) = z3(1,:);

Vlast = zeros(1, nSite);
lastProb = zeros(1, nSite);

iswet = z3(1,:)>=1;
isdry = z3(1,:)<=-1;
isnormal = z3(1,:)>-1 & z3(1,:)<1;

newPeriod(iswet | isdry) = true;
x3(1, isnormal) = 0;
useX12(1, isnormal) = true;

% Integrate the PDSI model over time
for t = 2:nTime
    
    % Get the starting status of each site
    normal = isnormal;
    wet = iswet;
    dry = isdry;
    
    % Reset status for the next month
    isnormal(:) = false;
    iswet(:) = false;
    isdry(:) = false;
    
    % Update the coefficients
    x1(t, :) = max(0, 0.897 * x1(t-1,:) + z3(t,:));
    x2(t, :) = min(0, 0.897 * x2(t-1,:) + z3(t,:));
    x3(t, :) = 0.897 * x3(t-1,:) + z3(t,:);
    
    %% Normal month
    
    % Staying normal
    staynormal = normal & x1(t,:)<1 & x2(t,:)>-1;
    useX12(t, staynormal) = true;
    isnormal(staynormal) = true;
    x3(t, staynormal) = 0;
    
    % Change to wet period
    nextwet = normal & x1(t,:)>=1;
    x3(t, nextwet) = x1(t, nextwet);
    iswet(nextwet) = true;
    newPeriod(nextwet) = true;
    inSubperiod(nextwet) = false;
    nump(nextwet) = 0;
    
    % Change to dry period
    nextdry = normal & x2(t,:)<=-1;
    x3(t, nextdry) = x2(t, nextdry);
    isdry(nextdry) = true;
    newPeriod(nextdry) = true;
    inSubperiod(nextdry) = false;
    nump(nextdry) = 0;
    
    %% New periods
    
    % Reset X1/X2 for new wet/dry periods    
    newwet = wet & newPeriod;
    x1(t, newwet) = 0;
    newdry = dry & newPeriod;
    x2(t, newdry) = 0;
    newPeriod(newwet | newdry) = false;
    
    %% Wet and dry subperiods
    
    % Get sign for dry vs wet calculations
    sign = ones(1, nSite);
    sign(wet) = -1;
    
    % Calculate effective wetness
    effWet = z(t,:) + sign*0.15;
    
    % Get subperiod indices
    nosub = ~inSubperiod & ((wet & effWet>=0) | (dry & effWet<=0));
    newsub = ~inSubperiod & ((wet & effWet<0) | (dry & effWet>0));
    oldsub = inSubperiod;
    insub = newsub | oldsub;
    
    % Get the probability that the period has ended
    probEnd(t, nosub) = 0;
    
    Ze = -2.691 * x3(t-1,:) - sign*1.5;
    V = effWet;
    Q = Ze;
    V(oldsub) = V(oldsub) + Vlast(oldsub);
    Q(oldsub) = Q(oldsub) + Vlast(oldsub);
    probEnd(t, insub) = V(insub) ./ Q(insub);
    
    probEnd(t, probEnd(t,:)<0) = 0;
    probEnd(t, probEnd(t,:)>1) = 1;
    
    % Update the subperiods
    inSubperiod(newsub) = true;
    
    %% Wet and dry periods, next state
    
    % Wet/Dry period has no dry/wet subperiod
    full = (wet|dry) & probEnd(t,:)==0;
    subend = full & lastProb>0;
    fullwet = full & wet;
    fulldry = full & dry;
    
    iswet(fullwet) = true;
    isdry(fulldry) = true;
    x1(t, fullwet) = 0;
    x2(t, fulldry) = 0;
    
    x1(t, fulldry & subend) = 0;
    x2(t, fullwet & subend) = 0;
    inSubperiod(subend) = false;
    V(subend) = 0;
    
    % Change to normal
    change = (wet|dry) & probEnd(t,:)==1;
    nextnormal = change & ((wet & x2(t,:)>-1) | (dry & x1(t,:)<1));
    isnormal(nextnormal) = true;
    x3(t, nextnormal) = 0;
    if any(nextnormal)
        sites = find(nextnormal);
        for k = 1:numel(sites)
            s = sites(k);
            useX12(t-nump(s):t, s) = true;
        end
    end
    
    % Change to opposite state
    flip = change & ((wet & x2(t,:)<=-1) | (dry & x1(t,:)>=1));
    newPeriod(flip) = true;
    inSubperiod(flip) = false;
    
    nextdry = flip & wet;
    isdry(nextdry) = true;
    x3(t, nextdry) = x2(t, nextdry);
    if any(nextdry)
        sites = find(nextdry);
        for k = 1:numel(sites)
            s = sites(k);
            useX2(t-nump(s):t, s) = true;
        end
    end
    
    nextwet = flip & dry;
    iswet(nextwet) = true;
    x3(t, nextwet) = x1(t, nextwet);
    if any(nextwet)
        sites = find(nextwet);
        for k = 1:numel(sites)
            s = sites(k);
            useX1(t-nump(s):t, s) = true;
        end
    end
    
    % Persist in subperiod
    persist = (wet|dry) & ~full & ~change;
    staywet = persist & wet;
    staydry = persist & dry;
    iswet(staywet) = true;
    isdry(staydry) = true;
    
    % Update nump
    nump(full) = 0;
    nump(change) = 0;    
    nump(persist) = nump(persist) + 1;
    
    % Save V for probability calculations in the next month
    Vlast = V;
    lastProb = probEnd(t,:);
end

% Get the standard PDSI
x = x3;
x(useX1) = x1(useX1);
x(useX2) = x2(useX2);

maxwet = abs(x1) >= abs(x2);
maxX1 = useX12 & maxwet;
maxX2 = useX12 & ~maxwet;
x(maxX1) = x1(maxX1);
x(maxX2) = x2(maxX2);

% Get the modified PDSI
modify = probEnd>0 & probEnd<1;
moddry = modify & x3<0;
modwet = modify & x3>0;

xm = x;
xm(moddry) = probEnd(moddry) .* x1(moddry) + (1-probEnd(moddry)) .* x3(moddry);
xm(modwet) = probEnd(modwet) .* x2(modwet) + (1-probEnd(modwet)) .* x3(modwet);

end