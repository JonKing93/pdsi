function[s] = soilMoisture(P, PE, awcs, awcu, ssi, sui, type, showprogress)
%% Calculates soil moisture
%
% s = soilMoisture(P, PE, awcs, awcu, ssi, sui, 1)
% Calculates soil moisture in the surface and upper layers.
%
% s = soilMoisture(P, PE, awcs, awcu, ssi, sui, 2)
% Calculates et, r, pr, ro, pro, loss, and ploss
%
% s = soilMoisture(P, PE, awvs, awvu, ssi, sui, 3)
% Calculates pr, pr, and ploss
%
% ----- Inputs -----
%
% P: A matrix of monthly precipitation data (in inches). Each column is a
%    site and each row is one monthly time step. (nMonths x nSite)
%
% PE: A matrix of monthly potential evapotranspiration data (in inches). Each
%    column is a site and each row is one monthly time step. (nMonths x nSite)
%
% awcs: Available water capacity in the surface layer (inches). If a
%    scalar, uses the same value for each site. Otherwise, a vector with
%    one element per site.
%
% awcu: Available water capacity in the underlying layer (inches). If a
%    scalar, uses the same value for each site. Otherwise, a vector with
%    one element per site.
%
% ssi: Starting soil moisture in surface layer (inches).  If a scalar, uses
%    uses the same value for each site. Otherwise, a vector with one
%    element per site.
%
% sui: Starting soil moisture in underlying layer (inches). If a scalar,
%    uses the same value for each site. Otherwise a vector with one element
%    per site.
%
% type: Switch for the outputs required
%    1 - Only need soil moisture in the final year
%    2 - Need terms for CAFEC coefficients
%    3 - Need pro, pr, and ploss for Z indices
%
% showprogress: A scalar logical indicating whether to show a progress bar
%
% ----- Outputs -----
%
% s: An output structure. All fields are (nMonths x nSite)
%    r: Recharge, combined layers
%    pr: Potential recharge
%    ro: Runoff
%    pro: Potential runoff
%    loss: Loss
%    ploss: Potential loss
%    et: Estimated actual evapotranspiration
%    ssi: Initial surface soil moisture in each time step
%    sui: Initial underlying soil moisture in each time step

% ----- Written By -----
% Original script by Dave Meko.
% Updated and vectorized by Jonathan King, 2020.

% Default progressbar
if ~exist('showprogress','var') || isempty(showprogress)
    showprogress = false;
end

% Error check data types of inputs
args = {P, PE, awcs, awcu, ssi, sui};
names = {'P','PE','awcs','awcu','ssi','sui'};
for k = 1:numel(args)
    assert(isnumeric(args{k}), sprintf('%s must be numeric', names{k}));
    assert(~any(isnan(args{k}) | isinf(args{k}) | ~isreal(args{k}), 'all'), ...
        sprintf('%s cannot contain NaN, Inf, or complex values', names{k}));
    assert(all(args{k}>=0, 'all'), sprintf('%s cannot have negative values', names{k}));
end

% Error check P and PE. Get sizes
assert(ismatrix(P), 'P must be a matrix');
assert(ismatrix(PE), 'PE must be a matrix');
assert(isequal(size(P), size(PE)), 'P and PE must have the same size');
[nMonths, nSite] = size(P);

% Error check water capacity and initial soil moisture sizes. Ensure water
% capacity vectors are row vectors
for k = 3:numel(args)
    assert(isscalar(args{k}) | (isvector(args{k}) & numel(args{k})==nSite), ...
        sprintf('%s must either be a scalar or a vector with one element per site (%.f)', names{k}, nSite));
end
awcs = awcs(:)';
awcu = awcu(:)';

% Preallocate outputs
ss = NaN(nMonths+1, nSite);
su = NaN(nMonths+1, nSite);
if type == 2
    r = NaN(nMonths, nSite);
    losses = NaN(nMonths, nSite);
end
if type>1
    ro = NaN(nMonths, nSite);
end

% Preallocate slices
dels = NaN(1, nSite);
delu = NaN(1, nSite);
es = NaN(1, nSite);
if type == 2
    rs = NaN(1, nSite);
    ru = NaN(1, nSite);
    eu = NaN(1, nSite);
end
if type > 1
    ro_m = NaN(1, nSite);
end

% Compute water balance deficit and water capacity
d = PE - P;
awc = awcu + awcs;
if isscalar(awc)
    awc = repmat(awc, [1, nSite]);
end

% Initialize soil moistures
ss(1,:) = ssi;
su(1,:) = sui;

% Optionally display progress bar
if showprogress
    percent = 0;
    str = 'Running soil moisture model: ';
    message = sprintf('%s%.f%%', str, percent);
    step = ceil(nMonths/100);
    f = waitbar(0, message);
end

% Run the moisture model over each month. Start by getting the empty
% capacity in each time step
for m = 1:nMonths
    sempty = awcs - ss(m,:);
    uempty = awcu - su(m,:);
    
    % Slice the rows (improves scaling for large data sets)
    d_m = d(m,:);
    ss_m = ss(m,:);
    su_m = su(m,:);
    
    % Surface calculations
    loss = d_m >= 0;
    dels(loss) = -min(d_m(loss), ss_m(loss));
    dels(~loss) = min(sempty(~loss), -d_m(~loss));
    es(loss) = -dels(loss);
    es(~loss) = 0;
    if type==2
        rs(loss) = 0;
        rs(~loss) = dels(~loss);
    end
    if type > 1
        ro_m(loss) = 0;
    end
    
    % Underlayer calculations
    excess = -d_m - dels;
    gain = ~loss & excess>0;
    delu(gain) = min(uempty(gain), excess(gain));
    eu(~gain) = (d_m(~gain) - es(~gain)) .* su_m(~gain) ./ awc(~gain);
    eu(~gain) = min(eu(~gain), su_m(~gain));
    delu(~gain) = -max(eu(~gain), 0);
    if type == 2
        eu(gain) = 0;
        ru(gain) = delu(gain);
        ru(~gain) = 0;
    end
    if type > 1
        ro_m(gain) = max( excess(gain)-uempty(gain), 0 );
        ro_m(~gain) = 0;
    end
    
    % Update the soil moisture in each layer, runoff, and recharge
    ss(m+1,:) = ss(m,:) + dels;
    su(m+1,:) = su(m,:) + delu;
    if type == 2
        r(m,:) = rs + ru;
        losses(m,:) = es + eu;
    end
    if type > 1
        ro(m,:) = ro_m;
    end
    
    % Update the waitbar
    if showprogress && (mod(m,step)==0 || m==nMonths)
        percent = 100 * m/nMonths;
        message = sprintf('%s%.f%%', str, percent);
        waitbar(m/nMonths, f, message);
    end
end
if showprogress
    close(f);
end

% Evapotranspiration
if type == 2
    et = P - ro - diff(ss,1) - diff(su,1);
end

% Remove the last moisture entry
ss = ss(1:end-1,:);
su = su(1:end-1,:);

% Potential recharge, loss and runoff
if type > 1
    pro = ss + su;
    pr = awc + pro;
    plosss = min(PE, ss);
    plossu = (PE - plosss) .* (su ./ awc);
    ploss = plosss + plossu;
end

% Build the output structure
if type == 1
    s = struct('ssi',ss,'sui',su);
elseif type == 2
    s = struct('et',et,'r',r,'pr',pr,'ro',ro,'pro',pro,'loss',losses,'ploss',ploss);
elseif type == 3
    s = struct('pro',pro,'pr',pr,'ploss',ploss);
end

end
