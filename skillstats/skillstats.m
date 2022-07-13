function S = skillstats(r, f)
%SKILLSTATS Calculate skill metric statistics for a dataset vs a reference
%
% S = skillstats(r, f)
%
% This function calculates a variety of statistics useful in assessing the
% skill of a model dataset relative to a reference dataset.  These
% statistics can be used to produce skill summary diagrams like Taylor
% diagrams (Taylor 2001, J Geophys Res Atmos 106:7183?7192) or target
% diagrams (Jolliff et al. 2009, J Mar Syst 76:64?82). 
%
% Input variables:
%
%   r:  n x 1 array, values of reference (observed) variable
%
%   f:  n x m array, values of modeled/predicted variable, for m different
%       models 
%
% Output variables:
%
%   S:  1 x 1 structure with the following fields; all are 1 x m+1 arrays,
%       where the first element corresponds to the reference data and 2:end
%       correspond to the models.
%
%       std:        standard deviation, normalized to n.
%
%       cor:        correlation coefficient
%
%       rmsd:       root mean squared difference
%
%       crmsd:      centered pattern (i.e. unbiased) root mean squared
%                   difference 
%
%       bias:       bias, i.e. average error
%
%       stdnorm:    normalized standard deviation
%
%       rmsdnorm:   root mean squared difference, normalized to standard
%                   deviation of reference data
%
%       crmsdnorm:  centered root mean squared difference, normalized to
%                   standard deviation of reference data 
%
%       aae:        average absolute error
%
%       ri:         reliability index (factor by which model differs from
%                   reference... note that this metric falls apart if any
%                   of the modeled values are 0).
%
%       mef:        modeling efficiency (skill relative to average of
%                   observations, 1 = perfect, 0 = same as averaging obs,
%                   <1 = worse than just averaging observations) 

% Copyright 2013-2017 Kelly Kearney


narginchk(2,2);
validateattributes(r, {'numeric'}, {'column'}, 'skillstats', 'r', 1);
validateattributes(f, {'numeric'}, {'2d'}, 'skillstats', 'f', 2);

[npt, nf] = size(f);
if ~isvector(r) || ~isequal(npt, length(r))
    error('Length of r (input 1) must match number of rows in f (input 2)');
end

isn = isnan(r) | any(isnan(f),2);
if any(isn)
    warning('NaNs found; corresponding samples removed from both datasets');
    r = r(~isn);
    f = f(~isn,:);
end

% Standard deviation

S.std = [std(r,1) std(f,1)];

% Correlation coef

fbar = mean(f);
rbar = mean(r);

S.cor(1) = 1;
for ii = 1:nf
    S.cor(ii+1) = (sum((f(:,ii) - fbar(ii)) .* (r - rbar))./npt)./(S.std(1).*S.std(ii+1));
end

% RMSD

S.rmsd = sqrt(sum((bsxfun(@minus, [r f], r)).^2)/npt);

% Centered pattern RMSD

S.crmsd = zeros(1,nf+1);
S.crmsd(1) = 0;
for ii = 1:nf
    S.crmsd(ii+1) = sqrt(sum(((f(:,ii) - fbar(ii)) - (r - rbar)).^2)./npt);
end

% Bias

S.bias = [0 bsxfun(@minus, mean(f,1), mean(r))];

% Normalized stats

S.stdnorm   = S.std./S.std(1);
S.rmsdnorm  = S.rmsd./S.std(1);
S.crmsdnorm = S.crmsd./S.std(1);
S.biasnorm  = S.bias./S.std(1);

% Extras from Stow et al 2009

S.aae = sum(abs(bsxfun(@minus, [r f], r)))./npt;

S.ri = exp(sqrt(sum((log(bsxfun(@rdivide, r, [r f]))).^2))./npt);

S.mef = (sum((r - rbar).^2) - sum((bsxfun(@minus, [r f], r)).^2))./(sum((r - rbar).^2));

