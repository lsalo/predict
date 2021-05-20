function thick = getFaultThickness()
% Get fault core thickness.
%
% Key references:
%   Childs et al., GSLSP (2007)
%   Childs et al., JSG (2009)
%
% INPUT:
%   fault: An instance of Fault
%   Nsim:  Number of realizations
%
% MODEL:
%   Fault rock thickness is known to primarily depend on
%   displacement (see papers above and references therein).
%   We consider that D/T is between 10 and 1000 and most of the
%   data is around D/T = 100, which agrees with the published
%   literature. Hence, we sample from a beta distribution with
%   parameters a, b = 5. This gives a bounded distribution,
%   symmetric, and with very few data in the upper and lower
%   20%. We scale each beta random sample to the interval
%   [log10(10), log10(1000)].
%
% OUTPUT:
%  Structure thick with fields:
%       type: distribution type (beta)
%       param: shape parameters of the beta distribution
%       range: thickness range (function of fault displacement)
%       fcn:  function to compute thickness, by passing a scalar fault
%             displacement value or a Nx1 array of displacement values.
%
% EXAMPLE:
%       D     = repelem(100, 1000, 1);
%       thick = getFaultThickness();
%       vals  = thick.fcn(D);
%
%--------------------------------------------------------------

thick.rangeDT = [10, 1000];
limLog = log10(thick.rangeDT);
[a, b] = deal(5);
thick.type = 'beta';
thick.param = [a, b];
thick.range = @(D) [D./thick.rangeDT(2) D./thick.rangeDT(1)];
thick.fcn = @(D) D ./ ( 10.^(limLog(1) + ...
                        betarnd(a, b, numel(D), 1).*diff(limLog)) );
end