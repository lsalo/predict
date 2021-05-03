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
%   Although other controls exist, such as clay content (i.e.
%   rheology) and juxtapositions in the deformed section, it is
%   difficult to include them rigorously due to lack of data.
%   Hence, we consider only the dependency of fault thickness
%   (T) on fault displacement (D).
%   We consider that D/T is between 10 and 1000 and most of the
%   data is around D/T = 100, which agrees with the published
%   literature. Hence, we sample from a beta distribution with
%   parameters a, b = 5. This gives a bounded distribution,
%   symmetric, and with very few data in the upper and lower
%   20%. We scale each beta random sample to the interval
%   [log10(10), log10(1000)].
%
% OUTPUT:
%  fcn:  Fault core thickness anonymous fcn
%  type: distribution type (beta)
%
% EXAMPLE:
%       myFault = myFault.getFaultThickness
%
%--------------------------------------------------------------

thick.rangeDT = [10, 1000];
limLog = log10([10, 1000]);
[a, b] = deal(5);
thick.fcn = @(D) D ./ ( 10.^(limLog(1) + ...
                        betarnd(a, b, numel(D), 1).*diff(limLog)) );
thick.type = 'beta';
thick.param = [a, b];
end