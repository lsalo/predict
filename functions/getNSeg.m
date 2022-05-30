function nSeg = getNSeg(vcl, isClayVcl, zf)
% Get number of along-strike segments.
%
% Key references:
%   Noorsalehi-Garakani et al., JSG (2013)
%   Kettermann et al., Solid Earth (2016)
%   Kettermann et al., JGR: Solid Earth (2017)
%
% INPUT:
%   vcl: tbd (1xN array)
%
% MODELS:
%       Function of average section clay content (Vcl) and faulting depth 
%       (zf). The model is based on number of along-strike smear segments
%       obtained from excavated faults in the lab and field. The number of
%       along-strike segments is obtained via counting NH (the number of
%       holes) along strike-parallel horizontal lines with L = D. Then, 
%       nSeg equals 2NH, since each hole is defined by continuous smear at
%       each side.
%          (1) Include the effect of zf by modifying the nSeg endpoints 
%              obtained from papers above. 
%          (2) Compute the mode for the triangular distribution that will
%              be used to sample a nSeg value from the range determined by
%              the nSeg endpoints. The mode moves closer to the nSeg_min
%              (less segmented) as clay content increases. This is done to 
%              agree with observations that thicker clay layers generate 
%              more  continuous smears.
%          (3) Sample a value of nSeg consistent with Vcl and zf
%
% OUTPUT:
% nSeg structure with the fields below.
%       type:  distribution type (triangular)
%       param: mode of the triangular distribution 
%       range: nSeg endpoints 
%       fcn:   function to compute n values of nSeg consistent with inputs
%              for each layer.
%
% EXAMPLE:
%   vcl = [0.1, 0.3, 0.5, 0.7, 0.9];
%   isClayVcl = 0.4;
%   zf = 500;
%   nSeg = getNSeg(vcl, isClayVcl, zf)
%   vals = nSeg.fcn(n samples)
%
%--------------------------------------------------------------

% Initialize
nSeg.type  = [];
nSeg.param = [];
nSeg.range = [];
nSeg.dist  = [];
nSeg.fcn   = [];

% Mean vals
zfm = mean(zf);
vclm = mean(vcl(vcl > isClayVcl));
assert(zfm <= 3000)
assert(all([vclm >= isClayVcl, vclm <= 1]))

% nSeg min and max based on data compiled from papers above (zf < 1m)
endpoints = [1 16];

% Modify endpoints to account for zf. Strong changes for sediments
% faulted at very shallow depths ( < 500 m). Deeper, the changes become
% less pronounced.
if zfm > 0
    if  zfm <= 500
        endpoints(2) = endpoints(2) - zfm/250;
    elseif zfm <= 3000
        v = zfm - 500;
        endpoints(2) = endpoints(2) - 500/250 - v/500;
    end
end
endpoints = round(endpoints);
nSeg.range = endpoints;
    
% Compute mode based on average vcl
vcl_lim = [isClayVcl 1];
peak = endpoints(1) + (vcl_lim(2) - vclm)/(diff(vcl_lim)) * ...
                      (endpoints(2) - endpoints(1));
nSeg.param = peak;

% Generate fcn
nSeg.type = 'tri';
nSeg.dist = makedist('Triangular', 'a', endpoints(1), 'b', peak, ...
                     'c', endpoints(2));
nSeg.fcn = @(x) random(nSeg.dist, x, 1);
end

        