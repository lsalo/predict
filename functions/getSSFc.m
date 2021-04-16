function [SSFc, SSFcBounds] = getSSFc(vcl, isClayVcl, zf, thick, ...
                                      faultDisp, idHW)
% Get critical shale smear factor (SSFc), i.e. the value at which a given 
% smear becomes discontinuous. Note that an object fault with field disp
% must be passed as well.
%
% Key references:
%   Childs et al., GSLSP (2007)
%   Grant, Petroleum Geoscience (2017)
%
% INPUT:
%   fault: Usually, an instance of Fault with property Throw in m.
%   FS: An instance of FaultedSection.
%
% MODELS:
%   Vcl < smear threshold: NaN (it does not apply).
%
%   Vcl >= smear threshold: Function of clay content (Vcl), faulting depth 
%                           (zf) and layer thickness (T). The model builds 
%                           on the algorithm presented by
%                           Grant, Petroleum Geoscience, (2017), and 
%                           extends it by incorporating zf and T:
%          (1) Get base SSFc endpoints (min and max) for each clay content 
%              range according to Fig. 4a in Grant (2017). These curves are 
%              most representative to faulting at 500m depth (see Fig. 3b 
%              and 4b in paper.
%          (2) Include the effect of zf by modifying the endpoints obtained 
%              for the corresponding Vcl range to which each material 
%              belongs. This is done to agree with the experimental results 
%              of Giger et al., AAPG bull. (2013) and 
%              Ciftci et al., AAPG bull. (2013).
%          (3) Compute the mode for the triangular distribution that will 
%              be used to sample a SSFc value from the range determined by 
%              the SSFc endpoints. The mode moves closer to the SSFc_max 
%              endpoint as layer thickness increases. This is done to agree
%              with observations that thicker clay layers generate more 
%              continuous smears. This could be related to increased water 
%              content in the clay, which has been documented in the 
%              laboratory to generate longer (more continuous) smears (e.g.
%              Sperrevik et al. Petroleum Geoscience, 2000).
%          (4) Sample a value of SSFc consistent with Vcl, zf and thickness 
%              of each layer in the faulted stratigraphy.
%
% OUTPUT:
%  SSFc: Critical shale smear factor for each fault material derived from 
%        the footwall and hangingwall stratigraphy (field vals) and
%        min and max bounds for each material (field bounds).
%
% EXAMPLE:
%       SSFc = getSSFc(myFaultedSection, fault)
%
%--------------------------------------------------------------

% Initialize
SSFc = nan(1, numel(vcl));
id = find(vcl >= isClayVcl);

% If smearing clays are present:
if sum(id) > 0
    vclc   = vcl(id);
    thickc = thick(id);
    if nargin > 5
        zf = [repelem(zf(1), sum(id < idHW(1))), ...
              repelem(zf(2), sum(id >= idHW(1)))];
    else
       zf = repelem(zf, numel(id)); 
    end
    
    % 1. SSFc_min and SSFc_max (endpoints) for each layer
    %    Values from Grant (2017), representative for shallow
    %    faulting at around 500m (see Fig. 3b and Giger et al.,
    %    2013).
    endpoints = zeros(2, numel(vclc));
    endpoints(:, vclc<=0.5) = repmat([2; 5], 1, sum(vclc<=0.5));
    endpoints(:, all([vclc>0.5; vclc<=0.6])) = repmat([3; 7], 1, ...
                                            sum(all([vclc>0.5; vclc<=0.6])));
    endpoints(:, all([vclc>0.6; vclc<=0.7])) = repmat([5; 10], 1, ...
                                            sum(all([vclc>0.6; vclc<=0.7])));
    endpoints(:, vclc>0.7) = repmat([7; 12], 1, sum(vclc>0.7));
    
    % 2. Modify endpoints to account for zf. Strong changes
    % for sediments faulted at very shallow depths (<500m) vs
    % mid depth (1-1.5km). Deeper, the changes become less
    % pronounced with depth.
    endpoints(:, zf<=500) = endpoints(:, zf<=500) - (500 - zf(zf<=500))/250;
    endpoints(:, all([zf>500; zf<=1500])) = ...
                            endpoints(:, all([zf>500; zf<=1500])) + ...
                            (zf(all([zf>500; zf<=1500])) - 500)/250;
    endpoints(:, zf>1500) = endpoints(:, zf>1500) + 4 + ...
                            (zf(zf>1500) - 1500)/1000;
    
    % Assign to output
    SSFcBounds = endpoints;
    
    % 3. Compute mode of triangular distribution to sample from.
    thickMax = faultDisp;
    thickMin = faultDisp/50;            % algorithm limit accounting for
                                        % resoultion when meshing fault.
    peak = (1 - min([repelem(1, numel(thickc)); ...
                     ((thickMax - thickc)/(thickMax - thickMin))])) ...
           .* (endpoints(2, :) - endpoints(1, :)) + endpoints(1, :);
    
    
    % 4. Generate samples
    for k = 1:numel(vclc)
        triDist = makedist('Triangular', 'a', endpoints(1, k), ...
                           'b', peak(1, k), 'c', endpoints(2, k));
        SSFc(id(k)) = random(triDist, 1, 1);
    end   
end

end

        