function SSFc = getSSFc(vcl, isClayVcl, zf, thick, faultDisp, idHW)
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
N = numel(vcl);
SSFc.type  = cell(1, N);
SSFc.param = cell(1, N);
SSFc.range = cell(1, N);
SSFc.dist = cell(1, N);
SSFc.fcn   = cell(1, N);

% Expand faulting depths
if nargin > 5
    zf = [repelem(zf(1), idHW(1)-1), repelem(zf(2), numel(idHW))];
else
    zf = repelem(zf, numel(vcl));
end
assert(numel(zf) == numel(vcl));

% Maximum and minimum available thickness for triangular mode placement
thickMax = faultDisp;
thickMin = faultDisp/50;            % algorithm limit accounting for
                                    % resoultion when meshing fault. 

for n = 1:N
    if vcl(n) >= isClayVcl    
        % 1. SSFc_min and SSFc_max (endpoints) for each layer
        %    Values from Grant (2017), representative for shallow
        %    faulting at around 500m (see Fig. 3b and Giger et al.,
        %    2013).
        if vcl(n) <= 0.5
            endpoints = [2, 5];
        elseif all([vcl(n)>0.5; vcl(n)<=0.6])
            endpoints = [3, 7];
        elseif all([vcl(n)>0.6; vcl(n)<=0.7])
            endpoints = [5, 10];
        else
            endpoints = [7, 12];
        end

        % 2. Modify endpoints to account for zf. Strong changes
        % for sediments faulted at very shallow depths (<500m) vs
        % mid depth (1-1.5km). Deeper, the changes become less
        % pronounced with depth.
        if zf(n) <= 500
            endpoints = endpoints - (500 - zf(n))/250;
        elseif all([zf(n) > 500; zf(n) <= 1500])
            endpoints = endpoints + (zf(n) - 500)/250;
        else 
            endpoints = endpoints + ((zf(n) - 1500)/1000 + 4);
        end

        % Assign to output
        SSFc.range{n} = endpoints;

        % 3. Compute mode of triangular distribution to sample from.
        assert(thick(n) <= thickMax)
        peak = (1 - min([1; ((thickMax - thick(n))/(thickMax - thickMin))])) ...
               .* (endpoints(2) - endpoints(1)) + endpoints(1);
        SSFc.param{n} = peak;

        % 4. Generate fcn
        SSFc.type{n} = 'tri';
        SSFc.dist{n} = makedist('Triangular', 'a', endpoints(1), 'b', peak, ...
                                'c', endpoints(2));
        SSFc.fcn{n} = @(x) random(SSFc.dist{n}, x, 1);
    end
end

end

        