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
vclm = mean(vcl(isClayVcl));

% nSeg min and max based on data compiled from papers above

% Modify endpoints to account for zf. Strong changes for sediments
% faulted at very shallow depths ( < 500 m). Deeper, the changes become
% less pronounced.

% Compute mode based on average vcl

% Generate fcn.

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

        