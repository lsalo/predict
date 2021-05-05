function permAniso = getAnisotropyRatio(vcl, zf, clayMine, idHW)
% Get permeability anisotropy [m^2] for smear fault materials. Note that an 
% object fault with fields disp (displacement in meters) and thick 
% (thickness in meters) must be passed as well.
%
% INPUT:
%   fault: usually, an instance of Fault with fields disp and thick.
%   FS:   an instance of FaultedSection.
%   silt (optional): If 0 (default), it is considered that no silt-sized 
%                    particles are present within the clay fraction of any 
%                    of the faulted materials. If 1, we then account for 
%                    the potential presence of silt particles within the 
%                    clay fraction. This leads to much smaller anisotropy 
%                    ratios (i.e. 1-3 instead of the typical 5-15). See 
%                    Daigle & Dugan, WRR (2011) for further details.
%   isUndercompacted (optional): pass 1 for undercompacted shales.
%
% MODELS:
%   Vcl < smear threshold: 1 (Perm_across = Perm_along).
%
%   Vcl >= smear threshold: Function of predominant clay mineral and 
%                           porosity.
%          We use the model of Daigle & Dugan, WRR (2011).
%          (1) Get grain aspect ratio (m) based on dominant clay mineral 
%              and proportion of minerals in Vcl i.e. only clay, or silt 
%              likely).
%          (2) Compute porosity at faulting depth. We use a compaction 
%              curve (Depth to porosity range). By default, we use the 
%              Baldwin-Butler, AAPG Bull. (1985) +-5% limits, which cover 
%              several other curves for "normally compacted argillaceous 
%              sediments". We sample randomly according to a uniform 
%              probability distribution within these limits. In case of 
%              thick shales and/or likely undercompaction due to e.g.
%              overpressure, the Dickinson, AAPG Bull. (1953) curve is 
%              used. In that case, the isUndercompacted argument must be 
%              passed with value 1.
%          (3) Compute theta as a result of compaction and faulting (i.e. 
%              shear strain).
%          (4) Compute permeability anisotropy ratio 
%              (Perm_along / Perm_across) following Eq. (19) in the paper.
%
% OUTPUT:
%  permAniso: Permeability anisotropy of each fault material derived from 
%             each layer in the footwall and hangingwall.
%  for fcn, gamma must be a scalar, poro can be a 1xN array or a scalar. 
%
% EXAMPLE:
%   Compute permeability anisotropy for material derived from each layer, 
%   with silt particles likely within the clay fraction:
%       silt = 1;
%       permAniso = getAnisotropyRatio(myFaultedSection, fault, silt)
%
%--------------------------------------------------------------

% Initialize
N = numel(vcl);
permAniso.type  = cell(1, N);
permAniso.param = cell(1, N);
permAniso.range = cell(1, N);
permAniso.fcn   = cell(1, N);

if nargin > 3
    % Expand faulting depths
    zf = [repelem(zf(1), idHW(1)-1), repelem(zf(2), numel(idHW))];
    
    % (1) Grain aspect ratio
    assert(numel(clayMine) == 2);           % predominant in FW and HW
    m = zeros(2, 1);
    if strcmp(clayMine{1}, 'kao') || ...
            strcmp(clayMine{1}, 'Ill')           % kaolinite, illite
        m(1) = 10;                          % Daigle & Dugan
    elseif strcmp(clayMine{1}, 'sme')       % smectite group
        m(1) = 300;                         % can get to 100-500
    elseif strcmp(clayMine{1}, 'mic')       % micas, chlorite
        m(1) = 100;                         % can get to 50-200
    end
    
    if strcmp(clayMine{2}, 'kao') || strcmp(clayMine{2}, 'Ill')
        m(2) = 10;
    elseif strcmp(clayMine{2}, 'sme')
        m(2) = 300;
    elseif strcmp(clayMine{2}, 'mic')
        m(2) = 100;
    end
    
else
    zf = repelem(zf, numel(vcl));
    
    if strcmp(clayMine, 'kao') || strcmp(clayMine, 'Ill')
        m = repelem(10, 2, 1);
    elseif strcmp(clayMine, 'sme')
        m = repelem(300, 2, 1);
    elseif strcmp(clayMine, 'mic')
        m = repelem(100, 2, 1);
    end
    
end
assert(numel(zf) == numel(vcl));
 
for n=1:N     
    if vcl(n) > 0.95
        warning(['Vcl of layer(s) with vcl > 0.98 set to 0.95 for ' ...
                 'perm anisotropy calculation.']);
       vcl(n) = 0.95;
    end
    
    % Grain aspect ratio
    f = [vcl(n), 1-vcl(n)];
    if nargin > 3
        if n < idHW(1)
            meq = sqrt(1/(f(1)/m(1)^2 + f(2)/1));   % m equivalent. Eq. (13) in
                                                    % Daigle & Dugan (2011).
        else
            meq = sqrt(1/(f(1)/m(2)^2 + f(2)/1));
        end
    else
        meq = sqrt(1/(f(1)/m(1)^2 + f(2)/1));
    end

    % Clay grain orientation after compaction and shearing (i.e. 45 degrees 
    % is totally random, 0 degrees is completely flat leading to max 
    % anisotropy).
    % First, we should rigorously estimate the uniaxial strain. This 
    % requires some data, and, also, it turns out that the grain 
    % orientation before shearing is not that important once shear strain 
    % is above 5 (always true here, since shear strain is fault 
    % displacement over fault thickness, and that ratio is considered to be 
    % at least 10). Even if that were not the case, if the clay contains
    % some silt (very likely) it  also does not matter much in terms of the 
    % final  anisotropy. So, we provide a very rough approximation only.
    epsv   = 0.3149*exp(0.0002304*zf(n)) -0.3117*exp(-0.001417*zf(n));
    theta0 = 45;                             % [deg]
    theta1 = atand((1-epsv)*tand(theta0));   % after compaction
    theta  = @(g) acotd(g + cotd(theta1));   % after shearing, g = gamma
    
    
    % Permeability anisotropy ratio (k_along / k_across)
    permAniso.type = 'det';         % deterministic (assuming known poro)
    
    num = @(g, porov) 1 + ((8*meq/9)*cosd(theta(g)) + ...
                              (2/pi)*sind(theta(g))) ./ ...
                             (3*pi./(8*(1-porov)) - 0.5);
    den = @(g, porov) 1 + ((8*meq/9)*sind(theta(g)) + ...
                              (2/pi)*cosd(theta(g))) ./ ...
                             (3*pi./(8*(1-porov)) - 0.5);
    permAniso.fcn{n}  = @(gamma, porovals) (num(gamma, porovals) ./ ...
                                            den(gamma, porovals)).^2;    
end

end

        