function permAniso = getAnisotropyRatio(vcl, isClayVcl, zf, clayMine, ...
                                        gamma, silt, porozf, idHW)
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
%
% EXAMPLE:
%   Compute permeability anisotropy for material derived from each layer, 
%   with silt particles likely within the clay fraction:
%       silt = 1;
%       permAniso = getAnisotropyRatio(myFaultedSection, fault, silt)
%
%--------------------------------------------------------------
permAniso = ones(1, numel(vcl));
id = find(vcl >= isClayVcl);
if nargin > 7
    assert(numel(zf) == 2)
    zf = [repelem(zf(1), sum(id < idHW(1))), ...
          repelem(zf(2), sum(id >= idHW(1)))];
elseif numel(zf) == 1
    zf = repelem(zf, numel(id));
end

if ~isempty(id)
    % (1) Grain aspect ratio
    if nargin > 7
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
        if strcmp(clayMine, 'kao') || strcmp(clayMine, 'Ill')
            m = repelem(10, 2, 1);
        elseif strcmp(clayMine, 'sme')
            m = repelem(300, 2, 1);
        elseif strcmp(clayMine, 'mic')
            m = repelem(100, 2, 1);
        end
    end
    
    if ~isempty(silt) && silt == 1
        f = [0.9 0.1];
        m = [sqrt(1/(f(1)/m(1)^2 + f(2)/1)); ...
             sqrt(1/(f(1)/m(2)^2 + f(2)/1))];  % m equivalent. Eq. (13) in
                                               % Daigle & Dugan (2011).
    end
    
    if nargin > 7
       m = [repelem(m(1), sum(id < idHW(1))), ...
             repelem(m(2), sum(id >= idHW(1)))]; 
    else
       m = repelem(m(1), numel(id));
    end
    
    
    % (2) Porosity at faulting depth
    if size(porozf, 2) == 2
        poro = porozf(id, 1) + rand(numel(id), 1).*(porozf(id, 2) - porozf(id, 1));
    elseif size(porozf, 2) == 1
        poro = porozf(id);
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
    epsv   = 0.3149*exp(0.0002304*zf) -0.3117*exp(-0.001417*zf);
    theta0 = 45;                            % [deg]
    theta1 = atand((1-epsv)*tand(theta0));  % after compaction
    theta  = acotd(gamma + cotd(theta1));   % after shearing
    
    
    % Permeability anisotropy ratio (k_along / k_across)
    num = 1 + ((8*m./9).*cosd(theta) + ...
          (2/pi)*sind(theta))./(3*pi./(8*(1-poro')) - 0.5);
    den = 1 + ((8*m./9).*sind(theta) + ...
          (2/pi)*cosd(theta))./(3*pi./(8*(1-poro')) - 0.5);
    permAniso(id)  = (num./den).^2;
    
end

end

        