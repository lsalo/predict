function perm = getPermeability(vcl, isclayVcl, zf, zmax, cap, idHW)
% Get permeability [m^2] across each faulted material.
%
% INPUT:
%   FS: an instance of FaultedSection.
%   cap (optional): max perm val (mD) for sand-based material.
%   isUndercompacted (optional): pass 1 for undercomp. shales.
%
% MODELS:
%   Vcl < smear threshold: Function of Vcl, zmax and zf.
%         (1) We use Sperrevik et al., NPSSP (2002)'s model, based on Vcl, 
%             max burial depth and faulting depth. This model provides a 
%             value for k_across.
%         (2) We pick a random term from a uniform probability distribution 
%             bounded by k_across -+ 1 order of magnitude (OM).
%
%   Vcl >= smear threshold: Function of Vcl and zmax. Same approach as in 
%                           Grant, GSLSP (2019), but using different data.
%          (1) Compaction curve (Depth to porosity range). By default, we 
%              use the Baldwin-Butler, AAPG Bull. (1985) +-5% limits, which 
%              cover several other curves for "normally compacted
%              argillaceous sediments". In case of thick shales and/or 
%              likely undercompaction due to e.g. overpressure, the 
%              Dickinson, AAPG Bull. (1953) curve is used. In that case, 
%              the isUndercompacted argument must be passed with value 1.
%          (2) Poro vs perm correlation (get perm range). We use the model 
%              by Yang and Aplin., MPG (2010) to obtain the 
%              bedding-perpendicular permeability based on clay content and 
%              porosity. We obtain two bounds based on (1).
%          (3) Sample from the perm range using a uniform probability 
%              distribution.
%
% NOTE:
%   We allow for user-defined permeability values. Hence, if the footwall 
%   or hangingwall objects have Perm values, we will use those instead. See 
%   the documentation of the Stratigraphy class for more information.
%
% OUTPUT:
%  perm: Permeability of the fault material derived from each layer in the 
%        stratigraphy object. We output a function of N for each layer. N 
%        permeability samples from the corresponding perm distribution can
%        then be obtained by passing N (see example below).
%
% EXAMPLE:
%   Compute permeability distribution functions for material derived from 
%   each layer:
%       perm = getPermeability(myFaultedSection)
%   Get N samples for material derived from k layer, where indices k run 
%   from 1 (bottom of footwall), to obj.FW.NumLayers + obj.HW.NumLayers 
%   (top of HW):
%       N = 100;
%       k = 3;
%       samples = perm{k}(N)  (units = m^2)
%
%--------------------------------------------------------------

% Utilities, initialize and check if values are passed
md_to_m2 = 9.869232667160129e-16;
N = numel(vcl);

perm.type  = cell(1, N);
perm.param = cell(1, N);
perm.range = cell(1, N);
perm.fcn   = cell(N, 1);

% (I) Vcl < threshold for smearing
if nargin < 6
    zf = repelem(zf(1), N);
else
    zf = [repelem(zf(1), idHW(1)-1), repelem(zf(2), numel(idHW))];
end

for n=1:N
    if vcl(n) < isclayVcl    
        % Sand Perm (Sperrevik et al., 2002)
        a    = [8*10^4, 19.4, 0.00403, 0.0055, 12.5];
        permSand  = (a(1)*exp(-(a(2)*vcl(n) + a(3)*zmax(n) + ...
                    (a(4)*zf(n) - a(5)).*(1-vcl(n)).^7)));    % [mD]
        permSand = [0.5.*permSand, 5.*permSand];              % min & max
        perm.type{n} = 'unif';

        % Cap values
        if ~isempty(cap) && cap ~= 0
            permSand(permSand > cap) = cap;
        end
        permSand = log10(permSand*md_to_m2);
        perm.range{n} = 10.^permSand;
        
        % Fcn
        perm.fcn{n} = @(x) 10.^(permSand(1) + rand(x, 1) * ...
                                (permSand(2) - permSand(1)));   % [m^2]
    
    else                % (II) Vcl >= threshold for smearing
        % Clay volume fraction to mass fraction
        % Unless the mineralogy of the sand and clay is very diff,
        % the density of sand and clay is dependent on their poro
        % difference only. This porosity difference will rarely be
        % larger than ~0.2, hence the clay volume fraction and mass
        % volume fraction are similar.
        sandPoro = 0.49 / (exp(zmax(n)/(3.7*1000)));
        rho_s = 2650;                           % "avg" density of solid part.
        rho_w = 1040;                           % "avg" formation water density
        rho_clay = @(p) p*rho_w + (1-p)*rho_s;  % p = poro
        rho_mat  = @(p) (1-vcl(n))*(sandPoro*rho_w + (1-sandPoro)*rho_s) + ...
                        vcl(n)*(rho_clay(p));
        mcl = @(p) vcl(n) * rho_clay(p) ./ rho_mat(p);
        
        % Get permeability fcns
        e = @(plim) plim ./ (1 - plim);         % plim = upper or lower bound              
        f1 = @(p) -69.59 -26.79*mcl(p) + 44.07*mcl(p).^0.5;
        f2 = @(p, plim) (-53.61 -80.03*mcl(p) + 132.78*mcl(p).^0.5)*e(plim).^1;
        f3 = @(p, plim) (86.61 + 81.91*mcl(p) -163.61*mcl(p).^0.5)*e(plim).^0.5;
        permc = @(p, plim) f1(p) + f2(p, plim) + f3(p, plim); % ln(perm [m^2])
        
        perm.type{n} = 'unif';
        perm.fcn{n} = @(p, pmin, pmax) exp(permc(p, pmin) + ...
                                           rand(numel(p), 1) .* ...
                                           (permc(p, pmax) - permc(p, pmin)));
        perm.range{n} = @(pmin, pmax) exp([permc(pmin, pmin) ...
                                           permc(pmax, pmax)]);
                                       
    end

end
% Overwrite perm for layers with passed perm
% if nargin == 7
%     idPassed = find([~isnan(FS.FW.Perm) ~isnan(FS.HW.Perm)]);
%     permPassed = [FS.FW.Perm(~isnan(FS.FW.Perm)), ...
%                   FS.HW.Perm(~isnan(FS.HW.Perm))];
%     for k = 1:numel(idPassed)
%         perm{idPassed(k)} = @(N) (permPassed(k) * md_to_m2) * ones(N, 1);
%     end
% end

end
