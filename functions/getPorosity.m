function [poro] = getPorosity(vcl, isclayVcl, zf, zmax, zclay, isUndercompacted, idHW)
%
%
%

ids = find(vcl < isclayVcl);
idc = find(vcl >= isclayVcl);
poro = zeros(numel(vcl), 2);
if nargin < 7
    zf = repelem(zf, numel(vcl));
else
    zf  = [repelem(zf(1), idHW(1)-1), ...
           repelem(zf(2), numel(idHW))];
end

if numel(ids) > 0           % ideal packing model (Revil et al., JGR, 2002)
    b = 6*10^(-8);                              % [1/Pa] compaction coeff of the sand end-member
    phi_0 = 0.49 + rand(numel(ids), 1).*0.1;    % [-] depositional poro of sand end-member
    phi_r = 0.2556*exp(-5.028*10^(-4).*zf(ids)');% [-] residual porosity of sand end-member
    rho_g = 2650;                               % [kg/m^3] bulk density of sand grains
    rho_w = 1050;                               % ["]      buld density of pure water
    g = 9.806;                                  % [m/s^2] gravitational acceleration
    z = zmax(ids)';                             % maximum depth
    zm = 1 ./ ( (1-phi_r).*g*b*(rho_g - rho_w) );   % characteristic length
    num = phi_0 - phi_r + (1-phi_0) .* phi_r .* exp(z ./ zm);
    den = phi_0 - phi_r + (1-phi_0) .* exp(z ./ zm);
    poro(ids, :) = repmat(num ./ den, 1, 2);
    % TO ADD: vcl should be included (not endmember poro always). Just use
    % a generic (e.g. kaolinite) for poro clay and compute using ideal
    % packing.
    
end

if numel(idc) > 0                % compaction curves
    if strcmp(zclay, 'zf')
        z = zf(idc);
    elseif strcmp(zclay, 'zmax')
        z = zmax(idc);
    end
    
    if isUndercompacted == 0
        poro(idc, 1) = (1 - nthroot(z./(6.02*1000), 6.35)) - 0.05;
        poro(idc, 2) = poro(idc, 1) + 0.1;
    elseif isUndercompacted == 1
        poro(idc, :) = 1 - nthroot(z./(15*1000), 8);
    else
        error('Inputs not supported. Check documentation.')
    end
end



end