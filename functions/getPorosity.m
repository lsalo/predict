function [poro] = getPorosity(vcl, isclayVcl, z, zeval, zf, isUndercompacted, idHW)
%
%
%

% Initialize
N          = numel(vcl);
poro.type  = cell(1, N);
poro.param = cell(1, N);
poro.range = cell(1, N);
poro.fcn   = cell(1, N);

% Expand faulting detph
if nargin < 7
    zf = repelem(zf(1), numel(vcl));
else
    zf  = [repelem(zf(1), idHW(1)-1), ...
           repelem(zf(2), numel(idHW))];
end
if strcmp(zeval, 'zf')
    z = zf;
end

for n=1:N
    % Depth at which porosity should be computed
    if strcmp(zeval, 'zf')
        poro.param{n} = 'zf';
    elseif strcmp(zeval, 'zmax')
        poro.param{n} = 'zmax';
    end
    
    if vcl(n) < isclayVcl   % ideal packing model (Revil et al., JGR, 2002)
        % Sand end member (Eq 9)
        b = 6*10^(-8);                              % [1/Pa] compaction coeff of the sand end-member
        phi_0 = 0.49 + rand(1, 1)*0.1;              % [-] depositional poro of sand end-member
        phi_r = 0.2556*exp(-5.028*10^(-4)*zf(n));   % [-] residual porosity of sand end-member
        rho_g = 2650;                               % [kg/m^3] bulk density of sand grains
        rho_w = 1050;                               % ["]      bulk density of pure water
        g = 9.806;                                  % [m/s^2] gravitational acceleration
        zm = 1 / ( (1-phi_r)*g*b*(rho_g - rho_w) ); % characteristic length
        num = phi_0 - phi_r + (1-phi_0) * phi_r * exp(z(n) / zm);
        den = phi_0 - phi_r + (1-phi_0) * exp(z(n) / zm);
        sPoro = num /den;
       
        % Clay end-member (Eq 10). Here we use a generic kaolinite value,
        % bc this model makes sand porosity more accurate than not using the
        % clay fraction at all. However, most of these factors are
        % uncertain. Hence, we do not follow this process for the clay
        % layers(below), and use the full envelope to determine porosity 
        % and permeability ranges.
        phi0  = 0.55;                                   % depositional porosity
        bc    = 4*10^(-8);                              % 1/Pa "compaction coefficient of the shale end-member"
        zmc   = 1/(g*bc*(rho_g-rho_w));                 % characteristic length [m]
        cPoro = phi0./(phi0 + (1-phi0)*exp(z(n)./zmc)); % clay end-member porosity [-]
        
        % Porosity of an ideal sand-clay mixture (Eq. 11 and 12)
        if vcl(n) < sPoro
            poro.range{n} = repelem(sPoro - vcl(n)*(1-cPoro), 2);
        else
            poro.range{n} = repelem(vcl(n)*cPoro , 2);
        end        
        poro.type{n} = 'det';                           % deterministic
        poro.fcn{n} = @(x) repelem(poro.range{n}(1), x, 1);
        
    else        
        if isUndercompacted == 0
            poro.range{n} = [(1 - nthroot(z(n)./(6.02*1000), 6.35)) - 0.05, 0];
            poro.range{n}(2) = poro.range{n}(1) + 0.1;
            poro.type{n} = 'unif';
            poro.fcn{n} = @(x) poro.range{n}(1) + rand(x, 1) .* ...
                               (poro.range{n}(2) - poro.range{n}(1));
        elseif isUndercompacted == 1
            poro.range{n} = repelem(1 - nthroot(z(n)./(15*1000), 8), 1, 2);
            poro.type{n} = 'det';
            poro.fcn{n} = @(x) repelem(poro.range{n}(1), x, 1);
        else
            error('Inputs not supported. Check documentation.')
        end
    end
end



end