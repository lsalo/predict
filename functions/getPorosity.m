function [poro] = getPorosity(vcl, isclayVcl, zf, zmax, zclay, isUndercompacted, idHW)
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

for n=1:N
    if vcl(n) < isclayVcl   % ideal packing model (Revil et al., JGR, 2002)
        b = 6*10^(-8);                               % [1/Pa] compaction coeff of the sand end-member
        phi_0 = 0.49 + rand(1, 1).*0.1;              % [-] depositional poro of sand end-member
        phi_r = 0.2556*exp(-5.028*10^(-4).*zf(n)');  % [-] residual porosity of sand end-member
        rho_g = 2650;                                % [kg/m^3] bulk density of sand grains
        rho_w = 1050;                                % ["]      buld density of pure water
        g = 9.806;                                   % [m/s^2] gravitational acceleration
        z = zmax(n);                                 % maximum depth
        zm = 1 ./ ( (1-phi_r).*g*b*(rho_g - rho_w) );   % characteristic length
        num = phi_0 - phi_r + (1-phi_0) .* phi_r .* exp(z ./ zm);
        den = phi_0 - phi_r + (1-phi_0) .* exp(z ./ zm);
        poro.range{n} = repmat(num ./ den, 1, 2);
        % TO ADD: vcl should be included (not endmember poro always). Just use
        % a generic (e.g. kaolinite) for poro clay and compute using ideal
        % packing.
        
        poro.type{n} = 'det';       % deterministic
        poro.fcn{n} = @(x) repelem(poro.range{n}(1), x, 1);
        
    else
        if strcmp(zclay, 'zf')
            z = zf(n);
            poro.param{n} = 'zf';
        elseif strcmp(zclay, 'zmax')
            z = zmax(n);
            poro.param{n} = 'zmax';
        end
        
        if isUndercompacted == 0
            poro.range{n} = [(1 - nthroot(z./(6.02*1000), 6.35)) - 0.05, 0];
            poro.range{n}(2) = poro.range{n}(1) + 0.1;
            poro.type{n} = 'unif';
            poro.fcn{n} = @(x) poro.range{n}(1) + rand(x, 1) .* ...
                               (poro.range{n}(2) - poro.range{n}(1));
        elseif isUndercompacted == 1
            poro.range{n} = repelem(1 - nthroot(z./(15*1000), 8), 1, 2);
            poro.type{n} = 'det';
            poro.fcn{n} = @(x) repelem(poro.range{n}(1), x, 1);
        else
            error('Inputs not supported. Check documentation.')
        end
    end
end



end