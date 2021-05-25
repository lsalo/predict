function phi = getResidualFrictionAngle(vcl)
% Get residual friction angle (degrees) of each fault material.
%
% INPUT:
%  FS: An instance of FaultedSection.
%
% MODELS:
%   vcl < 0.2: We simply use random integers from a uniform discrete
%         distribution between 28 and 38. This is because the variation is
%         small enough (we could use 30 to 35) and also there is not a 
%         reliable model to predict it based on the Stratigraphy properties.
%
%   vcl >= 0.2: Function of Vcl.
%         (1) We use data in:
%               - Skempton, Geotechnique (1964)
%               - Mesri & Cepeda-Diaz, Geotechnique (1986)
%             to fit an exponential function.
%         (2) We incorporate a stochastic term following a beta 
%             distribution centered around the value obtained from the
%             regression model in (1). We use parameters a = 3 and b = 5 so 
%             that it's relatively narrow (i.e. not many values at 
%             extremes) and somewhat skewed to the smaller values (more 
%             datapoints).
%         *See details in utils/fitClayResFricSke64Mes86.m
%
% OUTPUT:
%  Structure phi with fields:
%       type: distribution type (unif or beta)
%       param: shape parameters of the beta distribution
%       range: thickness range 
%       fcn:  function to compute n values of phi consistent with each
%             value provided in vcl.
%
% EXAMPLE:
%   n = 1000;
%   vcl = [0.1, 0.3, 0.5, 0.7, 0.9];
%   phi = getResidualFrictionAngle(vcl)
%   vals = cell2mat(cellfun(@(x) x(n), phi.fcn, 'uniformOutput', false));
%--------------------------------------------------------------

% Initialize
N = numel(vcl);
phi.type  = cell(1, N);
phi.param = cell(1, N);
phi.range = cell(1, N);
phi.fcn   = cell(1, N);

for n=1:N
    if vcl(n) < 0.2         % Random with uniform P for Vcl < 0.2
        phi.type{n}  = 'unif';
        phi.range{n} = [28, 38];
        phi.fcn{n}   = @(x) phi.range{n}(1) + rand(x, 1).*diff(phi.range{n});
    
    else                    % Vcl > 0.2
        phi.type{n}  = 'beta';
        a = 3; b = 5;
        phi.param{n} = [a, b];

        % Add fit bounds (utils/fitClayResFricSke64Mes86.m)
        ci = [27.6374  -12.5232   25.7510   -1.9403; ...
              49.6431   -6.5168   33.4577   -1.5849];
        lower = ci(1,1).*exp(ci(1,2).*vcl(n)) + ci(1,3).*exp(ci(1,4).*vcl(n));
        upper = ci(2,1)*exp(ci(2,2)*vcl(n)) + ci(2,3)*exp(ci(2,4)*vcl(n));
        phi.range{n} = [lower, upper];
        
        % Random data and variable phi
        phi.fcn{n} = @(x) lower + betarnd(a, b, x, 1).*(upper - lower);
    end
end

% Overwrite (if values were provided in Stratigraphy objects)
% if nargin > 1 && ( ~isempty(FS.FW.ResFric) || ~isempty(FS.HW.ResFric) )
%     assert(numel(FS.FW.ResFric) == FS.FW.NumLayers && ...
%            numel(FS.HW.ResFric) == FS.HW.NumLayers, ...
%            strcat('if ResFric is provided, both FW and HW', ...
%            ' must be 1xN arrays, where N is the', ...
%            ' number of layers in each object. NaN', ...
%            ' entries are allowed.'))
%     
%     isVal = [~isnan(FS.FW.ResFric), ~isnan(FS.HW.ResFric)];
%     phi(isVal) = [FS.FW.ResFric(~isnan(FS.FW.ResFric)), ...
%     FS.HW.ResFric(~isnan(FS.HW.ResFric))];
% end

end