function phi = getResidualFrictionAngle(vcl, FS)
% Get residual friction angle (degrees) of each fault material.
%
% INPUT:
%  FS: An instance of FaultedSection.
%
% MODELS:
%   vcl < 0.2: We simply use random integers from a uniform discrete
%         distribution between 25 and 35. This is because the variation is
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
% NOTE:
%   We allow for user-defined residual friction angles. Hence, if the 
%   footwall or hangingwall objects have ResFric values, we will use those 
%   instead. See the documentation of the Stratigraphy class for more 
%   information.
%
% OUTPUT:
%  Residual friction angle of each layer in obj (footwall layers first from
%  bottom to top, and hangingwall layers follow, also from bottom to top).
%
% EXAMPLE:
%   phi = getResidualFrictionAngle(myFaultedSection)
%--------------------------------------------------------------

% Initialize
phi = zeros(1,numel(vcl));

% (1) Random with uniform P for Vcl < 0.2
ids = vcl < 0.2;
N = sum(ids);
phi(ids) = randi([28, 38], 1, N);

% (2) Vcl > 0.2
a = 3;
b = 5;
N = sum(~ids);
vclc = vcl(vcl >= 0.2);

% Add fit bounds (utils/fitClayResFricSke64Mes86.m)
ci = [27.6374  -12.5232   25.7510   -1.9403; ...
      49.6431   -6.5168   33.4577   -1.5849];
lower = ci(1,1).*exp(ci(1,2).*vclc) + ci(1,3).*exp(ci(1,4).*vclc);
upper = ci(2,1)*exp(ci(2,2)*vclc) + ci(2,3)*exp(ci(2,4)*vclc);

% Random data and variable phi
rVals = betarnd(a, b, 1, N);
phi(~ids) = lower + rVals.*(upper - lower);

% Overwrite (if values were provided in Stratigraphy objects)
if nargin > 1 && ( ~isempty(FS.FW.ResFric) || ~isempty(FS.HW.ResFric) )
    assert(numel(FS.FW.ResFric) == FS.FW.NumLayers && ...
           numel(FS.HW.ResFric) == FS.HW.NumLayers, ...
           strcat('if ResFric is provided, both FW and HW', ...
           ' must be 1xN arrays, where N is the', ...
           ' number of layers in each object. NaN', ...
           ' entries are allowed.'))
    
    isVal = [~isnan(FS.FW.ResFric), ~isnan(FS.HW.ResFric)];
    phi(isVal) = [FS.FW.ResFric(~isnan(FS.FW.ResFric)), ...
    FS.HW.ResFric(~isnan(FS.HW.ResFric))];
end

end