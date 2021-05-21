%% Example 0: Single stratigraphic case + analysis
%
% Use parfor instead of for when running several simulations.
% 

clear
close all force

%% Mrst modules. 
% Run startup (enough for grid) + add required for flow upscaling. With 
% either development or release version of MRST.
%
% We use the following MRST utilities:
%   * merge_options.m -->
%   *
%   * 
mrstModule add mrst-gui coarsegrid upscaling incomp mpfa


%% Define model and upscale permeability
% Mandatory Input parameters
%           {[FW], [HW]}
thickness = {repelem(25, 1, 4), [5 10 15 10 20 10 10 5 15]};
vcl       = {[0.1 0.4 0.2 0.5], ...
             [0.3 0.6 0.1 0.7 0.2 0.8 0.3 0.9 0.1]};
dip       = [0, 0];
faultDip  = 70;
Nsim      = 1000;               % Number of simulations/realizations

% Optional Input parameters
zf      = [1000, 1000];         % [m]
maxPerm = 1000;                 % [mD]
rho     = 0.7;                  % Corr. coeff. for multivariate distributions

% Flow upscaling options
U.useAcceleration = 1;          % requires MEX and AMGCL setup
U.method          = 'mpfa';     % 'tpfa' recommended if useAcc. = 0
U.outflux         = 0;          % compare outflux of fine and upsc. models
U.ARcheck         = 0;          % check if Perm obtained with grid with 
                                % Aspect Ratio of only 3 gives same Perm

% FW and HW
footwall = Stratigraphy(thickness{1}, vcl{1}, 'Dip', dip(1), ...
                        'DepthFaulting', zf(1), 'DepthBurial', repelem(2500, 1, 4));
hangingwall = Stratigraphy(thickness{2}, vcl{2}, 'Dip', dip(2), 'IsHW', 1, ...
                           'NumLayersFW', footwall.NumLayers, ...
                           'DepthFaulting', zf(2), 'DepthBurial', repelem(2500, 9));

% Instantiate FaultedSection object (Strati in Faulted Section)
mySect = FaultedSection(footwall, hangingwall, faultDip, 'maxPerm', maxPerm);

% Get material property distributions
mySect = mySect.getMatPropDistr();

% Generate fault object with properties for each realization
faults = cell(Nsim, 1);
smears = cell(Nsim, 1);
tic
parfor n=1:Nsim    % parfor allowed if you have the parallel computing toolbox
    myFault = Fault(mySect, faultDip);
    
    % Get material property samples
    myFault = myFault.getMaterialProperties(mySect, 'corrCoef', rho);
    
    % Generate smear object with T, Tap, L, Lmax
    smear = Smear(mySect.Vcl, mySect.IsClayVcl, mySect.Thick, ...
                  mySect.Tap, mySect.DepthFaulting, myFault, 1, mySect);
    
    % Compute upscaled permeability distribution
    myFault = myFault.upscaleSmearPerm(mySect, smear, U);
    
    % Save result
    faults{n} = myFault;
    smears{n} = smear;
    if mod(n, 100) == 0
        disp(['Simulation ' num2str(n) ' / ' num2str(Nsim) ' completed.'])
    end
end
toc


%% Output analysis
% Visualize Strati, fault thickness of 1st realization
mySect.plotStrati(faults{1}.MatProps.thick, faultDip);  

% Histograms for each MatProp (all sims, we select one stratigraphic layer)
% This should plot for all realizations that contain the given id.
%layerId = 1;                                            
%plotMatPropsHist(faults, smears, mySect, layerId) 

% MatProps correlations
[R, pval] = plotMatPropsCorr(faults, mySect, 3);

% General fault materials and perm view
%plotId = selectSimId('minX', faults, Nsim);                % simulation index
%faults{plotId}.plotMaterials(mySect) 

% Plot upscaled Poro and Perm (all sims, 3 directions)
plotUpscaledPerm(faults)
