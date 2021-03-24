%% Required number of simulations 
%  Goal: to obtain perm distros that are representative of the full 
%        parameter uncertainty. This code incrementally adds the number of
%        realizations (N) to compare the output permeability distributions 
%        as N increases.
%% Example 0: Single stratigraphic case + analysis
%
% Use parfor instead of for when running several simulations.
% 

clear
close all

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
thickness = {[20 10 20 10 30 10], [20 10 20 10 30 10]};
vcl       = {repmat([0.05 0.4 0.1 0.5 0.15 0.6], 1, 1), ...
             repmat([0.2, 0.7, 0.25, 0.8, 0.3, 0.9], 1, 1)};
dip       = [0, 0];
faultDip  = 60;

% Optional Input parameters
zf   = [500, 500];              % m
maxPerm = [];                   % mD
siltInClay = true;              % consider 25% of silt fraction in clay     

% Flow upscaling options
U.useAcceleration = 1;          % requires MEX and AMGCL setup
U.method          = 'mpfa';     % 'tpfa' recommended if useAcc. = 0
U.outflux         = 0;          % compare outflux of fine and upsc. models
U.ARcheck         = 0;          % check if Perm obtained with grid with 
                                % Aspect Ratio of only 3 gives same Perm

% FW and HW
footwall = Stratigraphy(thickness{1}, vcl{1}, dip(1), ...
                        'DepthFaulting', zf(1));
hangingwall = Stratigraphy(thickness{2}, vcl{2}, dip(2), 'IsHW', 1, ...
                           'NumLayersFW', footwall.NumLayers, ...
                           'DepthFaulting', zf(2));

% Strati in Faulted Section
mySect = FaultedSection(footwall, hangingwall);

% Prepare loop
Nsim = [10 20 50 100 200 500 1000];% 2000 5000];
perms = cell(numel(Nsim), 1);
for k=1:numel(Nsim) 
perm = nan(Nsim(k), 3);

parfor n=1:Nsim(k)    % parallel computing toolbox required for parfor
    myFault = Fault(mySect, faultDip);
    
    % Get dependent variables
    myFault = myFault.getMaterialProperties(mySect, 'maxPerm', maxPerm, ...
                                            'siltInClay', siltInClay);
    
    % Generate smear object with T, Tap, L, Lmax
    Tap = getApparentThick(mySect, myFault.Dip);
    smear = Smear(mySect.Vcl, mySect.IsClayVcl, mySect.Thick, Tap, ...
                  mySect.DepthFaulting, myFault, 1, mySect);
    
    % Compute upscaled permeability distribution
    myFault = myFault.upscaleSmearPerm(mySect, smear, U);
    
    % Save result
    perm(n, :) = myFault.Perm;
end
perms{k} = perm;

disp(['Simulation ' num2str(k) ' (' num2str(Nsim(k)) ' realizations) done.']) 
disp([num2str(numel(Nsim) - k) ' simulations remaining.'])
end
toc


%% Output analysis
% Compare perm distros
