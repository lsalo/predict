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
Nsim = [10 20 50 100 200];% 500 1000 2000 5000];
k_md = cell(numel(Nsim), 1);
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
k_md{k} = perm ./ (milli*darcy);

disp(['Simulation ' num2str(k) ' (' num2str(Nsim(k)) ' realizations) done.']) 
disp([num2str(numel(Nsim) - k) ' simulations remaining.'])
end
toc


%% Output analysis
% Compare perm distros

% Plot utils
sz = [14, 12];
latx = {'Interpreter','latex'};

% Hist params
nbins = 25;
lims = [floor(log10(min(min(k_md{end})))), ceil(log10(max(max(k_md{end}))))];
edges = logspace(lims(1), lims(2), nbins);
colrs = repmat(hsv(numel(Nsim)), 1, 1);

f1 = figure(1);
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
for k=1:numel(Nsim)
   nexttile(1)
   hold on
   histogram(k_md{k}(:, 1), edges, 'Normalization', 'probability', ...
             'DisplayStyle', 'stairs', 'EdgeColor', colrs(k, :), ...
             'DisplayName', ['$N_\mathrm{sim}$ = ' num2str(Nsim(k))])
   xlabel('$k_{xx}$ [mD]', latx{:}, 'fontSize', sz(2))
   ylabel('P [-]', latx{:}, 'fontSize', sz(2))
   xlim([10^lims(1) 10^lims(2)])
   ylim([0 1])
   grid on
   set(gca,'XScale','log')
   leg = legend(latx{:}, 'fontSize', sz(2), 'location', 'northwest');
   set(leg.BoxFace, 'ColorType','truecoloralpha', ...
       'ColorData', uint8(255*[1;1;1;.6])); 
   
   %nexttile(2)
   
   
   %nexttile(3)
end
set(f1, 'position', [200, 200, 700, 350]);
