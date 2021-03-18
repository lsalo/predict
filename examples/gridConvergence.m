%% Grid convergence analysis

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
thickness = {[20 10 20 10 30 10; 20 10 20 10 30 10], ...
             [25 25 25 25; 25 25 25 25]};
vcl       = {[0.05 0.4 0.1 0.5 0.15 0.6; 0.2, 0.7, 0.25, 0.8, 0.3, 0.9], ...
             [0.8 0.3 0.5 0; 0.3, 0.7, 0.15, 0.6]};
dip       = {[0, 0], [0, 0]};
faultDip  = [45, 60];

% Optional Input parameters
nl   = [size(vcl{1}, 2), size(vcl{2}, 2)];   % just for convenience here
zf   = {[100, 100], [500, 500]};    % m
zmax = {[repelem(2000, 1, nl(1)); repelem(2000, 1, nl(1))], ...
        [repelem(1000, 1, nl(2)); repelem(1000, 1, nl(2))]};
cm = {'kao', 'sme'};            % predominant clay mineral
maxPerm = [];                   % cap max perm? [mD]
siltInClay = [true, false];     % consider 25% of silt fraction in clay?

% Flow upscaling options
U.useAcceleration = 1;          % 1 requires MEX and AMGCL setup
U.method          = 'mpfa';     % 'tpfa' recommended if useAcc. = 0
U.outflux         = 0;          % compare outflux of fine and upsc. models
U.ARcheck         = 0;          % check if Perm obtained with grid with 
                                % Aspect Ratio of only 3 gives same Perm
                                
% Initialize target grid cell sizes and result containers
gridRes = repmat([2, 20; 1, 10; 0.5, 5; 0.2, 2; 0.1, 1; 0.05, 0.5], 1, 1, 2);         
Ngrid = size(gridRes, 1);
Nstrat = numel(vcl);
poro  = zeros(Ngrid, Nstrat);
perm  = zeros(Ngrid, 3, Nstrat);
cellDim = zeros(Ngrid, 2, Nstrat);

tic
for k=1:Nstrat
% FW and HW
footwall = Stratigraphy(thickness{k}(1,:), vcl{k}(1,:), dip{k}(1), ...
                        'DepthFaulting', zf{k}(1), ...
                        'DepthBurial', zmax{k}(1,:), 'ClayMine', cm{k});
hangingwall = Stratigraphy(thickness{k}(2,:), vcl{k}(2,:), dip{k}(2), ...
                           'IsHW', 1, 'NumLayersFW', footwall.NumLayers, ...
                           'DepthFaulting', zf{k}(2), ...
                           'DepthBurial', zmax{k}(2,:), 'ClayMine', cm{k});

% Strati in Faulted Section
mySect = FaultedSection(footwall, hangingwall);

    % Run grid loop for each strati case
    for n=1:Ngrid
    % Reset random seed (same values in MatProps and smear placement)
    rng('default');

    % Generate fault object with properties for each realization
    myFault = Fault(mySect, faultDip(k), gridRes(n, :, k));

    % Get dependent variables
    myFault = myFault.getMaterialProperties(mySect, 'maxPerm', maxPerm, ...
        'siltInClay', siltInClay(k));

    % Generate smear object with T, Tap, L, Lmax
    Tap = getApparentThick(mySect, myFault.Dip);
    smear = Smear(mySect.Vcl, mySect.IsClayVcl, mySect.Thick, Tap, ...
                  mySect.DepthFaulting, myFault, 1, mySect);

    % Compute upscaled permeability distribution
    myFault = myFault.upscaleSmearPerm(mySect, smear, U);

    % Save result
    poro(n, k) = myFault.Poro;
    perm(n, :, k) = myFault.Perm;
    cellDim(n, :, k) = myFault.Grid.CellDim;
    end
    
disp(['Case ' num2str(k) 'out of ' num2str(Nstrat) ' completed.'])
end
toc

%% Output analysis

% Plot permeabilities vs 1/h
hL = cellDim(:,2);
A  = cellDim(:,1).*cellDim(:,2);
latx = {'Interpreter','latex'};
k_md = perm./(milli*darcy);
limy = [floor(log10(fix(min(min(k_md))))), ...
        ceil(log10(fix(max(max(k_md)))))];

f1 = figure(1);
tiledlayout(Nstrat, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile
plot(1./hL, k_md(:,1), '-ok', 'MarkerSize', 4)
hold on
plot(1/hL(5), k_md(5,1), '-ok', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6)
hold off
grid on
xlabel('$1/h_\mathrm{L}$ [m$^{-1}$]', latx{:}, 'fontSize', 12)
ylabel('$\hat{k}_{xx}$ [mD]', latx{:}, 'fontSize', 12)
xlim([0.01 12])
xticks([0.01 0.1 1 10])
xticklabels({'0.01', '0.1' '1' '10'})
set(gca,'XScale','log', 'YScale', 'log')
ylim(10.^limy)
yticks(10.^(limy(1):limy(2)))


nexttile
plot(1./hL, k_md(:,2), '-or', 'MarkerSize', 4)
hold on
plot(1/hL(5), k_md(5,2), '-or', 'MarkerFaceColor', [255, 125, 125]/255, 'MarkerSize', 6)
hold off
grid on
xlabel('$1/h_\mathrm{L}$ [m$^{-1}$]', latx{:}, 'fontSize', 12)
ylabel('$\hat{k}_{yy}$ [mD]', latx{:}, 'fontSize', 12)
xlim([0.01 12])
xticks([0.01 0.1 1 10])
xticklabels({'0.01', '0.1' '1' '10'})
set(gca,'XScale','log', 'YScale', 'log')
ylim(10.^limy)
yticks(10.^(limy(1):limy(2)))

nexttile
plot(1./hL, k_md(:,3), '-ob')
hold on
plot(1/hL(5), k_md(5,3), '-ob', 'MarkerFaceColor', [125, 125, 255]/255, 'MarkerSize', 6)
hold off
grid on
xlabel('$1/h_\mathrm{L}$ [m$^{-1}$]', latx{:}, 'fontSize', 12)
ylabel('$\hat{k}_{zz}$ [mD]', latx{:}, 'fontSize', 12)
xlim([0.01 12])
xticks([0.01 0.1 1 10])
xticklabels({'0.01', '0.1' '1' '10'})
set(gca,'XScale','log', 'YScale', 'log')
ylim(10.^limy)
yticks(10.^(limy(1):limy(2)))
set(f1, 'position', [500, 500, 700, 200]);