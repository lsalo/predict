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
name      = {'A', 'B', 'C', 'D', 'E'};
thickness = {[repelem(10, 1, 10); repelem(10, 1, 10)], ...
             [25 25 25 25; 25 25 25 25], [50 50; 50 50], ...
             [20 10 20 10 30 10; 20 10 20 10 30 10], ...
             [20 40 20 20; 40 20 20 20]};
vcl       = {[repmat([0.2 0.6], 1, 5); repmat([0.5 0.3], 1, 5)], ...
             [0.8 0.3 0.5 0; 0.3, 0.7, 0.15, 0.6], [0.5 0.1; 0.5 0.1], ...
             [0.05 0.4 0.1 0.5 0.15 0.6; 0.2, 0.7, 0.25, 0.8, 0.3, 0.9], ...
             [0 0.4 0.3 0.6; 0.1 0.4 0.2 0.6]};
dip       = {[0, 0], [0, 0], [10, 20], [10, 5], [0, 0]};
faultDip  = [45, 60, 65, 75, 55];

% Optional Input parameters
nl   = [size(vcl{1}, 2), size(vcl{2}, 2), size(vcl{3}, 2), ...
        size(vcl{4}, 2), size(vcl{5}, 2)];   % just for convenience here
zf   = {[100, 100], [500, 500], [1000, 1000], [100, 100], [2000, 2000]};    % m
zmax = {[repelem(800, 1, nl(1)); repelem(800, 1, nl(1))], ...
        [repelem(1000, 1, nl(2)); repelem(1000, 1, nl(2))], ...
        [repelem(1000, 1, nl(3)); repelem(1000, 1, nl(3))], ...
        [repelem(2000, 1, nl(4)); repelem(2000, 1, nl(4))], ...
        [repelem(3000, 1, nl(5)); repelem(3000, 1, nl(5))]};
cm = {'kao', 'sme', 'ill', 'mic', 'kao'};    % predominant clay mineral
maxPerm = [];                   % cap max perm? [mD]
siltInClay = [true, false, false, true, true];     % silt fraction in clay?

% Flow upscaling options
U.useAcceleration = 1;          % 1 requires MEX and AMGCL setup
U.method          = 'mpfa';     % 'tpfa' recommended if useAcc. = 0
U.outflux         = 0;          % compare outflux of fine and upsc. models
U.ARcheck         = 0;          % check if Perm obtained with grid with 
                                % Aspect Ratio of only 3 gives same Perm
                                
% Initialize target grid cell sizes and result containers
Nstrat = numel(vcl);
gridRes = repmat([1, 10; 0.5, 5; 0.2, 2; 0.1, 1; 0.05, 0.5; ...
                  0.02, 0.2; 0.01, 0.1], 1, 1, Nstrat);         
Ngrid = size(gridRes, 1);
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
    rng(5)

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
    %disp(myFault.MatProps)
    disp(['Grid ' num2str(n) ' / ' num2str(Ngrid) ' completed.'])
    end
 
disp('-------------------------------------------------------------------')
disp(['Case ' num2str(k) ' out of ' num2str(Nstrat) ' finished.'])
disp('-------------------------------------------------------------------')
end
toc

%% Output analysis

% Plot permeabilities vs 1/h
hL = cellDim(:,2);
A  = cellDim(:,1).*cellDim(:,2);
k_md = perm./(milli*darcy);

% Plotting utilities
sz = [14, 12];
latx = {'Interpreter','latex'};
%colrs = [128 0 0; 0 130 200; 255 225 128; 0 0 128; 0 0 0] ./ 255;

f1 = figure(1);
tiledlayout(3, Nstrat, 'Padding', 'compact', 'TileSpacing', 'compact');
for k=1:Nstrat
limy = [floor(log10(min(min(k_md(:,1,k))))), ...
        ceil(log10(max(max(k_md(:,1,k))))); ...
        floor(log10(min(min(k_md(:,2,k))))), ...
        ceil(log10(max(max(k_md(:,2,k))))); ...
        floor(log10(min(min(k_md(:,3,k))))), ...
        ceil(log10(max(max(k_md(:,3,k)))))];
    
nexttile(k)
hold on
plot(1./hL, k_md(:,1,k), '-o', 'color', 'k', 'MarkerSize', 4, ...
     'DisplayName', name{k})
plot(1/hL(4), k_md(4,1,k), '-o', 'color', 'k', 'MarkerFaceColor', ...
     [0.5 0.5 0.5], 'MarkerSize', 6, 'HandleVisibility','off')
hold off
grid on
if k == 1
    xlabel('$1/h_\mathrm{L}$ [m$^{-1}$]', latx{:}, 'fontSize', 12)
    ylabel('$k_{xx}$ [mD]', latx{:}, 'fontSize', 12)
end
xlim([0.05 12])
xticks([0.1 1 10])
xticklabels({'0.1' '1' '10'})
set(gca,'XScale','log', 'YScale', 'log')
ylim(10.^limy(1, :))
yticks(10.^(limy(1, 1):limy(1, 2)))
%leg = legend(latx{:}, 'fontSize', sz(2), 'location', 'northwest');
%set(leg.BoxFace, 'ColorType','truecoloralpha', ...
%    'ColorData', uint8(255*[1;1;1;.6]));
title(name{k}, latx{:}, 'fontSize', sz(2));

nexttile(k + Nstrat)
hold on
plot(1./hL, k_md(:,2,k), '-o', 'color','r', 'MarkerSize', 4)
plot(1/hL(4), k_md(4,2,k), '-o', 'color', 'r', ...
     'MarkerFaceColor', [255, 125, 125]/255, 'MarkerSize', 6)
hold off
grid on
if k==1
    %xlabel('$1/h_\mathrm{L}$ [m$^{-1}$]', latx{:}, 'fontSize', 12)
    ylabel('$k_{yy}$ [mD]', latx{:}, 'fontSize', 12)
end
xlim([0.05 12])
xticks([0.1 1 10])
xticklabels({'0.1' '1' '10'})
set(gca,'XScale','log', 'YScale', 'log')
ylim(10.^limy(2, :))
yticks(10.^(limy(2, 1):limy(2, 2)))


nexttile(k + 2*Nstrat)
hold on
plot(1./hL, k_md(:,3,k), '-o', 'color', 'b', 'MarkerSize', 4)
plot(1/hL(4), k_md(4,3,k), '-o', 'color', 'b', ...
     'MarkerFaceColor', [125, 125, 255]/255, 'MarkerSize', 6)
hold off
grid on
if k == 1
    %xlabel('$1/h_\mathrm{L}$ [m$^{-1}$]', latx{:}, 'fontSize', 12)
    ylabel('$k_{zz}$ [mD]', latx{:}, 'fontSize', 12)
end
xlim([0.05 12])
xticks([0.1 1 10])
xticklabels({'0.1' '1' '10'})
set(gca,'XScale','log', 'YScale', 'log')
ylim(10.^limy(3, :))
yticks(10.^(limy(3, 1):limy(3, 2)))
end
set(f1, 'position', [200, 200, 700, 450]);