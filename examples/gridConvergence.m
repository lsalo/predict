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
thickness = {{repelem(10, 1, 10), repelem(10, 1, 10)}, ...
             {[25 25 25 25], [25 25 25 25]}, ...
             {[50 50], [50 50]}, ...
             {[5 10 15 10 20 10 10 5 15], [20 10 20 10 30 10]}, ...
             {[20 30 30 20]; [40 20 10 10 10 10]}};
vcl       = {{repmat([0.2 0.6], 1, 5), repmat([0.5 0.3], 1, 5)}, ...
             {[0.8 0.3 0.5 0], [0.3, 0.7, 0.15, 0.6]}, ...
             {[0.5 0.1], [0.5 0.1]}, ...
             {[0.3 0.6 0.1 0.7 0.2 0.8 0.3 0.9 0.1], [0.2, 0.7, 0.25, 0.8, 0.3, 0.9]}, ...
             {[0 0.4 0.3 0.6], [0.1 0.4 0.2 0.6 0.1 0.5]}};
dip       = {[0, 0], [0, 0], [10, 20], [5, -5], [0, 0]};
faultDip  = [50, 60, 75, 85, 55];

% Optional Input parameters
nl   = [numel(vcl{1}{1}), numel(vcl{1}{2}); ...
        numel(vcl{2}{1}), numel(vcl{2}{2}); ...
        numel(vcl{3}{1}), numel(vcl{3}{2}); ...
        numel(vcl{4}{1}), numel(vcl{4}{2}); ...
        numel(vcl{5}{1}), numel(vcl{5}{2})];  % just for convenience here
zf   = {[100, 100], [500, 500], [1000, 1000], [100, 100], [2000, 2000]};    % m
zmax = {{repelem(800, 1, nl(1,1)), repelem(800, 1, nl(1,2))}, ...
        {repelem(1000, 1, nl(2,1)), repelem(1000, 1, nl(2,2))}, ...
        {repelem(1000, 1, nl(3,1)); repelem(1000, 1, nl(3,2))}, ...
        {repelem(2000, 1, nl(4,1)); repelem(2000, 1, nl(4,2))}, ...
        {repelem(3000, 1, nl(5,1)); repelem(3000, 1, nl(5,2))}};
cm = {'kao', 'sme', 'ill', 'mic', 'kao'};   % predominant clay mineral
maxPerm = 5000;                             % cap max perm? [mD]
rho = 0.7;                                  % Corr. coeff. for multiv. distr.

% Flow upscaling options
U.useAcceleration = 1;          % 1 requires MEX and AMGCL setup
U.method          = 'mpfa';     % 'tpfa' recommended if useAcc. = 0
U.outflux         = 0;          % compare outflux of fine and upsc. models
U.ARcheck         = 0;          % check if Perm obtained with grid with 
                                % Aspect Ratio of only 3 gives same Perm
                                
% Initialize target grid cell sizes and result containers
Nstrat = numel(vcl);
gridRes = repmat([0.5, 5; 0.2, 2; 0.1, 1; 0.05, 0.5; ...
                  0.02, 0.2; 0.01, 0.1], 1, 1, Nstrat);         
Ngrid = size(gridRes, 1);
poro  = zeros(Ngrid, Nstrat);
perm  = zeros(Ngrid, 3, Nstrat);
cellDim = zeros(Ngrid, 2, Nstrat);

tic
parfor k=1:Nstrat
faultDipIt = faultDip(k);
% FW and HW
footwall = Stratigraphy(thickness{k}{1}, vcl{k}{1}, ...
                        'Dip', dip{k}(1), 'DepthFaulting', zf{k}(1), ...
                        'DepthBurial', zmax{k}{1}, 'ClayMine', cm{k});
hangingwall = Stratigraphy(thickness{k}{2}, vcl{k}{2}, ...
                           'Dip', dip{k}(2), 'IsHW', 1, ...
                           'NumLayersFW', footwall.NumLayers, ...
                           'DepthFaulting', zf{k}(2), ...
                           'DepthBurial', zmax{k}{2}, 'ClayMine', cm{k});

% Strati in Faulted Section
mySect = FaultedSection(footwall, hangingwall, faultDipIt, 'maxPerm', maxPerm);

% Get material property distributions
mySect = mySect.getMatPropDistr();

    % Run grid loop for each strati case
    for n=1:Ngrid
    % Reset random seed (same values in MatProps and smear placement for
    % successive grids with different sizes)
    rng(10)

    % Generate fault object with properties for each realization
    myFault = Fault(mySect, faultDipIt, gridRes(n, :, k));                 

    % Get dependent variables
    myFault = myFault.getMaterialProperties(mySect, 'corrCoef', rho);   

    % Generate smear object with T, Tap, L, Lmax
    smear = Smear(mySect, myFault, 1);

    % Compute upscaled permeability distribution
    myFault = myFault.upscaleSmearPerm(mySect, smear, U);
    
    % Save result
    poro(n, k) = myFault.Poro;
    perm(n, :, k) = myFault.Perm;
    cellDim(n, :, k) = myFault.Grid.cellDim;
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
hL = cellDim(:,2,1);
A  = cellDim(:,1,1).*cellDim(:,2,1);
k_md = perm./(milli*darcy);

% Plotting utilities
sz = [14, 12];
latx = {'Interpreter','latex'};
%colrs = [128 0 0; 0 130 200; 255 225 128; 0 0 128; 0 0 0] ./ 255;
idRes = 3;

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
plot(1/hL(idRes), k_md(idRes,1,k), '-o', 'color', 'k', 'MarkerFaceColor', ...
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
plot(1/hL(idRes), k_md(idRes,2,k), '-o', 'color', 'r', ...
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
plot(1/hL(idRes), k_md(idRes,3,k), '-o', 'color', 'b', ...
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