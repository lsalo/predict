%% Analysis of cell aspect ratio for flow-based upscaling
%
% This analysis is very RAM intensive, so start by modifying the 'arTarget'
% variable to a minimum of 25 and then diminish progressively based on
% available RAM.
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
U.useAcceleration = 1;          % requires MEX and AMGCL setup
U.method          = 'mpfa';     % 'tpfa' recommended if useAcc. = 0
U.outflux         = 0;          % compare outflux of fine and upsc. models
U.ARcheck         = 0;          % check if Perm obtained with grid with 
                                % Aspect Ratio of only 3 gives same Perm
                                
% Prepare
Nstrat = numel(vcl);
arTarget = [NaN, 50, 25, 10, 5];
ar = zeros(Nstrat, numel(arTarget));
Ngrid = numel(arTarget);
perm = nan(Ngrid, 2, Nstrat);
cellDim = zeros(Ngrid, 2, Nstrat);
tic
for k=1:Nstrat
    disp(['Stratigraphic case ' num2str(k) ' / ' num2str(Nstrat)])
    
    % FW and HW
    footwall = Stratigraphy(thickness{k}{1}, vcl{k}{1}, 'Dip', dip{k}(1), ...
                            'DepthFaulting', zf{k}(1), ...
                            'DepthBurial', zmax{k}{1}, 'ClayMine', cm{k});
    hangingwall = Stratigraphy(thickness{k}{2}, vcl{k}{2}, 'Dip', dip{k}(2), ...
                               'IsHW', 1, 'NumLayersFW', footwall.NumLayers, ...
                               'DepthFaulting', zf{k}(2), ...
                               'DepthBurial', zmax{k}{2}, 'ClayMine', cm{k});
    
    % Strati in Faulted Section
    mySect = FaultedSection(footwall, hangingwall, faultDip(k), ...
                            'maxPerm', maxPerm);
    
    % Get material property distributions
    mySect = mySect.getMatPropDistr();
    
    % Generate fault object with properties for each realization
    myFault = Fault(mySect, faultDip(k));                 

    % Get dependent variables
    myFault = myFault.getMaterialProperties(mySect, 'corrCoef', rho); 
    
    itnum = 0;
    while myFault.MatProps.thick < myFault.Disp/500  % avoid thinest faults
        if itnum == 0
            disp('Material properties recalculated.')
        end
        myFault = myFault.getMaterialProperties(mySect, 'corrCoef', rho);
        itnum = itnum + 1;
        disp(['iteration number = ' num2str(itnum)])
    end
    
    % Generate smear object
    smear = Smear(mySect, myFault, 1);
    
    % Compute upscaled permeability distribution
    myFault = myFault.upscaleSmearPerm(mySect, smear, U);
    
    for n=1:Ngrid
        if n == 1
            perm(n, :, k) = [myFault.Perm(1) myFault.Perm(end)];
            cellDim(n, :, k) = myFault.Grid.cellDim;
    
            % Base grid
            G = makeFaultGrid(myFault.MatProps.thick, myFault.Disp, ...
                              myFault.Grid.targetCellDim);
            ar(k, n) = G.CellDim(2)/G.CellDim(1);
            L = max(G.faces.centroids) - min(G.faces.centroids);
            Dp{1} = 5*barsa;
            disp(['Base grid aspect ratio = ' num2str(ar(k, n))])

        elseif arTarget(n) < ar(k, 1)
            % Iteration grid
            Lz = G.cartDims(2)*G.CellDim(2);
            nels = [G.cartDims(1) round(Lz/(arTarget(n)*G.CellDim(1)))];
            cellDim(n, :, k) = [myFault.MatProps.thick/nels(1), ...
                                myFault.Disp/nels(2)];
            ar(k, n) = cellDim(n, end, k) / cellDim(n, 1, k);
            G2 = computeGeometry(cartGrid([nels(1), nels(2)], ...
                                          [myFault.MatProps.thick, myFault.Disp]));
            disp(['--> Iteration grid with cell aspect ratio = ' ...
                  num2str(round(ar(k, n), 2)) ' created.'])
            disp(['    This grid has a total of ' num2str(G2.cells.num) ...
                  ' cells. Now computing cell distances...'])

            % assign perm based on cell centroid distance to closest cell
            distx = pdist2(G2.cells.centroids(:,1), G.cells.centroids(:,1));
            distz = pdist2(G2.cells.centroids(:,end), G.cells.centroids(:,end));
            dist = sqrt(distx.^2 + distz.^2);
            disp(['    Distances computed. Now calculating MPFA perm with ' ...
                  'finer z mesh...'])

            % Compute 2nd perm
            [~, mapG2toG] = min(dist,[],2);
            rock2.perm = myFault.Grid.perm(mapG2toG, :);
            hTmp2 = computeMultiPointTrans(G2, rock2, 'invertBlocks', 'mex');
            psolver = @(state0, G, fluid, bc) incompMPFA(state0, G, hTmp2, ...
                                                         fluid, 'bc', bc);
            fluid = initSingleFluid('mu' , 1, 'rho', 1);
            perm2 = diag(myupscalePermeabilityFixed(G2, Dp{1}, psolver, fluid, L));
            disp(['    Finer z mesh perm [x, z], in m^2, is [' ...
                  num2str(perm2(1)) ' ' num2str(perm2(end)) ']'])
            disp(['    Ratio Perm(finerGrid) / Perm(originalGrid) [x, z] is [' ...
                  num2str(perm2(1)/myFault.Perm(1)) ' ' ...
                  num2str(perm2(end)/myFault.Perm(end)) ']'])

            % Save perm result
            perm2 = [perm2(1), perm2(4)];
            perm(n, :, k) = perm2; clear perm2
        end
        
        
        disp(['    Simulation ' num2str(n) ' out of ' num2str(Ngrid) ...
              ' completed.'])
    end
    disp('---------------------------------------------------------------')
end
toc


%% Output analysis

% Cell permeability assignment check
% f99 = figure(99);
% tiledlayout(1, 4, 'Padding', 'compact', 'TileSpacing', 'compact');
% nexttile
% plotCellData(G, log10(myFault.Grid.Perm(:,end)./(milli*darcy)), 'facea', 1, 'edgea', 0);
% colormap(copper); title('G, Z Perm')
% xlim([0 myFault.MatProps.Thick]); ylim([0 myFault.Disp]); colorbar;
% nexttile
% plotCellData(G2, log10(rock2.perm(:,end)./(milli*darcy)), 'facea', 1, 'edgea', 0);
% colormap(copper); title('G2, Z Perm')
% xlim([0 myFault.MatProps.Thick]); ylim([0 myFault.Disp]); colorbar;
% nexttile
% plotCellData(G, log10(myFault.Grid.Perm(:,end)./(milli*darcy)), 'facea', 1, 'edgea', 0.2);
% colormap(copper); title('G, Z Perm')
% xlim([0 myFault.MatProps.Thick]); axis equal; ylim([20 25]);
% nexttile
% plotCellData(G2, log10(rock2.perm(:,end)./(milli*darcy)), 'facea', 1, 'edgea', 0.2);
% colormap(copper); title('G2, Z Perm')
% xlim([0 myFault.MatProps.Thick]); axis equal; ylim([20 25]); 
% set(f99, 'position', [300, 300, 800, 600]);

% Plotting utilities
sz = [14, 12];
latx = {'Interpreter','latex'};
%colrs = [128 0 0; 0 130 200; 255 225 128; 0 0 128; 0 0 0] ./ 255;
markrs = ['o', 's', 'd', 'p', 'h'];
limx = [1 max([max(max(ar)), 100])];

f1 = figure(1);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
for k=1:Nstrat
k_md = perm(:,:,k)./(milli*darcy);
idnan = isnan(k_md(:,1));
k_md(idnan, :) = [];

nexttile(1)
hold on
plot(ar(k, ~idnan), k_md(:,1)/k_md(1,1), '-', 'marker', markrs(k), 'color', ...
     'k', 'MarkerSize', 4, 'DisplayName', name{k})
plot(ar(k, 1), 1, '-', 'marker', markrs(k), 'color', 'k', 'MarkerFaceColor', ...
     [0.5 0.5 0.5], 'MarkerSize', 6, 'HandleVisibility','off')
hold off
if k == Nstrat
    grid on
    xlabel('AR [-]', latx{:}, 'fontSize', 12)
    ylabel('$k_{xx}^{\mathrm{ref.}} / k_{xx}^{\mathrm{base}}$ [-]', latx{:}, 'fontSize', 12)
    xlim(limx)
    xticks([1 10 100])
    xticklabels({'1' '10' '100'})
    set(gca,'XScale','log')
    ylim([1 1.4])
    yticks(1:.1:1.4)
    leg = legend(latx{:}, 'fontSize', sz(2), 'location', 'northeast');
    set(leg.BoxFace, 'ColorType','truecoloralpha', ...
        'ColorData', uint8(255*[1;1;1;.6])); 
end

nexttile(2)
hold on
plot(ar(k, ~idnan), k_md(:,2)/k_md(1,2), '-', 'marker', markrs(k), 'color', ...
     'b', 'MarkerSize', 4)
plot(ar(k, 1), 1, '-', 'marker', markrs(k), 'color', 'b', ...
     'MarkerFaceColor', [125, 125, 255]/255, 'MarkerSize', 6)
hold off
if k == Nstrat
    grid on
    %xlabel('$1/h_\mathrm{L}$ [m$^{-1}$]', latx{:}, 'fontSize', 12)
    ylabel('$k_{zz}^{\mathrm{ref.}} / k_{zz}^{\mathrm{base}}$ [-]', latx{:}, 'fontSize', 12)
    xlim(limx)
    xticks([1 10 100])
    xticklabels({'1' '10' '100'})
    set(gca,'XScale','log')
    ylim([1 1.4])
    yticks(1:.1:1.4)
end
end
set(f1, 'position', [500, 500, 500, 250]);