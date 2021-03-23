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
U.useAcceleration = 1;          % requires MEX and AMGCL setup
U.method          = 'mpfa';     % 'tpfa' recommended if useAcc. = 0
U.outflux         = 0;          % compare outflux of fine and upsc. models
U.ARcheck         = 0;          % check if Perm obtained with grid with 
                                % Aspect Ratio of only 3 gives same Perm
                                
% Prepare
Nstrat = numel(vcl);
ar = [NaN, 50, 25]; %10 5 2]
Ngrid = numel(ar);
perm = zeros(Ngrid, 3, Nstrat);
cellDim = zeros(Ngrid, 2, Nstrat);
tic
for k=1:Nstrat
    disp(['Stratigraphic case ' num2str(k) ' / ' num2str(Nstrat)])
    
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
    
    % Generate fault object with properties for each realization
    myFault = Fault(mySect, faultDip(k));
    
    % Get dependent variables
    myFault = myFault.getMaterialProperties(mySect, 'maxPerm', maxPerm, ...
                                            'siltInClay', siltInClay(k));
    
    % Generate smear object with T, Tap, L, Lmax
    Tap = getApparentThick(mySect, myFault.Dip);
    smear = Smear(mySect.Vcl, mySect.IsClayVcl, mySect.Thick, Tap, ...
                  mySect.DepthFaulting, myFault, 1, mySect);
    
    % Compute upscaled permeability distribution
    myFault = myFault.upscaleSmearPerm(mySect, smear, U);
    
    for n=1:Ngrid
        if n == 1
            perm(n, :, k) = myFault.Perm;
            cellDim(n, :, k) = myFault.Grid.CellDim;
    
            % Base grid
            G = makeFaultGrid(myFault.MatProps.Thick, myFault.Disp, ...
                              myFault.Grid.TargetCellDim);
            ar(n) = G.CellDim(2)/G.CellDim(1);
            L = max(G.faces.centroids) - min(G.faces.centroids);
            Dp{1} = 5*barsa;

        elseif ar(n) < ar(1)
            % Iteration grid
            Lz = G.cartDims(2)*G.CellDim(2);
            nels = [G.cartDims(1) round(Lz/(ar(n)*G.CellDim(1)))];
            G2 = computeGeometry(cartGrid([nels(1), nels(2)], ...
                                          [myFault.MatProps.Thick, myFault.Disp]));
            disp(['--> Iteration grid with cell aspect ratio = ' num2str(ar(n)) ...
                  ' created.'])
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
            rock2.perm = myFault.Grid.Perm(mapG2toG, :);
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

            % Save result
            perm2 = [perm2(1), perm2(2), perm2(4)];
            perm(n, :, k) = perm2; clear perm2
            cellDim(n, :, k) = [myFault.MatProps.Thick/nels(1), ...
                                myFault.Disp/nels(2)];
        end
        
        disp(['Simulation ' num2str(n) 'out of ' num2str(Ngrid) ' completed.'])
    end
    
    disp('---------------------------------------------------------------')
end
toc


%% Output analysis

% Plot Perm in each dir. vs aspect ratio
