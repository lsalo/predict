%% Example 0: Single stratigraphic case + analysis
% This is a complete introductory example. It shows how to load the appropriate 
% MRST modules, define the inputs according to a given faulted stratigraphy, and 
% generate the output permeability distributions. A comprehensive analysis of 
% the results is also shown.
% 
% We first make sure that the workspace is clean:
clear
close all force

%% 1. Load Required MRST Modules
% First, navigate to the mrst folder and run |startup.m|. We can then load the 
% appropriate modules for generating MRST grids and upscale the permeability:
mrstModule add mrst-gui coarsegrid upscaling incomp mpfa mimetic
mrstVerbose off     % set to on for more insight in the command window

%% 2. Define Model and Upscale Permeability

% 2.1 Mandatory input parameters
% Footwall first and hangingwall next, e.g. {[footwall, FW], [hangingwall, HW]}. 
% We need to define the layer thickness, clay content, layer dip angle, fault 
% dip angle, faulting depth, and burial depth. Further details about input parameter 
% formatting, etc can always be checked from the documentation in the classes 
% and functions.
thickness = {repelem(25, 1, 4), [20 30 25 25]};                % [m]
vcl       = {[0.1 0.4 0.2 0.5], ...
             [0.3 0.6 0.1 0.7]};                        % fraction [-]
dip       = [0, -5];                                                        % [deg.]
faultDip  = 70;                                                             % [deg.]
zf        = [500, 500];                                                     % [FW, HW], [m]
zmax      = {repelem(2000, numel(vcl{1})), repelem(2000, numel(vcl{2}))};   % {FW, HW}
dim       = 3;                    % dimensions (2 = 2D, 3 = 3D)

% 2.2 Optional input parameters
% In this case, we indicate a maximum fault material permeability of and a correlation 
% coefficient for dependent variables:
maxPerm = 1000;                 % [mD]
rho     = 0.6;                  % Corr. coeff. for multivariate distributions

% 2.3 Flow upscaling options and number of simulations
U.useAcceleration = 1;          % 1 requires MEX setup, 0 otherwise (slower for MPFA).
U.method          = 'tpfa';     % 'tpfa' recommended for 3D
U.coarseDims      = [1 1 1];    % number of cells [x, y, z] in coarse grid
Nsim              = 1000;       % Number of 3D simulations/realizations

% 2.4 Define Stratigraphy and FaultedSection objects
% Organize the input parameters in HW and FW, and use that info to create a 
% FaultedSection object which contains all required information.
% FW and HW
footwall = Stratigraphy(thickness{1}, vcl{1}, 'Dip', dip(1), ...
                        'DepthFaulting', zf(1), 'DepthBurial', zmax{1});
hangingwall = Stratigraphy(thickness{2}, vcl{2}, 'Dip', dip(2), 'IsHW', 1, ...
                           'NumLayersFW', footwall.NumLayers, ...
                           'DepthFaulting', zf(2), 'DepthBurial', zmax{2});

% Instantiate FaultedSection object (Strati in Faulted Section)
mySect = FaultedSection(footwall, hangingwall, faultDip, 'maxPerm', maxPerm);

% 2.5 Get material distributions
% We use the inputs to constrain the ranges and distributions for each of the 
% intermediate variables.
% Get material property distributions
mySect = mySect.getMatPropDistr();

% 2.6 Generate intermediate variable samples, calculate smear dimensions 
%     and upscale permeability.
% We create two container variables (faults and smears) where we'll save all 
% data for each realization. For each realization, the code defines a Fault object, 
% generates intermediate variable samples, calculates the smear dimensions, and, 
% within upscaleSmearPerm, generates a fault material distribution consistent 
% with the inputs and upscales the permeability.
% Generate fault object with properties for each realization
assert(dim==3);
faultSections = cell(Nsim, 1);
smears = cell(Nsim, 1);
faults = cell(Nsim, 1);
upscaledPerm = zeros(Nsim, 3);
D = sum(mySect.Tap(mySect.FW.Id));      % displacement
tstart = tic;
parfor n=1:Nsim    % parfor allowed if you have the parallel computing toolbox
    % Instantiate fault section and get segmentation for this realization
    myFaultSection = Fault(mySect, faultDip);
    myFault = ExtrudedFault(myFaultSection, mySect);
    myFault = myFault.getSegmentationLength(mySect, 4);
    G = [];
    for k=1:numel(myFault.SegLen)
        % Get material property (intermediate variable) samples, and fix
        % fault thickness of current realization (3D only).
        myFaultSection = myFaultSection.getMaterialProperties(mySect, 'corrCoef', rho);
        myFaultSection.MatProps.thick = myFault.Thick;
        if isempty(G)
            G = makeFaultGrid(myFault.Thick, myFault.Disp, ...
                              myFault.Length, myFault.SegLen);
        end
        
        % Generate smear object with T, Tap, L, Lmax
        smear = Smear(mySect, myFaultSection, G, 1);
        
        % Place fault materials and assign cell-based properties in 2D section
        myFaultSection = myFaultSection.placeMaterials(mySect, smear, G);
        
        % Extrude 2D section to fill current segment
        myFault = myFault.assignExtrudedVals(G, myFaultSection, k);
        
        % Save results
        faultSections{n}{k} = myFaultSection;
        smears{n}{k} = smear;
    end
    
    % Compute 3D upscaled permeability distribution
    myFault = myFault.upscaleProps(G, U);
    
    % Save results
    faults{n} = myFault;
    if mod(n, 50) == 0
        disp(['Simulation ' num2str(n) ' / ' num2str(Nsim) ' completed.'])
    end
end
telapsed = toc(tstart);

%% 3. Output Analysis
% 3.1 Visualize stratigraphy and fault (with thickness corresponding to 1st realization)
%mySect.plotStrati(faults{1}.MatProps.thick, faultDip);  

% % 3.2 Visualize intermediate variables
% % We define a given parent material (id from 1 to n of materials in stratigraphy), 
% % and generate histograms and correlation matrix plots.
% layerId = 4;                                            
% plotMatPropsHist(faults, smears, mySect, layerId) 
% % MatProps correlations
% [R, P] = plotMatPropsCorr(faults, mySect, layerId);
% 
% % 3.3 Visualize fault materials
% % Visualization for one realization. Choice can be 'randm' (random), 'maxX' 
% % (realization with maximum upscaled permeability in across the fault), 'minX', 
% % 'maxZ' or 'minZ'.
% % General fault materials and perm view
% plotId = selectSimId('randm', faults, Nsim);                % simulation index
% faults{plotId}.plotMaterials(mySect, G0) 
% 
% % 3.4. Visualize upscaled permeability
% % Plot upscaled permeability distributions (all simulations)
% plotUpscaledPerm(faults)

% Histograms
latx = {'Interpreter', 'latex'};
sz = [14, 12];
K = log10(upscaledPerm/(milli*darcy));
nbins = 25;
logMinP = min(min(K));
logMaxP = max(max(K));
edges = linspace(fix(logMinP)-1, fix(logMaxP)+1, nbins);

fh = figure(5);
tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
labls = ["$\log_{10}(k_{xx}$ [mD])", ...
    "$\log_{10}(k_{yy}$ [mD])", ...
    "$\log_{10}(k_{zz}$ [mD])"];
nexttile(1)
histogram(K(:, 1), edges, 'Normalization', 'probability', ...
    'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 1)
xlabel(labls(1), latx{:}, 'fontSize', sz(2))
ylabel('P [-]', latx{:}, 'fontSize', sz(2))
xlim([fix(logMinP)-1 fix(logMaxP)+1])
ylim([0 0.6]); yticks(0:.2:.6)
grid on
%xticks(10.^(fix(logMinP)-1:2:fix(logMaxP)+1))

nexttile(2)
rr = [255, 125, 125]/255;
histogram(K(:, 2), edges, 'Normalization', 'probability', ...
    'FaceColor', rr, 'FaceAlpha', 1)
xlabel(labls(2), latx{:}, 'fontSize', sz(2))
%ylabel('P [-]', latx{:}, 'fontSize', sz(2))
xlim([fix(logMinP)-1 fix(logMaxP)+1])
ylim([0 0.6]); yticks(0:.2:.6)
grid on
%xticks(10.^(fix(logMinP)-1:2:fix(logMaxP)+1))
hold off

nexttile(3)
bb = [125, 125, 255]/255;
histogram(K(:, 3), edges, 'Normalization', 'probability', ...
    'FaceColor', bb, 'FaceAlpha', 1)
xlabel(labls(3), latx{:}, 'fontSize', sz(2))
%ylabel('P [-]', latx{:}, 'fontSize', sz(2))
xlim([fix(logMinP)-1 fix(logMaxP)+1])
ylim([0 0.6]); yticks(0:.2:.6)
grid on
%xticks(10.^(fix(logMinP)-1:2:fix(logMaxP)+1))
set(fh, 'position', [200, 200, 150, 350])