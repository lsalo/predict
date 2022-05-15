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
mrstModule add mrst-gui coarsegrid upscaling incomp mpfa
mrstVerbose on

%% 2. Define Model and Upscale Permeability

% 2.1 Mandatory input parameters
% Footwall first and hangingwall next, e.g. {[footwall, FW], [hangingwall, HW]}. 
% We need to define the layer thickness, clay content, layer dip angle, fault 
% dip angle, faulting depth, and burial depth. Further details about input parameter 
% formatting, etc can always be checked from the documentation in the classes 
% and functions.
thickness = {repelem(25, 1, 4), [50 50]};                % [m]
vcl       = {[0.1 0.7 0.2 0.8], ...
             [0.3 0.5]};                        % fraction [-]
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
U.method          = 'mpfa';     % 'tpfa' recommended if useAcceleration = 0
U.outflux         = 0;          % compare outflux of fine and upscaled model
U.ARcheck         = 0;          % check if Perm obtained with grid with aspect ratio of 
                                % only 5 gives same output.
Nsim              = 100;        % Number of 3D simulations/realizations

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

% 2.6 Get base grid
% We generate a base grid with arbitrary thickness, to be modified at each
% realization (much faster than generating n grids from scratch). Note that
% the 3D grid is an extruded grid, where the x dimension is along-strike.
D = sum(mySect.Tap(mySect.FW.Id));
L  = mySect.MatPropDistr.length.fcn(D);    %  equal to disp for now
T0 = 1;
disp('Constructing initial grid...')
if dim == 2,        G0 = makeFaultGrid(T0, D);
elseif dim == 3,    G0 = makeFaultGrid(T0, D, L);
end

% 2.7 Generate intermediate variable samples, calculate smear dimensions and upscale permeability
% We create two container variables (faults and smears) where we'll save all 
% data for each realization. For each realization, the code defines a Fault object, 
% generates intermediate variable samples, calculates the smear dimensions, and, 
% within upscaleSmearPerm, generates a fault material distribution consistent 
% with the inputs and upscales the permeability.
% Generate fault object with properties for each realization
faults = cell(Nsim, 1);
smears = cell(Nsim, 1);
tstart = tic;
assert(dim==3);
parfor n=1:Nsim    % parfor allowed if you have the parallel computing toolbox
    
    % TBD: get segmentation for this realization
    segLen = [5 10 5 15 10 20 5 10 20]; 
    nSeg = numel(segLen);
    
    extrudedPerm = zeros(G.cells.num, 6);   %[kxx, kxy, kxz, kyy, kyz, kzz]
    for k=1:nSeg
        myFaultSection = Fault(mySect, faultDip, dim);
        
        % Get material property (intermediate variable) samples, and fix
        % along-strike thickness of current realization if 3D.
        myFaultSection = myFaultSection.getMaterialProperties(mySect, 'corrCoef', rho);
        if k==1
            thick3D = myFaultSection.MatProps.thick;
            % Update grid dimensions with sampled fault thickness
            G = updateGrid(G0, myFaultSection.MatProps.thick);
        else
            myFaultSection.MatProps.thick = thick3D;
        end
        
        % Generate smear object with T, Tap, L, Lmax
        smear = Smear(mySect, myFaultSection, G, 1);
        
        % Place fault materials and assign cell-based properties in 2D
        % section
        myFaultSection = myFaultSection.placeMaterials(mySect, smear, G);
        
        % Extrude 2D section to fill current segment
        % TBD: we want to know material type, issmear, etc for each cell in
        %      extruded grid, so change extrudedPerm for extrudedVals.
        extrudedPerm = assignExtrudedPerm(extrudedPerm, myFaultSection, ...
                                          segLen(k), G.cellDim(1));
        
        % Save results
        faults{n, k} = myFaultSection;
        smears{n, k} = smear;
    end
    
    % Compute 3D upscaled permeability distribution
    myFault = myFaultSection.upscaleProps(G, U);
    
    if mod(n, 100) == 0
        D(['Simulation ' num2str(n) ' / ' num2str(Nsim) ' completed.'])
    end
end
telapsed = toc(tstart);

%% 3. Output Analysis
% 3.1 Visualize stratigraphy and fault (with thickness corresponding to 1st realization)
mySect.plotStrati(faults{1}.MatProps.thick, faultDip);  

% 3.2 Visualize intermediate variables
% We define a given parent material (id from 1 to n of materials in stratigraphy), 
% and generate histograms and correlation matrix plots.
layerId = 4;                                            
plotMatPropsHist(faults, smears, mySect, layerId) 
% MatProps correlations
[R, P] = plotMatPropsCorr(faults, mySect, layerId);

% 3.3 Visualize fault materials
% Visualization for one realization. Choice can be 'randm' (random), 'maxX' 
% (realization with maximum upscaled permeability in across the fault), 'minX', 
% 'maxZ' or 'minZ'.
% General fault materials and perm view
plotId = selectSimId('randm', faults, Nsim);                % simulation index
faults{plotId}.plotMaterials(mySect, G0) 

% 3.4. Visualize upscaled permeability
% Plot upscaled permeability distributions (all simulations)
plotUpscaledPerm(faults)