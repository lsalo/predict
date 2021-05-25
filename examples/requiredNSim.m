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
name = {'A', 'B', 'C', 'D', 'E'};           % strati base names.
fname = 'requiredNSim_5manualStrati_100k';  % [] (empty) to not save data.

% Mandatory Input parameters
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
cm = {'kao', 'sme', 'ill', 'mic', 'kao'};       % predominant clay mineral
maxPerm = 5000;                                 % cap max perm? [mD]
rho = 0.7;                                      % Corr. coeff. for multiv. distr.    

% Flow upscaling options
U.useAcceleration = 1;          % requires MEX and AMGCL setup
U.method          = 'mpfa';     % 'tpfa' recommended if useAcc. = 0
U.outflux         = 0;          % compare outflux of fine and upsc. models
U.ARcheck         = 0;          % check if Perm obtained with grid with 
                                % Aspect Ratio of only 3 gives same Perm

% Prepare loop
nStrat = numel(vcl);
nSim = [10 100 500 1000 5000 1*10^4 5*10^4 10^5];
k_md = cell(nStrat, numel(nSim));
tic
for j=1:nStrat           % For each stratigraphy
    disp(['Stratigraphy ' num2str(j) ' / ' num2str(nStrat) '.'])
    fDipIt = faultDip(j);
    
    % FW and HW
    footwall = Stratigraphy(thickness{j}{1}, vcl{j}{1}, 'Dip', dip{j}(1), ...
                            'DepthFaulting', zf{j}(1), ...
                            'DepthBurial', zmax{j}{1}, 'ClayMine', cm{j});
    hangingwall = Stratigraphy(thickness{j}{2}, vcl{j}{2}, 'Dip', dip{j}(2), ...
                               'IsHW', 1, 'NumLayersFW', footwall.NumLayers, ...
                               'DepthFaulting', zf{j}(2), ...
                               'DepthBurial', zmax{j}{2}, 'ClayMine', cm{j});
    
    % Strati in Faulted Section
    mySect = FaultedSection(footwall, hangingwall, fDipIt, ...
                            'maxPerm', maxPerm);
    
    % Get material property distributions
    mySect = mySect.getMatPropDistr();
    
    % Loop over multiple simulation numbers
    for k=1:numel(nSim)
        nSimIt = nSim(k);
        perm = nan(nSimIt, 3);
        
        disp(['Started simulation set ' num2str(k) ' / ' num2str(numel(nSim)) ...
              ' (' num2str(nSimIt) ' realizations)...'])
        parfor n=1:nSimIt    % loop for each realization in a given set
            
            % Generate fault object with properties for each realization
            myFault = Fault(mySect, fDipIt);
            
            % Get dependent variables
            myFault = myFault.getMaterialProperties(mySect, 'corrCoef', rho);
            
            % Generate smear object with T, Tap, L, Lmax
            smear = Smear(mySect, myFault, 1);
            
            % Compute upscaled permeability distribution
            myFault = myFault.upscaleSmearPerm(mySect, smear, U);
            
            % Save result
            perm(n, :) = myFault.Perm;
            
            if mod(n, 5000) == 0
                disp([num2str(n) ' realizations out of ' num2str(nSimIt), ...
                      ' completed.'])
            end
        end
        k_md{j, k} = perm ./ (milli*darcy);
    
        disp(['Simulation ' num2str(k) ' done. ' num2str(nStrat - j) ...
              ' stratigraphies remaining.'])
        disp('-----------------------------------------------------------')
    end
    
    disp(['Stratigraphy ' num2str(j) ' finished.'])
    disp('***************************************************************')
end
telapsed = toc;

% Save data?
if ~isempty(fname)
    disp(['ATTENTION: data saved in: ' pwd ' with filename ' fname])    
    save([fname '.mat']) %,'-v7.3') % larger than 2GB
end


%% Output analysis
% Compare perm distros

% Plot utils
sz = [14, 12];
latx = {'Interpreter','latex'};

% Plot
nbins = 50;
colrs = [0.35 0.35 0.35; 1 0 0; 0 0 1];
plotId = [1 4 6 8];
fh = figure(1);
tiledlayout(numel(plotId), nStrat, 'Padding', 'compact', 'TileSpacing', 'compact');
tileIds = reshape(1:nStrat*numel(plotId), nStrat, numel(plotId))';
for j=1:nStrat
    limx = [floor(log10(min(min(k_md{j, end})))), ...
            ceil(log10(max(max(k_md{j, end}))))];
    edges = logspace(limx(1), limx(2), nbins);
    probs = [histcounts(k_md{j, end}(:, 1), edges); ...
             histcounts(k_md{j, end}(:, 2), edges); ...
             histcounts(k_md{j, end}(:, 3), edges)]' ./ nSim(end);
    limy = floor(log10(min(probs(probs > 0))));
    for k=1:numel(plotId)
        nexttile(tileIds(k, j))
        hold on
        histogram(k_md{j, plotId(k)}(:, 1), edges, 'Normalization', 'probability', ...
            'FaceColor', colrs(1, :), 'EdgeColor', colrs(1,:), ...
            'FaceAlpha', 1, 'DisplayName', '$k_{xx}$')
        histogram(k_md{j, plotId(k)}(:, 2), edges, 'Normalization', 'probability', ...
            'FaceColor', colrs(2, :), 'EdgeColor', 'none', ...
            'FaceAlpha', .7, 'DisplayName', '$k_{yy}$')
        histogram(k_md{j, plotId(k)}(:, 3), edges, 'Normalization', 'probability', ...
            'FaceColor', colrs(3, :), 'EdgeColor', 'none', ...
            'FaceAlpha', .5, 'DisplayName', '$k_{zz}$')
        if k == 1
            xlabel('$k$ [mD]', latx{:}, 'fontSize', sz(2))
            ylabel('P [-]', latx{:}, 'fontSize', sz(2))
            title([name{j} ', $N_\mathrm{sim}$ = ' num2str(nSim(plotId(k)))], latx{:}, ...
                   'fontSize', sz(2))
            if j == 1
                leg = legend(latx{:}, 'fontSize', sz(2), 'location', 'northwest');
                set(leg.BoxFace, 'ColorType','truecoloralpha', ...
                    'ColorData', uint8(255*[1;1;1;.6]));
            end
        else
            if nSim(plotId(k)) < 10000
                title(num2str(nSim(plotId(k))), latx{:}, 'fontSize', sz(2))
            else
                title(num2str(nSim(plotId(k)), 3), latx{:}, 'fontSize', sz(2))
            end
        end
        xlim([10^limx(1, 1) 10^limx(1, 2)])
        xticks(10.^(limx(1, 1):2:limx(1, 2)));
        %ylim([10^limy 1])
        %yticks(10.^(limy:0))
        ylim([0 0.5])
        yticks(0:.1:.5)
        grid on
        %set(gca,'XScale','log','YScale','log')
        set(gca, 'XScale', 'log');
    end    
end
set(fh, 'position', [200, 200, 900, 600]);