classdef Fault
    %
    % SUMMARY:
    %   Define a FaultedSection object, with corresponding properties,
    %   based on the footwall and hangingwall stratigraphy. It does not
    %   depend on fault characteristics, just Stratigraphy objects.
    %
    %
    % SYNOPSIS:
    %   mySection = FaultedSection(footwall, hangingwall)
    %
    %
    % DESCRIPTION:
    % 
    %   TBD
    %  
    %   Decreasing porosity and perm with increased shear strain is
    %   consistent with literature, though more pronounced for
    %   low-clay content and measured at high confinement (Crawford
    %   et al., JGR, 2008).
    % 
    % REQUIRED PARAMETERS:
    %   footwall:       Stratigraphy object with corresponding property 
    %                   values for the footwall of the fault.
    %
    %   hangingwall:    Stratigraphy object with corresponding property 
    %                   values for the footwall of the fault.
    %
    %
    % OPTIONAL PARAMETERS:
    %   'property' - Set property to the specified value.
    %
    %
    % RETURNS:
    %   Class instance.
    %   
    % ____________________________________________________________________
    
    
    properties
        Dip                             % Dip angle of the fault [deg.]    
        MatProps                        % Material props (phi, SSFc, k, k')
        MatMap                          % Material mapping to each grid cell
    end
    
    properties (SetAccess = protected)
       Disp                             % Fault displacement [m]
       Poro                             % Final upscaled Poro (1, 1)
       Perm                             % Final upscaled Perm (1, 3)
       Vcl                              % Final upscaled Vcl (1, 1)
       Grid                             % contains resolution, cell Poro and Perm (G.cells.num, Nsim)
    end
    
    properties (SetAccess = protected, Dependent)
       Throw                            % Fault throw [m]
    end
    
    properties (SetAccess = protected, Hidden)
       Delta    % dip of clay smears within fault
       Alpha    % Small angle between the fault dip direction & the smear dip
       Zeta     % Angle between a line perpendicular to the smear dip 
                % direction and a diagonal from the top left or top right 
                % of fault to bottom right or bottom left of fault.
    end
    
    methods
        function obj = Fault(FS, dip)
            %
            % We instantiate a Fault object with the fundamental 
            % dimensions.
            %
            
            % Geometry
            obj.Dip    = dip;
            obj.Disp   = sum(FS.Tap(FS.FW.Id));   
        end
        
        function throw = get.Throw(obj)
            %
            % Get vertical throw
            %
            throw = obj.Disp*sind(obj.Dip); 
        end
        
        function obj = getMaterialProperties(obj, FS, varargin)
            % Get material parameters
            %
            % INPUT:
            %   obj: An instance of Fault
            %   FS:  An instance of FaultedSection
            %
            
            % Optional inputs
            opt.corrCoef = 0.6;                   % correlation coefficient
            opt = merge_options_relaxed(opt, varargin{:});
            
            % Convenience
            idc  = FS.Vcl >= FS.IsClayVcl;
            zf   = FS.DepthFaulting;
            
            
            % -------------------------------------------------------------
            % Fault thickness, poro, and Perm across (correlated)
            % -------------------------------------------------------------
            n = 1;
            p = opt.corrCoef;
            if opt.corrCoef > 0
                U = copularnd('Gaussian', [1 p p; p 1 p; p p 1], n);                
            else 
                U = copularnd('Gaussian', [1 -p -p; -p 1 -p; -p -p 1], n); 
            end
            
            % Assert that marginal distros are all beta/unif
            assert(strcmp(FS.MatPropDistr.thick.type, 'beta'));
            %iskunif = cellfun(@(x) strcmp(x, 'unif'), ...
            %                  FS.MatPropDistr.perm.type);
                           
            % Transform bivariate samples to marginals and ranges
            % Fault thickness
            X1 = betainv(U(1), FS.MatPropDistr.thick.param(1), ...
                         FS.MatPropDistr.thick.param(2));
            thickRangeLog = log10(FS.MatPropDistr.thick.range(obj.Disp));
            obj.MatProps.thick = 10^(thickRangeLog(1) + ...
                                     diff(thickRangeLog) .* X1);
            
            % Sand Perm and Poro
            sandPermRange = cell2mat(FS.MatPropDistr.perm.range(~idc)');
            obj.MatProps.permRange(:, ~idc) = sandPermRange';
            sandPermRangeLog = log10(sandPermRange);
            obj.MatProps.perm(~idc) = 10.^(sandPermRangeLog(:, 1) + ...
                                       diff(sandPermRangeLog, [], 2) .* U(3))';
            sandPoroRange = cell2mat(FS.MatPropDistr.poro.range(~idc)');
            obj.MatProps.poroRange(:, ~idc) = sandPoroRange';
            obj.MatProps.poro(~idc) = (sandPoroRange(:, 1) + ...
                                      diff(sandPoroRange, [], 2) .* U(2))';
            
            % Clay perm depends on poro, so we correlate poro too.
            clayPoroRange = cell2mat(FS.MatPropDistr.poro.range(idc)');
            obj.MatProps.poroRange(:, idc) = clayPoroRange';
            obj.MatProps.poro(idc) = (clayPoroRange(:, 1) + ...
                                     diff(clayPoroRange, [], 2) .* U(2))';
            if isa(FS.MatPropDistr.perm.range{find(idc, 1)}, 'function_handle')
                clayPermRange = cell2mat(cellfun(@(x, id) x(clayPoroRange(id,1), ...
                                                            clayPoroRange(id,2)), ...
                                                 FS.MatPropDistr.perm.range(idc), ...
                                                 num2cell((1:sum(idc))), ...
                                                 'UniformOutput', false)');
            elseif all(~isa(FS.MatPropDistr.perm.range(idc), 'function_handle'))
                clayPermRange = cell2mat(FS.MatPropDistr.perm.range(idc)');
                warning(['Since clay permeability was passed, ' ...
                         'not correlated to poro (nor fT)'])
            else
                error(['For now, the number of provided clay perm vals ' ...
                       'is all or none. Add loop for only some provided.'])
            end
            obj.MatProps.permRange(:, idc) = clayPermRange';
            clayPermRangeLog = log10(clayPermRange);
            obj.MatProps.perm(idc) = 10.^(clayPermRangeLog(:, 1) + ...
                                   diff(clayPermRangeLog, [], 2) .* U(3))';
            
            
            % ------------------------------------------------------------
            % Perm Anisotropy Ratio
            % ------------------------------------------------------------
            % poro at zf
            poroDist = getPorosity(FS.Vcl, FS.IsClayVcl, zf, 'zf', ...
                                   zf, FS.IsUndercompacted, FS.HW.Id);
            poroAtZf = cellfun(@(x) x(1), poroDist.fcn);
            
            % Perm anisotropy ratio
            gamma = obj.Disp / obj.MatProps.thick;
            obj.MatProps.permAnisoRatio = ...
                            cellfun(@(x, id) x(gamma, poroAtZf(id)), ...
                                    FS.MatPropDistr.permAnisoRatio.fcn, ...
                                    num2cell((1:FS.HW.Id(end))));
            
            % ------------------------------------------------------------
            % ResFric and SSFc (correlated)
            % ------------------------------------------------------------
            n = sum(idc);
            
            % 1. Sands (independent since SSFc only for smearing sources)
            obj.MatProps.resFric(~idc) = cellfun(@(x) x(1), ...
                                        FS.MatPropDistr.resFric.fcn(~idc));
            obj.MatProps.ssfc(~idc) = nan;
                                    
            % 2. Clays (correlated)     
            % Generate bivariate samples using a Gaussian copula
            if opt.corrCoef > 0
                U = copularnd('Gaussian', [1 -opt.corrCoef; ...
                                           -opt.corrCoef 1], n);                
            else 
                U = copularnd('Gaussian', [1 opt.corrCoef; ...
                                           opt.corrCoef 1], n); 
            end

            % Assert that resFric marginal distros are all beta
            assert(all(cellfun(@(x) strcmp(x, 'beta'), ...
                               FS.MatPropDistr.resFric.type(idc))));
            
            % Transform bivariate samples to marginals and ranges
            % ResFric
            resFricParam = cell2mat(FS.MatPropDistr.resFric.param(idc)');
            X1 = betainv(U(:, 1), resFricParam(:, 1), resFricParam(:, 2));
            resFricRange = cell2mat(FS.MatPropDistr.resFric.range(idc)');
            obj.MatProps.resFric(idc) = (resFricRange(:, 1) + ...
                                        diff(resFricRange, [], 2) .* X1)'; 
            
            % SSFc (this should work directly, even if the distributions
            % for each layer were different!)
            obj.MatProps.ssfc(idc) = cellfun(@(x, id) icdf(x, U(id, 2)), ...
                                     FS.MatPropDistr.ssfc.dist(idc), ...
                                     num2cell((1:sum(idc)))); 

           % ------------------------------------------------------------
           % Angles
           % ------------------------------------------------------------
            [obj.Alpha, obj.Delta, obj.Zeta] = ...
                                    getFaultAngles(obj.Throw, obj.Dip, ...
                                                   obj.MatProps.thick);
        end
        
        function obj = placeMaterials(obj, FS, smear, G)
            %
            % Place materials in the fault zone, assign permeabilities to
            % each fault material. See documentation in used functions 
            % below for details.
            %
            obj.Grid.cellDim = G.cellDim;
            
            % Mapping matrix (material in each grid cell)
            % MatMap makes this object somewhat heavy on RAM.
            obj.MatMap = faultMaterialMap(G, FS, smear);
            
            % Smear placement (object simulation)
            if any(obj.MatMap.Psmear < 1)
                    tol = 0.025;
                    obj.MatMap = placeSmearObjects(obj.MatMap, smear, FS, ...
                                                   G, tol, 0);
                else
                    if isempty(obj.MatMap.Psmear)
                        disp('No smear: P(smear) = 0')
                    %else
                    %    disp('Continuous smears: P(smear) = 1')
                    end
                    obj.MatMap.vals = transpose(obj.MatMap.vals);
                    obj.MatMap.P = [obj.MatMap.Psmear; ...  % desired (calculated)
                                    obj.MatMap.Psmear];     % obtained
            end
            
            % Assign vcl, porosity and permeability          
            [obj.Grid.poro, obj.Grid.perm, obj.Grid.permy, ...
             obj.Grid.vcl, obj.Grid.units] = setGridPoroPerm(obj, G, FS);
        end
        
        function obj = upscaleProps(obj, G, U)
            %
            % Upscale Vcl, Poro and Perm
            %
            
            % Upscale Vcl (additive)
            obj.Vcl = mean(obj.Grid.vcl);
            
            % Upscale Porosity (additive)
            % --------------------------------------------------------
            % If upscaled grid had more than one cell we would do:
            %     mrstModule add coarsegrid
            %     p = partitionCartGrid(G.cartDims, [1 1]);
            %     CG = generateCoarseGrid(G, p);
            %     Poro = accumarray(p, poroG)./accumarray(p,1);
            % Since it only has one cell:
            % --------------------------------------------------------
            obj.Poro = mean(obj.Grid.poro);
            
            % Upscale Permeability (not an additive property)
            obj.Perm = computeCoarsePerm(G, obj.Grid.perm, ...
                                         obj.Grid.permy, U, ...
                                         obj.Disp, obj.MatProps.thick);
            
        end
        
        function plotMaterials(obj, FS, G0)
           %
           %
           %
           
           % Generate Grid (Must be same as grid used within Fault)
           G = updateGrid(G0, obj.MatProps.thick);
           
           % utils
           M = obj.MatMap;
           latx = {'Interpreter', 'latex'};
           sz = [14, 11];
           %ydim = max(G.faces.centroids(:,2));
           rock.poro = obj.Grid.poro;
           rock.perm = obj.Grid.perm;
           idGrid = reshape(transpose(flipud(M.vals)), G.cells.num, 1);
           
           % SUBPLOTS
           % 1. Parent Ids
           hh = figure(randi(1000, 1));
           tiledlayout(1, 6, 'Padding', 'tight', 'TileSpacing', 'tight');
           nexttile
           set(gca, 'colormap', hot(max(M.unit)));
           plotToolbar(G, reshape(transpose(flipud(M.units)), G.cells.num, 1), ...
                       'EdgeColor', [0.2 0.2 0.2], 'EdgeAlpha', 0.1);
           xlim([0 obj.MatProps.thick]); ylim([0 obj.Disp]); c = colorbar;
           set(c,'YTick', sort(unique(M.unit)));
           xticks([0 obj.MatProps.thick])
           xticklabels([0 round(obj.MatProps.thick, 2)])
           xlabel('$x$ [m]', latx{:}); ylabel('$z$ [m]', latx{:})
           set(gca,'fontSize', 10)
           title('Parent Id', latx{:}, 'fontSize', sz(2));
           
           
           % 2. Smear vs sand
           nexttile
           layerTops = M.layerTop;
           layerBots = M.layerBot;
           layerCtrs = M.layerCenter;
           cmap = copper(2);
           if unique(idGrid) == 0
               cmap = cmap(2, :);
               set(gca, 'colormap', cmap);
           elseif unique(idGrid) == 1
               cmap = cmap(1, :);
               set(gca, 'colormap', cmap);
           else
               set(gca, 'colormap', flipud(cmap));
           end 
           plotToolbar(G, idGrid, 'EdgeColor', [0.2 0.2 0.2], ...
                       'EdgeAlpha', 0);
           hold on
           for n=FS.FW.Id
               plot(0,layerTops(n), '+r', 'MarkerSize', 4, 'LineWidth', 1)
               plot(0,layerCtrs(n), 'or', 'MarkerSize', 3, ...
                    'LineWidth', 1, 'MarkerFaceColor', 'r')
           end
           plot(0, layerBots(1), '+r', 'MarkerSize', 4, 'LineWidth', 1)
           for n=FS.HW.Id
               plot(obj.MatProps.thick, layerTops(n), '+r', ...
                    'MarkerSize', 4, 'LineWidth', 1)
               plot(obj.MatProps.thick, layerCtrs(n), 'or', ...
                    'MarkerSize', 3, 'LineWidth', 1, 'MarkerFaceColor', 'r')
           end
           plot(obj.MatProps.thick, layerBots(1), '+r', 'MarkerSize', ...
                4, 'LineWidth', 1)
           hold off
           xlim([0 obj.MatProps.thick]); ylim([0 obj.Disp]); ...
           axis off
           c = colorbar;
           set(c,'YTick', [0 1]);
           if unique(idGrid) == 0
               caxis([0 0.1]);
           elseif unique(idGrid) == 1
               caxis([0.9 1]);
           end
           %xlabel('$x$ [m]', latx{:}); ylabel('$z$ [m]', latx{:})
           set(gca,'fontSize', 10)
           title('Material', latx{:}, 'fontSize', sz(2));
           
           % 3. Cell and upscaled Vcl
           nexttile
           set(gca, 'colormap', flipud(copper))
           plotToolbar(G, obj.Grid.vcl, 'EdgeColor', [0.2 0.2 0.2], ...
                       'EdgeAlpha', 0);
           xlim([0 obj.MatProps.thick]); ylim([0 obj.Disp]); 
           axis off
           c = colorbar;
           caxis([0 1]);
           set(gca,'fontSize', 10)
           %c.Label.Interpreter = 'latex'; 
           %c.Label.String = '$n$ [-]';
           %c.Label.FontSize = 12;
           %xlabel('$x$ [m]', latx{:}); ylabel('$z$ [m]', latx{:})
           title(['f$_{V_\mathrm{cl}} =$ ' num2str(obj.Vcl, ' %1.2f') ' [-]'], latx{:}, ...
                 'fontSize', sz(2));
           
           
           % 4. Cell and upscaled porosity
           nexttile
           set(gca, 'colormap', copper)
           plotToolbar(G, rock.poro, 'EdgeColor', [0.2 0.2 0.2], ...
                       'EdgeAlpha', 0);
           xlim([0 obj.MatProps.thick]); ylim([0 obj.Disp]); 
           axis off
           c = colorbar;
           caxis([0 max([max(rock.poro), 0.25])]);
           set(gca,'fontSize', 10)
           %c.Label.Interpreter = 'latex'; 
           %c.Label.String = '$n$ [-]';
           %c.Label.FontSize = 12;
           %xlabel('$x$ [m]', latx{:}); ylabel('$z$ [m]', latx{:})
           title(['$n =$ ' num2str(obj.Poro, ' %1.2f') ' [-]'], latx{:}, ...
                 'fontSize', sz(2));
           
           % 5. Cell and upscaled permeability
           nexttile
           set(gca, 'colormap', copper)
           plotToolbar(G, log10(rock.perm(:,1)/(milli*darcy)), ...
                       'EdgeColor', [0.2 0.2 0.2], 'EdgeAlpha', 0);
           xlim([0 obj.MatProps.thick]); ylim([0 obj.Disp]); %axis off
           c = colorbar;
           %if ~any(obj.MatMap.isclay)
                caxis([min(log10(rock.perm(:,1)/(milli*darcy))) ...
                       max(log10(rock.perm(:,3)/(milli*darcy)))]);
           %else
           %    caxis([min(log10(rock.perm(:,1)/(milli*darcy))) 2]);
           %end
           c.Label.Interpreter = 'latex'; 
           c.Label.String = '$\log_{10} k_{xx}$ [mD]';
           c.Label.FontSize = 12;
           set(gca,'fontSize', 10)
           xlabel('$x$ [m]', latx{:}); ylabel('$z$ [m]', latx{:})
           xticks([0 obj.MatProps.thick])
           xticklabels([0 round(obj.MatProps.thick, 2)])
           val = obj.Perm(1)/(milli*darcy);
           if val < 1e-3
               val = val*1000;
               units = ' [$\mu$D]';
               form = ' %1.3f';
           elseif val > 999
               val = val/1000;
               units = ' [D]';
               form = ' %1.2f';
           else
               units = ' [mD]';
               form = ' %3.3f';
           end
           title(['$k_{xx} =$ ' num2str(val, form) units], latx{:}, ...
                 'fontSize', sz(2));
           
           % 6. kzz
           nexttile
           set(gca, 'colormap', copper)
           plotToolbar(G, log10(rock.perm(:,3)/(milli*darcy)), ...
                       'EdgeColor', [0.2 0.2 0.2], 'EdgeAlpha', 0);
           xlim([0 obj.MatProps.thick]); ylim([0 obj.Disp]); 
           axis off
           c = colorbar;
           %if ~any(obj.MatMap.isclay) % if there is sand
               caxis([min(log10(rock.perm(:,1)/(milli*darcy))) ...
                      max(log10(rock.perm(:,3)/(milli*darcy)))]);
           %else
           %   caxis([min(log10(rock.perm(:,1)/(milli*darcy))) 2]);
           %end
           c.Label.Interpreter = 'latex'; 
           c.Label.String = '$\log_{10} k_{zz}$ [mD]';
           c.Label.FontSize = 12;
           set(gca,'fontSize', 10)
           %xlabel('$x$ [m]', latx{:}, 'fontSize', 14); 
           %ylabel('$z$ [m]', latx{:}, 'fontSize', 14)
           val = obj.Perm(3)/(milli*darcy);
           if val < 1e-3
               val = val*1000;
               units = ' [$\mu$D]';
               form = ' %1.3f';
           elseif val > 999
               val = val/1000;
               units = ' [D]';
               form = ' %1.2f';
           else
               units = ' [mD]';
               form = ' %3.3f';
           end
           title(['$k_{zz} =$ ' num2str(val, form) units], latx{:}, ...
                 'fontSize', sz(2));
           set(hh, 'position', [200, 0, 900, 350]);       
           
           % MatProps 
           %Thick = repelem(obj.MatProps.thick, numel(obj.MatProps.resFric));
           %bounds = cell2mat(FS.MatPropDistr.ssfc.range(idc)');
           %SSFcMin = bounds(:, 1);
           %SSFcMax = bounds(:, 2);
           PoroMin = obj.MatProps.poroRange(1, :);
           PoroMax = obj.MatProps.poroRange(2, :);
           N = max(FS.HW.Id);
           PermMin = obj.MatProps.permRange(1, :);
           PermMax = obj.MatProps.permRange(2, :);
           tdata = table(obj.MatProps.resFric', obj.MatProps.ssfc', ...
                         PoroMin', PoroMax', PermMin', PermMax', ...
                         obj.MatProps.permAnisoRatio');
           tdata.Properties.VariableNames = {'ResFric', 'SSFc', ...
                                             'PoroMin', 'PoroMax', ...
                                             'PermMin', 'PermMax', ...
                                             'PermAniso'}; 
           fig = uifigure(randi(1000, 1), 'Position', [500 500 700 30*N]);
           uit = uitable(fig);
           uit.Position = [20 20 550 30*N-40];
           uit.Data = tdata;
           uit.RowName = 'numbered';
                     
           
        end
        
    end
end