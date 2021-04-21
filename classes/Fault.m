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
       Poro                             % Final upscaled Poro (Nsim, 1)
       Perm                             % Final upscaled Perm (Nsim, 3)
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
        function obj = Fault(FS, dip, targetCellDim)
            %
            % We instantiate a Fault object with the fundamental 
            % dimensions.
            %
            
            % Geometry
            obj.Dip    = dip;
            obj.Disp   = sum(FS.Tap(FS.FW.Id));   
            
            % User can pass different grid resolution
            if nargin < 3
                obj.Grid.TargetCellDim = [0.1, 1];     % [m], [thick., disp.]
            else
                assert(isa(targetCellDim, 'double') && numel(targetCellDim) == 2, ...
                       ['If grid resolution is passed, it must be a ', ...
                       'double array with 2 elements: ', ...
                       '[res. in thick. dir. (x), res. in displ. dir (y)]'])
                obj.Grid.TargetCellDim = targetCellDim;
            end
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
            opt.maxPerm = [];                   % mD
            opt.siltInClay = false;
            opt.isUndercompacted = false;
            opt = merge_options_relaxed(opt, varargin{:});
            
            % Shorten input names
            zf   = FS.DepthFaulting;
            zmax = [FS.FW.DepthBurial, FS.HW.DepthBurial];
            clayMine = {FS.FW.ClayMine, FS.HW.ClayMine};
            
            % Fault thickness
            obj.MatProps.Thick = getFaultThickness(obj.Disp, 1);
            
            % Angles
            [obj.Alpha, obj.Delta, obj.Zeta] = ...
                                    getFaultAngles(obj.Throw, obj.Dip, ...
                                                   obj.MatProps.Thick);
            
            % ResFric
            obj.MatProps.ResFric = getResidualFrictionAngle(FS.Vcl, FS);
            
            % SSFc and SSFc bounds
            [obj.MatProps.SSFc, obj.MatProps.SSFcBounds] = ...
                                getSSFc(FS.Vcl, FS.IsClayVcl, zf, ...
                                        FS.Thick, obj.Disp, FS.HW.Id);
                                              
            % Perm Anisotropy Ratio
            poroAtZf = getPorosity(FS.Vcl, FS.IsClayVcl, zf, zmax, ...
                                   'zf', opt.isUndercompacted, FS.HW.Id);
            gamma = obj.Disp / obj.MatProps.Thick;
            obj.MatProps.PermAnisoRatio = ...
                getAnisotropyRatio(FS.Vcl, FS.IsClayVcl, zf, clayMine, ...
                                   gamma, opt.siltInClay, poroAtZf, FS.HW.Id);
            % Porosity
            obj.MatProps.Poro = getPorosity(FS.Vcl, FS.IsClayVcl, ...         
                                            zf, zmax, 'zmax', ...
                                            opt.isUndercompacted, FS.HW.Id);
            % Perm
            obj.MatProps.Perm = getPermeability(FS.Vcl, FS.IsClayVcl, ...
                                                zf, zmax, opt.maxPerm, ...
                                                obj.MatProps.Poro, FS);
        end
        
        function obj = upscaleSmearPerm(obj, FS, smear, U)
            %
            %   
            %
            
            % Generate Grid
            G = makeFaultGrid(obj.MatProps.Thick, obj.Disp, ...
                              obj.Grid.TargetCellDim);
            obj.Grid.CellDim = G.CellDim;
            
            % Mapping matrix (material in each grid cell)
            % MatMap makes this somewhat heavy on RAM.
            obj.MatMap = faultMaterialMap(G, FS, smear.DomainLength, ...
                                          smear.ThickInFault, smear.Psmear);
            
            % Smear placement (object simulation)
            if any(obj.MatMap.Psmear < 1)
                tol = 0.001;
                obj.MatMap = placeSmearObjects(obj.MatMap, G, tol, ...
                                                smear.Length, ...
                                                smear.SegLenMax, ...
                                                smear.DomainLength, 0);
            else
                if isempty(obj.MatMap.Psmear)
                    disp('No smear: P(smear) = 0')
                else
                    disp('Continuous smears: P(smear) = 1')
                end
                obj.MatMap.vals = transpose(obj.MatMap.vals);
                obj.MatMap.P = [obj.MatMap.Psmear; ...  % desired (calculated)
                                obj.MatMap.Psmear];     % obtained
            end
            
            % Assign porosity and permeability
            [obj.Grid.Poro, obj.Grid.Perm, kz_loc] = ...
                                               setGridPoroPerm(obj, G, FS);
            
            % Upscale Porosity (additive)
            % --------------------------------------------------------
            % If upscaled grid had more than one cell we would do:
            %     mrstModule add coarsegrid
            %     p = partitionCartGrid(G.cartDims, [1 1]);
            %     CG = generateCoarseGrid(G, p);
            %     Poro = accumarray(p, poroG)./accumarray(p,1);
            % Since it only has one cell:
            % --------------------------------------------------------
            obj.Poro = mean(obj.Grid.Poro);
            
            % Upscale Permeability (not an additive property)
            obj.Perm = computeCoarsePerm(G, obj.Grid.Perm, kz_loc, U, ...
                                         obj.Disp, obj.MatProps.Thick);
        end
        
        function plotMaterials(obj, FS)
           %
           %
           %
           
           % Generate Grid (Must be same as grid used within Fault)
            G = makeFaultGrid(obj.MatProps.Thick, obj.Disp, ...
                              obj.Grid.TargetCellDim);
           
           % utils
           M = obj.MatMap;
           latx = {'Interpreter', 'latex'};
           sz = [14, 11];
           %ydim = max(G.faces.centroids(:,2));
           rock.poro = obj.Grid.Poro;
           rock.perm = obj.Grid.Perm;
           idGrid = reshape(transpose(flipud(M.vals)), G.cells.num, 1);
           
           % SUBPLOTS
           % 1. Parent Ids
           hh = figure(randi(1000, 1));
           tiledlayout(1, 5, 'Padding', 'compact', 'TileSpacing', 'none');
           nexttile
           set(gca, 'colormap', hot(max(M.unit)));
           plotToolbar(G, reshape(transpose(flipud(M.units)), G.cells.num, 1), ...
                       'EdgeColor', [0.2 0.2 0.2], 'EdgeAlpha', 0.1);
           xlim([0 obj.MatProps.Thick]); ylim([0 obj.Disp]); c = colorbar;
           set(c,'YTick', sort(unique(M.unit)));
           xticks([0 obj.MatProps.Thick])
           xticklabels([0 round(obj.MatProps.Thick, 2)])
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
               plot(obj.MatProps.Thick, layerTops(n), '+r', ...
                    'MarkerSize', 4, 'LineWidth', 1)
               plot(obj.MatProps.Thick, layerCtrs(n), 'or', ...
                    'MarkerSize', 3, 'LineWidth', 1, 'MarkerFaceColor', 'r')
           end
           plot(obj.MatProps.Thick, layerBots(1), '+r', 'MarkerSize', ...
                4, 'LineWidth', 1)
           hold off
           xlim([0 obj.MatProps.Thick]); ylim([0 obj.Disp]); ...
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
           
           
           % 3. Cell and upscaled porosity
           nexttile
           set(gca, 'colormap', copper)
           plotToolbar(G, rock.poro, 'EdgeColor', [0.2 0.2 0.2], ...
                       'EdgeAlpha', 0);
           xlim([0 obj.MatProps.Thick]); ylim([0 obj.Disp]); 
           axis off
           c = colorbar;
           caxis([0 max([max(rock.poro), 0.35])]);
           set(gca,'fontSize', 10)
           %c.Label.Interpreter = 'latex'; 
           %c.Label.String = '$n$ [-]';
           %c.Label.FontSize = 12;
           %xlabel('$x$ [m]', latx{:}); ylabel('$z$ [m]', latx{:})
           title(['$n =$ ' num2str(obj.Poro, ' %1.2f') ' [-]'], latx{:}, ...
                 'fontSize', sz(2));
           
           % 4. Cell and upscaled permeability
           nexttile
           set(gca, 'colormap', copper)
           plotToolbar(G, log10(rock.perm(:,1)/(milli*darcy)), ...
                       'EdgeColor', [0.2 0.2 0.2], 'EdgeAlpha', 0);
           xlim([0 obj.MatProps.Thick]); ylim([0 obj.Disp]); axis off
           c = colorbar;
           if ~any(obj.MatMap.isclay)
                caxis([min(log10(rock.perm(:,1)/(milli*darcy))) ...
                       max(log10(rock.perm(:,1)/(milli*darcy)))]);
           else
               caxis([min(log10(rock.perm(:,1)/(milli*darcy))) 1]);
           end
           c.Label.Interpreter = 'latex'; 
           c.Label.String = '$\log_{10} k_{xx}$ [mD]';
           c.Label.FontSize = 12;
           set(gca,'fontSize', 10)
           %xlabel('$x$ [m]', latx{:}); ylabel('$z$ [m]', latx{:})
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
           
           nexttile
           set(gca, 'colormap', copper)
           plotToolbar(G, log10(rock.perm(:,3)/(milli*darcy)), ...
                       'EdgeColor', [0.2 0.2 0.2], 'EdgeAlpha', 0);
           xlim([0 obj.MatProps.Thick]); ylim([0 obj.Disp]); 
           axis off
           c = colorbar;
           if ~any(obj.MatMap.isclay) % if there is sand
               caxis([min(log10(rock.perm(:,1)/(milli*darcy))) ...
                      max(log10(rock.perm(:,3)/(milli*darcy)))]);
           else
              caxis([min(log10(rock.perm(:,1)/(milli*darcy))) 1]);
           end
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
           set(hh, 'position', [200, 0, 650, 350]);       
           
           % MatProps 
           %Thick = repelem(obj.MatProps.Thick, numel(obj.MatProps.ResFric));
           SSFcMin = obj.MatProps.SSFcBounds(1, :);
           SSFcMax = obj.MatProps.SSFcBounds(2, :);
           PoroMin = obj.MatProps.Poro(:, 1);
           PoroMax = obj.MatProps.Poro(:, 2);
           N = max(FS.HW.Id);
           PermMin = zeros(N, 1);
           PermMax = zeros(N, 1);
           for n = 1:N
              samples = obj.MatProps.Perm{n}(5000);
              PermMin(n) = min(samples)/(milli*darcy);
              PermMax(n) = max(samples)/(milli*darcy);
           end
           tdata = table(obj.MatProps.ResFric', obj.MatProps.SSFc', ...
                         PoroMin, PoroMax, PermMin, PermMax, ...
                         obj.MatProps.PermAnisoRatio');
           tdata.Properties.VariableNames = {'ResFric', 'SSFc', ...
                                             'PoroMin', 'PoroMax', ...
                                             'PermMin', 'PermMax', ...
                                             'PermAniso'}; 
           fig = uifigure(randi(1000, 1), 'Position', [500 500 620 30*N]);
           uit = uitable(fig);
           uit.Position = [20 20 550 30*N-40];
           uit.Data = tdata;
           uit.RowName = 'numbered';
                     
           
        end
        
    end
end