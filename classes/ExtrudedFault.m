classdef ExtrudedFault
    %
    % SUMMARY:
    % 
    %   TBD
    %
    % SYNOPSIS:
    %   myFault = ExtrudedFault(TBD)
    %
    %
    % DESCRIPTION:
    % 
    %   TBD
    %
    % 
    % REQUIRED PARAMETERS:
    %   TBD
    %
    %   footwall:       Stratigraphy object with corresponding property 
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
        Thick                           % Dip-perpendicular fault thickness (m)
        SegLen                          % Strike-parallel length of material segments (m) 
    end
    
    properties (SetAccess = protected)
        Dip                              % Fault dip angle (deg)
        Disp                             % Fault displacement [m]
        Length                           % Fault length [m]
        Poro                             % Final upscaled Poro (nc, 1)
        Perm                             % Final upscaled Perm (nc, 3)
        Vcl                              % Final upscaled Vcl (nc, 1)
        Grid                             % contains cell quantities
    end
   
    
    methods
        function obj = ExtrudedFault(faultSection, FS)
            %
            % We instantiate an ExtrudedFault object with the fundamental 
            % dimensions and properties from any 2D faultSection
            %
            
            % Geometry
            obj.Dip    = faultSection.Dip;
            obj.Disp   = faultSection.Disp;  
            obj.Length = obj.Disp;
            obj.Thick = FS.MatPropDistr.thick.fcn(obj.Disp);
        end
        
        function obj = getSegmentationLength(obj, FS, nSeg)
            %
            % Get length of material segments in the strike direction
            %
            if nargin > 2   % prescribe number of segments of equal length
                obj.SegLen = repelem(round(obj.Length/nSeg), 1, nSeg);
            else            % Determine from material properties
                % TBD
            end
            
        end
        
        function obj = assignExtrudedVals(obj, G, faultSection, k)
            %
            % Assign cell-based values from 2D fault sections to 
            % 3D fault
            %
            
            [nx, ny, nz] = deal(G.cartDims(1), G.cartDims(2), G.cartDims(3));
            % Initial checks
            if ~isnan(G.cellDim(2))
                assert(mod(segLen/cellDim, 1)==0)
                nklayers = segLen/cellDim;                  % repeat 2D section for nklayers
            else
                assert(G.cartDims(2) == numel(obj.SegLen))  % 1 cell per along-strike segment  
                nklayers = 1;
            end
            if k == 1
                obj.Grid.isSmear = false(G.cells.num, 1);
                obj.Grid.units = zeros(G.cells.num, 1);
                obj.Grid.vcl = zeros(G.cells.num, 1);
                obj.Grid.poro = zeros(G.cells.num, 1);
                obj.Grid.perm = zeros(G.cells.num, 6);      %[kxx, kxy, kxz, kyy, kyz, kzz]
                idFirst0 = 1;
            else
                if ~isnan(G.cellDim(2))     % constant along-strike cell L, layeredGrid
                    assert(isfield(G, 'layerSize'))
                    idFirst0 = G.layerSize*sum(obj.SegLen(1:k-1)/G.cellDim(2)) + 1;
                else
                    assert(~isfield(G, 'layerSize'))            % cartGrid
                    idFirst0 = nx*(k-1) + 1;
                end
            end
            
            % Find cellIds for current layer(s)
            % cartGrid
            idCellsLayer1 = repmat((1:nx)', nz, 1) + idFirst0-1;
            idCellsLayer1 = idCellsLayer1 + repelem(nx*ny*(0:nz-1)', nx, 1);
            idCellsLayers = repmat(idCellsLayer1, 1, nklayers);
            if nklayers > 1
                idCellsLayers(:,2:end) = idCellsLayers(:,2:end) + ...
                    repelem(repelem(nx, numel(idCellsLayer1), 1), 1, nklayers-1) + ...
                    ((nx:nx:(nklayers-1)*nx) - nx);
                idCellsLayers = reshape(idCellsLayers, nx*nz*nklayers, 1);
            end
            
            % isSmear
            M = faultSection.MatMap;
            nCellsLayer = G.cartDims(1)*G.cartDims(3);
            obj.Grid.isSmear(idCellsLayers) = repmat(reshape(...
                  transpose(flipud(M.vals)), nCellsLayer, 1), nklayers, 1);
            
            % Parent Unit
            unitsv = faultSection.Grid.units;   % incl. unitInClayGaps
            obj.Grid.units(idCellsLayers) = repmat(unitsv, nklayers, 1);
            
            % Vcl
            vclv = faultSection.Grid.vcl;
            obj.Grid.vcl(idCellsLayers) = repmat(vclv, nklayers, 1);
            
            % Poro
            porov = faultSection.Grid.poro;
            obj.Grid.poro(idCellsLayers) = repmat(porov, nklayers, 1);
            
            % Perm
            % cartGrid
            kxx = faultSection.Grid.perm(:,1);
            kxz = faultSection.Grid.perm(:,2);
            kyy = faultSection.Grid.permy;
            kzz = faultSection.Grid.perm(:,3);
            idPerm3D = [1 3 4 6];                       % rotated about the y axis (along-strike)
            % layeredGrid (extruded)
            % kyy = faultSection.Grid.perm(:,1);
            % kyz = faultSection.Grid.perm(:,2);
            % kzz = faultSection.Grid.perm(:,3);
            % kxx = faultSection.Grid.permy;
            % idPerm3D = [1 4 5 6];                       % rotated about the x axis (along-strike)
            obj.Grid.perm(idCellsLayers, idPerm3D) = ...
                                 repmat([kxx, kxz, kyy, kzz], nklayers, 1);
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