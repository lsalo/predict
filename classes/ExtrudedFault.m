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
        
        function [obj, CG] = upscaleProps(obj, G, U)
            %
            %
            %
            
            % Upscale Vcl and Porosity (additive)
            p = partitionCartGrid(G.cartDims, U.coarseDims);
            if all(U.coarseDims == 1)
                obj.Vcl = mean(obj.Grid.vcl);
                obj.Poro = mean(obj.Grid.poro);
            else    % upscaled grid has more than one cell
                CG = generateCoarseGrid(G, p);
                % The following is correct for grids with uniform cell size
                obj.Vcl = accumarray(p, obj.Grid.vcl)./accumarray(p,1);
                obj.Poro = accumarray(p, obj.Grid.poro)./accumarray(p,1);
            end
            
            % Upscale Permeability (not an additive property)
            obj.Perm = computeCoarsePerm3D(G, obj.Grid.perm, U, p);
               
        end
        
        function plotMaterials(obj, faultSection, FS)
           %
           %
           %
           
           % utils
           M = faultSection.MatMap;
           latx = {'Interpreter', 'latex'};
           sz = [14, 11];
           
           % Create grid
           G = makeFaultGrid(obj.Thick, obj.Disp, obj.Length, obj.SegLen);
           
           % SUBPLOTS
           % 1. Parent Ids
           hh = figure(randi(1000, 1));
           tiledlayout(1, 3, 'Padding', 'tight', 'TileSpacing', 'tight');
           nexttile
           set(gca, 'colormap', hot(max(obj.Grid.units)));
           plotToolbar(G, obj.Grid.units, 'EdgeColor', [0.2 0.2 0.2], ...
                       'EdgeAlpha', 0.1);
           xlim([0 obj.Thick]); ylim([0 obj.Disp]); c = colorbar;
           set(c,'YTick', unique(obj.Grid.units));
           xticks([0 obj.Thick])
           xticklabels([0 round(obj.Thick, 1)])
           xlabel('$x$ [m]', latx{:}); 
           ylabel('$y$ [m]', latx{:})
           zlabel('$z$ [m]', latx{:})
           ax = gca;
           ax.DataAspectRatio = [0.1 1 1];
           ax.ZDir = 'normal';
           set(gca,'fontSize', 10)
           title('Parent Id', latx{:}, 'fontSize', sz(2));
           view([30 20])
           
           
           % 2. Smear vs sand
           nexttile
           layerTops = M.layerTop;
           layerBots = M.layerBot;
           layerCtrs = M.layerCenter;
           cmap = copper(2);
           if unique(obj.Grid.isSmear) == 0
               cmap = cmap(2, :);
               set(gca, 'colormap', cmap);
           elseif unique(obj.Grid.isSmear) == 1
               cmap = cmap(1, :);
               set(gca, 'colormap', cmap);
           else
               set(gca, 'colormap', flipud(cmap));
           end 
           plotToolbar(G, obj.Grid.isSmear, 'EdgeColor', [0.2 0.2 0.2], ...
                       'EdgeAlpha', 0.1);
           hold on
           for n=FS.FW.Id
               plot3([0 0],[0 obj.Length],[layerTops(n) layerTops(n)], ...
                     '-w', 'LineWidth', 1)
               plot3([0 0],[0 obj.Length],[layerCtrs(n) layerCtrs(n)], ...
                     '-', 'color', [0.7 0.7 0.7], 'LineWidth', 1)
           end
           plot3([0 0], [0 obj.Length], [layerBots(1) layerBots(1)], ...
                 '-w', 'LineWidth', 1)
           for n=FS.HW.Id
               plot3([obj.Thick obj.Thick], [0 obj.Length], ...
                    [layerTops(n) layerTops(n)], '-w', 'LineWidth', 1)
               plot3([obj.Thick obj.Thick], [0 obj.Length], ...
                     [layerCtrs(n) layerCtrs(n)], '-', ...
                     'color', [0.7 0.7 0.7], 'LineWidth', 1)
           end
           plot3([obj.Thick obj.Thick], [0 obj.Length], ...
                [layerBots(1) layerBots(1)], '-w', 'MarkerSize', ...
                4, 'LineWidth', 1)
           hold off
           xlim([0 obj.Thick]); ylim([0 obj.Length]); zlim([0 obj.Disp])
           axis off
           c = colorbar;
           set(c,'YTick', [0 1]);
           if unique(obj.Grid.isSmear) == 0
               caxis([0 0.1]);
           elseif unique(obj.Grid.isSmear) == 1
               caxis([0.9 1]);
           end
           ax = gca;
           ax.DataAspectRatio = [0.1 1 1];
           ax.ZDir = 'normal';
           set(gca,'fontSize', 10)
           view([30 20])
           title('Material (all)', latx{:}, 'fontSize', sz(2));
           
           % 3. Only smears, coloured by vcl
           nexttile
           cmap = copper;
           set(gca, 'colormap', flipud(cmap(1:128, :)));
           plotToolbar(G, obj.Grid.vcl, obj.Grid.isSmear, ...
                       'EdgeColor', [0.8 0.8 0.8], 'EdgeAlpha', 0.1);
           xlim([0 obj.Thick]); ylim([0 obj.Length]); zlim([0 obj.Disp])
           axis off
           c = colorbar;
           ax = gca;
           ax.DataAspectRatio = [0.1 1 1];
           ax.ZDir = 'normal';
           set(gca,'fontSize', 10)
           view([30 20])
           title('Vcl, smears only', latx{:}, 'fontSize', sz(2));
           set(hh, 'position', [200, 0, 900, 350]);
           
           % 4. Cell and upscaled Vcl
           hf = figure(randi(1000, 1));
           tiledlayout(1, 4, 'Padding', 'tight', 'TileSpacing', 'tight');
           nexttile
           set(gca, 'colormap', flipud(copper))
           plotToolbar(G, obj.Grid.vcl, 'EdgeColor', [0.2 0.2 0.2], ...
                       'EdgeAlpha', 0.1);
           xlim([0 obj.Thick]); ylim([0 obj.Disp]); 
           c = colorbar;
           caxis([min(obj.Grid.vcl) max(obj.Grid.vcl)]);
           set(gca,'fontSize', 10)
           %c.Label.Interpreter = 'latex'; 
           %c.Label.String = '$n$ [-]';
           %c.Label.FontSize = 12;
           xlabel('$x$ [m]', latx{:}); ylabel('$y$ [m]', latx{:});
           zlabel('$z$ [m]', latx{:})
           ax = gca;
           ax.DataAspectRatio = [0.1 1 1];
           ax.ZDir = 'normal';
           view([30 20])
           xticks([0 obj.Thick])
           xticklabels([0 round(obj.Thick, 1)])
           title(['f$_{V_\mathrm{cl}} =$ ' num2str(obj.Vcl, ' %1.2f') ' [-]'], latx{:}, ...
                  'fontSize', sz(2));
           
           % 5. Cell and upscaled porosity
           nexttile
           set(gca, 'colormap', copper)
           plotToolbar(G, obj.Grid.poro, 'EdgeColor', [0.2 0.2 0.2], ...
                       'EdgeAlpha', 0.1);
           xlim([0 obj.Thick]); ylim([0 obj.Disp]); 
           c = colorbar;
           caxis([0 max([max(obj.Grid.poro), 0.25])]);
           set(gca,'fontSize', 10)
           %c.Label.Interpreter = 'latex'; 
           %c.Label.String = '$n$ [-]';
           %c.Label.FontSize = 12;
           %xlabel('$x$ [m]', latx{:}); ylabel('$z$ [m]', latx{:})
           ax = gca;
           ax.DataAspectRatio = [0.1 1 1];
           ax.ZDir = 'normal';
           view([30 20])
           title(['$n =$ ' num2str(obj.Poro, ' %1.2f') ' [-]'], latx{:}, ...
                 'fontSize', sz(2));
           
           % 6. Cell and upscaled permeability
           nexttile
           cmap = copper;
           set(gca, 'colormap', cmap(1:128, :))
           plotToolbar(G, log10(obj.Grid.perm(:,1)/(milli*darcy)), ...
                       obj.Grid.isSmear, 'EdgeColor', [0.8 0.8 0.8], ...
                       'EdgeAlpha', 0.1);
           xlim([0 obj.Thick]); zlim([0 obj.Disp]); ylim([0 obj.Length]);
           c = colorbar;
           %caxis([min(log10(obj.Grid.perm(:,1)/(milli*darcy))) ...
           %       max(log10(obj.Grid.perm(:,4)/(milli*darcy)))]);
           c.Label.Interpreter = 'latex'; 
           c.Label.String = '$\log_{10} k_{xx}$ [mD] (Clay smears)';
           c.Label.FontSize = 12;
           set(gca,'fontSize', 10)  
           val = zeros(1,3);
           units = cell(1,3);
           form = cell(1,3);
           for n=1:3
               val(n) = obj.Perm(n)/(milli*darcy);
               if val(n) < 1e-3
                   val(n) = val(n)*1000;
                   units{n} = ' [$\mu$D]';
                   form{n} = ' %1.3f';
               elseif val(n) > 999
                   val(n) = val(n)/1000;
                   units{n} = ' [D]';
                   form{n} = ' %1.2f';
               else
                   units{n} = ' [mD]';
                   form{n} = ' %3.3f';
               end
           end
           ax = gca;
           ax.DataAspectRatio = [0.1 1 1];
           ax.ZDir = 'normal';
           view([30 20])
           title(['$k_{jj} =$ ' num2str(val(1), form{1}) units{1}, ...
                  ' $\vert$ ' num2str(val(2), form{2}) units{2}, ...
                  ' $\vert$ ' num2str(val(3), form{3}) units{3}], latx{:}, ...
                  'fontSize', sz(2));
           grid on
           
           % 7 Cell permeability (sands)
           nexttile
           cmap = copper;
           set(gca, 'colormap', cmap(156:end, :))
           plotToolbar(G, log10(obj.Grid.perm(:,1)/(milli*darcy)), ...
                       ~obj.Grid.isSmear, 'EdgeColor', [0.2 0.2 0.2], ...
                       'EdgeAlpha', 0.1);
           xlim([0 obj.Thick]); zlim([0 obj.Disp]); ylim([0 obj.Length]);
           c = colorbar;
           %caxis([min(log10(obj.Grid.perm(:,1)/(milli*darcy))) ...
           %       max(log10(obj.Grid.perm(:,4)/(milli*darcy)))]);
           c.Label.Interpreter = 'latex'; 
           c.Label.String = '$\log_{10} k_{xx}$ [mD] (Sand smears)';
           c.Label.FontSize = 12;
           set(gca,'fontSize', 10)
           ax = gca;
           ax.DataAspectRatio = [0.1 1 1];
           ax.ZDir = 'normal';
           view([30 20])
           grid on
           set(hf, 'position', [200, 0, 1200, 350]);
            
           
           % MatProps 
           %Thick = repelem(obj.MatProps.thick, numel(obj.MatProps.resFric));
           %bounds = cell2mat(FS.MatPropDistr.ssfc.range(idc)');
           %SSFcMin = bounds(:, 1);
           %SSFcMax = bounds(:, 2);
           PoroMin = faultSection.MatProps.poroRange(1, :);
           PoroMax = faultSection.MatProps.poroRange(2, :);
           N = max(FS.HW.Id);
           PermMin = faultSection.MatProps.permRange(1, :);
           PermMax = faultSection.MatProps.permRange(2, :);
           tdata = table(faultSection.MatProps.resFric', faultSection.MatProps.ssfc', ...
                         PoroMin', PoroMax', PermMin', PermMax', ...
                         faultSection.MatProps.permAnisoRatio');
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