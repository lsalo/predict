classdef FaultedSection
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
    
    
    properties (SetAccess = protected)
        HW              % Hangingwall object from Stratigraphy class
        FW              % Footwall object from Stratigraphy class
        Tap             % Apparent layer thickness on the fault
        Thick           % True layer thickness
        MatPropDistr    % Distros of material properties to sample from
        IsUndercompacted% Undercompacted section? Default: false
        MaxPerm         % Perm Cap for FZ materials? Default: [] (none)
        Failure         % shear (default) or hybrid (near-surface)
    end
    
    properties (Dependent, Hidden)
       ParentId 
       Vcl
       IsClayVcl
       DepthFaulting
       DepthBurial
    end
    
    methods
        function obj = FaultedSection(footwall, hangingwall, faultDip, varargin)
            % We instantiate the object with the FW and HW objects required
            % to construct it.
            %
            % Example: mySect = FaultedSection(footwall, hangingwall)
            %
            
            % Optional inputs
            opt.maxPerm = [];                   % mD
            opt.isUndercompacted = false;
            opt.failure = 'shear';
            opt = merge_options_relaxed(opt, varargin{:});
            obj.IsUndercompacted = opt.isUndercompacted;
            obj.MaxPerm = opt.maxPerm;
            obj.Failure = opt.failure;
            
            % Assign or compute apparent and true thicknesses
            if footwall.IsThickApp == 1
                TapFW = footwall.Thickness;
                TFW = getThick(footwall, hangingwall, faultDip);
            else
                TapFW = getApparentThick(footwall, hangingwall, faultDip);
                TFW = footwall.Thickness;
            end
            if hangingwall.IsThickApp == 1
                TapHW = hangingwall.Thickness;
                [~, THW] = getThick(footwall, hangingwall, faultDip);
            else
                [~, TapHW] = getApparentThick(footwall, hangingwall, ...
                                              faultDip);
                THW = hangingwall.Thickness;
            end
            
            % Checks
            assert(hangingwall.IsHW == 1)
            assert(abs(sum(TapFW) - sum(TapHW)) < 1e-4)
            assert(footwall.IsClayVcl == hangingwall.IsClayVcl);  % assumed by MatProps
            
            % Assign
            obj.FW = footwall;
            obj.HW = hangingwall;
            obj.Tap = [TapFW, TapHW];
            obj.Thick = [TFW, THW];
        end
        
        function parentId = get.ParentId(obj)
            % get parent Ids, ie the corresponding protolith Ids for our 
            % fault materials.
           parentId = [obj.FW.Id obj.HW.Id]; 
        end
        
        function vcl = get.Vcl(obj)
            % 
           vcl = [obj.FW.Vcl obj.HW.Vcl]; 
        end
        
        function isClayVcl = get.IsClayVcl(obj)
            % 
           assert(obj.FW.IsClayVcl == obj.HW.IsClayVcl)
           isClayVcl = obj.FW.IsClayVcl; 
        end
        
        function zf = get.DepthFaulting(obj)
            % 
           zf = [obj.FW.DepthFaulting obj.HW.DepthFaulting]; 
        end
        
        function zf = get.DepthBurial(obj)
            % 
           zf = [obj.FW.DepthBurial obj.HW.DepthBurial]; 
        end
        
        function obj = getMatPropDistr(obj)
           %
           % TBD
           %
            
           % Shorten input names
           zf   = obj.DepthFaulting;
           zmax = [obj.FW.DepthBurial, obj.HW.DepthBurial];
           clayMine = {obj.FW.ClayMine, obj.HW.ClayMine};
           Disp = sum(obj.Tap(obj.FW.Id));
           
           % Fault length (just for testing)
           obj.MatPropDistr.length.fcn = @(D) D;
           
           % Fault thickness
           obj.MatPropDistr.thick = getFaultThickness();
           
           % ResFric
           obj.MatPropDistr.resFric = getResidualFrictionAngle(obj.Vcl);
           
           % SSFc
           obj.MatPropDistr.ssfc = getSSFc(obj.Vcl, obj.IsClayVcl, zf, ...
                                           obj.Thick, Disp, obj.Failure, ...
                                           obj.HW.Id);
                
           % Perm Anisotropy Ratio
           obj.MatPropDistr.permAnisoRatio = ...
                    getAnisotropyRatio(obj.Vcl, zf, clayMine, obj.HW.Id);
                
           % Porosity at current depth
           obj.MatPropDistr.poro = getPorosity(obj.Vcl, obj.IsClayVcl, ...         
                                               zmax, 'zmax', zf, ...
                                               obj.IsUndercompacted, obj.HW.Id);
                                       
           % Perm
           if ~all([isempty(obj.FW.Perm) isempty(obj.HW.Perm)])
               obj.MatPropDistr.perm = getPermeability(obj.Vcl, obj.IsClayVcl, ...
                                                       zf, zmax, obj.MaxPerm, ...
                                                       obj.HW.Id, obj);
           else
               obj.MatPropDistr.perm = getPermeability(obj.Vcl, obj.IsClayVcl, ...
                                                       zf, zmax, obj.MaxPerm, ...
                                                       obj.HW.Id);
           end
            
        end
        
        function plotStrati(obj, faultThick, faultDip, unit, fontSize)
           %
           %
           % This plot considers that dip is constant for all layers in FW
           % and same for the HW (FW and HW dips may be different).
           
           % Multiplier
           if strcmp(unit, 'm')
               m = 1;
           elseif strcmp(unit, 'cm')
               m = 100;
           else
               error('Add corresponding multiplier')
           end
           
           % Measures
           dip = [obj.FW.Dip(1) obj.HW.Dip(1)];
           g = 90 - faultDip;
           if g == 0
               zFWf = [0 cumsum(obj.Tap(obj.FW.Id))];
               zHWf = [0 cumsum(obj.Tap(obj.HW.Id))];
           elseif all(dip == 0)
               zFWf = [0 cumsum(obj.Thick(obj.FW.Id))];
               zHWf = [0 cumsum(obj.Thick(obj.HW.Id))];
           elseif any(dip == 0)
               id = find(dip == 0);
               if id == 1
                   zFWf = [0 cumsum(obj.Thick(obj.FW.Id))];
                   zHWf = [0 cumsum(obj.Tap(obj.HW.Id)*cosd(g))];
               elseif id == 2
                   zFWf = [0 cumsum(obj.Tap(obj.FW.Id)*cosd(g))];
                   zHWf = [0 cumsum(obj.Thick(obj.HW.Id))];
               end
           else
               zFWf = [0 cumsum(obj.Tap(obj.FW.Id)*cosd(g))];
               zHWf = [0 cumsum(obj.Tap(obj.HW.Id)*cosd(g))];
           end
           zFWf = zFWf*m;
           zHWf = zHWf*m;
           xFWf = 0 - zFWf./tand(faultDip);
           faultThick = faultThick*m;
           xHWf = faultThick - zHWf./tand(faultDip);
           fcoord = [0, 0; faultThick, 0; xFWf(end), zFWf(end); ...
                     faultThick + xFWf(end), zFWf(end)];
           limx = [min(fcoord(:, 1)) - 10; max(fcoord(:,1)) + 10];
           zFW = zFWf + tand(dip(1))*(xFWf - limx(1));
           zHW = zHWf + tand(dip(2))*(xHWf - limx(2));
           
           % Utils
           colrs = flipud(copper(16));
           latx = {'Interpreter', 'latex'};
           if nargin < 5
               sz = [14, 12, 11];
               distTextVcl = [2, 8];
           elseif strcmp(fontSize, 'large')
               sz = [20, 20, 20]; 
               distTextVcl = [5, 11];
           end
           
           % Plot
           f = figure(1);
           hold on
           for n=1:numel(obj.FW.Thickness)
              x =  [limx(1), xFWf(n), xFWf(n+1), limx(1), limx(1)];
              z = [zFW(n), zFWf(n), zFWf(n+1), zFW(n+1), zFW(n)];
              if obj.FW.Vcl(n) >= obj.FW.IsClayVcl
                 idc = 1 + round((obj.FW.Vcl(n) - obj.FW.IsClayVcl) / ...
                                  (1 - obj.FW.IsClayVcl)*5);
                 colr = colrs(10+idc, :);
              else
                  idc = 1 + (5 - round((obj.FW.IsClayVcl - obj.FW.Vcl(n)) / ...
                                  (obj.FW.IsClayVcl - 0)*5));
                  colr = colrs(idc, :);
              end
              fill(x, z, colr); 
           end
              
           for n=1:numel(obj.HW.Thickness)
              x =  [limx(2), xHWf(n), xHWf(n+1), limx(2), limx(2)];
              z = [zHW(n), zHWf(n), zHWf(n+1), zHW(n+1), zHW(n)];
              if obj.HW.Vcl(n) >= obj.HW.IsClayVcl
                  idc = 1 + round((obj.HW.Vcl(n) - obj.HW.IsClayVcl) / ...
                                  (1 - obj.HW.IsClayVcl)*5);
                  colr = colrs(10+idc, :);
              else
                  idc = 1 + (5 - round((obj.HW.IsClayVcl - obj.HW.Vcl(n)) / ...
                                  (obj.HW.IsClayVcl - 0)*5));
                  colr = colrs(idc, :);
              end
              fill(x, z, colr); 
           end
           
           for n=1:numel(obj.FW.Thickness)
               if obj.FW.Vcl(n) >= obj.FW.IsClayVcl
                   colrtx = 'w';
               else
                   colrtx = 'k';
               end
               if n == 1
                  str = ['$V_\mathrm{cl}$ = ' num2str(obj.FW.Vcl(n))];
              else
                  str = num2str(obj.FW.Vcl(n));
              end
               text(limx(1) + distTextVcl(1), zFW(n) + (zFW(n+1) - zFW(n))/2, str, latx{:}, ...
                   'fontSize', sz(3), 'color', colrtx)
           end
           for n=1:numel(obj.HW.Thickness)
               if obj.HW.Vcl(n) >= obj.HW.IsClayVcl
                   colrtx = 'w';
               else
                   colrtx = 'k';
               end
               text(limx(2) - distTextVcl(2), zHW(n) + (zHW(n+1) - zHW(n))/2, ...
                   num2str(obj.HW.Vcl(n)), latx{:}, ...
                   'fontSize', sz(3), 'color', colrtx)
           end
           hold off
           %if strcmp(unit, 'm'), axis equal, end
           h = gca; h.XAxis.Visible = 'off';
           h.FontSize = sz(3);
           xlim([min([limx', xFWf, xHWf]), max([limx', xFWf, xHWf])])
           ylim([min([zFW, zHW, zFWf, zHWf]) max([zFW, zHW, zFWf, zHWf])])
           ylabel(['$z$ [' unit ']'], 'fontSize', sz(2), latx{:})
           title(['$z_\mathrm{f}$,  $z_\mathrm{max} =$ ' num2str(unique(obj.DepthFaulting)) ...
                  ', ' num2str(unique(obj.FW.DepthBurial)) ' m'], latx{:}, 'fontSize', sz(1))
           set(f, 'position', [200, 0, 300, 400]);
        end
    end
end

