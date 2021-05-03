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
    end
    
    properties (Dependent, Hidden)
       ParentId 
       Vcl
       IsClayVcl
       DepthFaulting
    end
    
    methods
        function obj = FaultedSection(footwall, hangingwall, faultDip)
            % We instantiate the object with the FW and HW objects required
            % to construct it.
            %
            % Example: mySect = FaultedSection(footwall, hangingwall)
            %
            
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
            assert(sum(TapFW) == sum(TapHW))
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
        
        function obj = getMatPropDistr(obj, varargin)
           %
           % TBD
           %
           
           % Optional inputs
           opt.maxPerm = [];                   % mD
           opt.siltInClay = false;
           opt.isUndercompacted = false;
           opt = merge_options_relaxed(opt, varargin{:});
            
           % Shorten input names
           zf   = obj.DepthFaulting;
           zmax = [obj.FW.DepthBurial, obj.HW.DepthBurial];
           clayMine = {obj.FW.ClayMine, obj.HW.ClayMine};
           Disp = sum(obj.Tap(obj.FW.Id));
           
           % Fault thickness
           obj.MatPropDistr.thick = getFaultThickness();
           
           % ResFric
           obj.MatProps.ResFric = getResidualFrictionAngle(obj.Vcl, FS);
           
           % SSFc and SSFc bounds
           obj.MatProps.SSFc = getSSFc(obj.Vcl, obj.IsClayVcl, zf, ...
                                       obj.Thick, Disp, obj.HW.Id);
            
            
        end
        
        function plotStrati(obj, faultThick, faultDip)
           %
           %
           % This plot considers that dip is constant for all layers in FW
           % and same for the HW (FW and HW dips may be different).
           
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
           xFWf = 0 - zFWf./tand(faultDip);
           xHWf = faultThick - zHWf./tand(faultDip);
           fcoord = [0, 0; faultThick, 0; xFWf(end), zFWf(end); ...
                     faultThick + xFWf(end), zFWf(end)];
           limx = [min(fcoord(:, 1)) - 10; max(fcoord(:,1)) + 10];
           zFW = zFWf + tand(dip(1))*(xFWf - limx(1));
           zHW = zHWf + tand(dip(2))*(xHWf - limx(2));
           
           % Utils
           colrs = flipud(copper(16));
           latx = {'Interpreter', 'latex'};
           sz = [14, 12];
           
           % Plot
           figure(1)
           hold on
           for n=1:numel(obj.FW.Thickness)
              x =  [limx(1), xFWf(n), xFWf(n+1), limx(1), limx(1)];
              z = [zFW(n), zFWf(n), zFWf(n+1), zFW(n+1), zFW(n)];
              if obj.FW.Vcl(n) >= obj.FW.IsClayVcl
                 idc = 1 + round((obj.FW.Vcl(n) - obj.FW.IsClayVcl) / ...
                                  (1 - obj.FW.IsClayVcl)*5);
                 colr = colrs(10+idc, :);
                 colrtx = 'w';
              else
                  idc = 1 + (5 - round((obj.FW.IsClayVcl - obj.FW.Vcl(n)) / ...
                                  (obj.FW.IsClayVcl - 0)*5));
                  colr = colrs(idc, :);
                  colrtx = 'k';
              end
              fill(x, z, colr); 
              if n == 1
                  str = ['$V_\mathrm{cl}$ = ' num2str(obj.FW.Vcl(n))];
              else
                  str = num2str(obj.FW.Vcl(n));
              end
              text(limx(1) + 3, zFW(n) + (zFW(n+1) - zFW(n))/2, str, latx{:}, ...
                   'fontSize', 11, 'color', colrtx)
           end
           for n=1:numel(obj.HW.Thickness)
              x =  [limx(2), xHWf(n), xHWf(n+1), limx(2), limx(2)];
              z = [zHW(n), zHWf(n), zHWf(n+1), zHW(n+1), zHW(n)];
              if obj.HW.Vcl(n) >= obj.HW.IsClayVcl
                  idc = 1 + round((obj.HW.Vcl(n) - obj.HW.IsClayVcl) / ...
                                  (1 - obj.HW.IsClayVcl)*5);
                  colr = colrs(10+idc, :);
                  colrtx = 'w';
              else
                  idc = 1 + (5 - round((obj.HW.IsClayVcl - obj.HW.Vcl(n)) / ...
                                  (obj.HW.IsClayVcl - 0)*5));
                  colr = colrs(idc, :);
                  colrtx = 'k';
              end
              fill(x, z, colr); 
              text(limx(2) - 10, zHW(n) + (zHW(n+1) - zHW(n))/2, ...
                   num2str(obj.HW.Vcl(n)), latx{:}, ...
                   'fontSize', 11, 'color', colrtx)
           end
           hold off
           axis equal
           h = gca; h.XAxis.Visible = 'off';
           xlim([min([limx', xFWf, xHWf]), max([limx', xFWf, xHWf])])
           ylim([min([zFW, zHW, zFWf, zHWf]) max([zFW, zHW, zFWf, zHWf])])
           ylabel('$z$ [m]', 'fontSize', sz(2), latx{:})
           title(['$z_\mathrm{f} =$ ' num2str(unique(obj.DepthFaulting)) ...
                  ' m $\mid$ $z_\mathrm{max} =$ ' num2str(unique(obj.FW.DepthBurial)) ...
                  ' m'], latx{:}, 'fontSize', sz(1))
        end
    end
end

