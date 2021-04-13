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
        HW              % Hangingwall object from Stratigraphy class.
        FW              % Footwall object from Stratigraphy class.
    end
    
    properties (Dependent, Hidden)
       ParentId 
       Vcl
       IsClayVcl
       Thick
       DepthFaulting
    end
    
    methods
        function obj = FaultedSection(footwall, hangingwall)
            % We instantiate the object with the FW and HW objects required
            % to construct it.
            %
            % Example: mySect = FaultedSection(footwall, hangingwall)
            %
            
            % Checks
            assert(hangingwall.IsHW == 1)
            assert(sum(footwall.Thickness) == sum(hangingwall.Thickness))
            assert(footwall.IsClayVcl == hangingwall.IsClayVcl);  % assumed by MatProps
            
            % Assign
            obj.FW = footwall;
            obj.HW = hangingwall;
        end
        
        function parentId = get.ParentId(obj)
            % get parent Ids, ie the corresponding protolith Ids for our 
            % fault materials.
           parentId = [obj.FW.Id obj.HW.Id]; 
        end
        
        function vcl = get.Vcl(obj)
            % get parent Ids, ie the corresponding protolith Ids for our 
            % fault materials.
           vcl = [obj.FW.Vcl obj.HW.Vcl]; 
        end
        
        function isClayVcl = get.IsClayVcl(obj)
            % get parent Ids, ie the corresponding protolith Ids for our 
            % fault materials.
           assert(obj.FW.IsClayVcl == obj.HW.IsClayVcl)
           isClayVcl = obj.FW.IsClayVcl; 
        end
        
        function thick = get.Thick(obj)
            % get parent Ids, ie the corresponding protolith Ids for our 
            % fault materials.
           thick = [obj.FW.Thickness obj.HW.Thickness]; 
        end
        
        function zf = get.DepthFaulting(obj)
            % get parent Ids, ie the corresponding protolith Ids for our 
            % fault materials.
           zf = [obj.FW.DepthFaulting obj.HW.DepthFaulting]; 
        end
        
        function Tap = getApparentThick(FS, faultDip)
            %
            % get apparent thickness of layers at fault cutoffs
            %
            Tap = [FS.FW.Thickness ./ (cosd(FS.FW.Dip)*sind(faultDip)), ...
                   FS.HW.Thickness ./ (cosd(FS.HW.Dip)*sind(faultDip))];
            
        end
        
        function plotStrati(obj, faultDip, layerDip)
           %
           %
           %
           
           % Measures
           throw = sum(obj.FW.Thickness);
           fT = 3;      % THIS IS JUST FOR THE PLOT
           fcoord = [0, 0; fT, 0; 0 - throw/tand(faultDip), throw; ...
                     fT - throw/tand(faultDip), throw];
           limx = [min(fcoord(:, 1)) - 10; max(fcoord(:,1)) + 10];
           if layerDip(1) >= 0
               zFWf = [0 cumsum(obj.FW.Thickness.*cosd(layerDip(1)))];
           else
               zFWf = [0 cumsum(obj.FW.Thickness./cosd(layerDip(1)))];
           end
           if layerDip(2) >= 0
               zHWf = [0 cumsum(obj.HW.Thickness.*cosd(layerDip(2)))];
           else
               zHWf = [0 cumsum(obj.HW.Thickness./cosd(layerDip(2)))];
           end
           xFWf = 0 - zFWf./tand(faultDip);
           xHWf = fT - zHWf./tand(faultDip);
           zFW = zFWf - tand(layerDip(1))*(xFWf - limx(1));
           zHW = zHWf - tand(layerDip(2))*(xHWf - limx(2));
           
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
           xlim(limx)
           ylim([0 max([zFW, zHW])])
           ylabel('$z$ [m]', 'fontSize', sz(2), latx{:})
           title(['$z_\mathrm{f} =$ ' num2str(unique(obj.DepthFaulting)) ...
                  ' m $\mid$ $z_\mathrm{max} =$ ' num2str(unique(obj.FW.DepthBurial)) ...
                  ' m'], latx{:}, 'fontSize', sz(1))
        end
    end
end

