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
        
        function plotStrati(obj, faultDip)
           %
           %
           %
           
           % Measures
           throw = sum(obj.FW.Thickness);
           fT = 3;      % THIS IS JUST FOR THE PLOT
           fcoord = [0, 0; fT, 0; 0 - throw/tand(faultDip), throw; ...
                     fT - throw/tand(faultDip), throw];
           limx = [min(fcoord(:, 1)) - 10; max(fcoord(:,1)) + 10];
           zFW = [0 cumsum(obj.FW.Thickness)];  
           zHW = [0 cumsum(obj.HW.Thickness)];
           xFWf = 0 - zFW./tand(faultDip);
           xHWf = fT - zHW./tand(faultDip);
           
           % Utils
           colrs = copper(2);
           latx = {'Interpreter', 'latex'};
           sz = [14, 12];
           
           % Plot
           figure(1)
           hold on
           for n=1:numel(obj.FW.Thickness)
              x =  [limx(1), xFWf(n), xFWf(n+1), limx(1), limx(1)];
              z = [zFW(n), zFW(n), zFW(n+1), zFW(n+1), zFW(n)];
              if obj.FW.Vcl(n) >= obj.FW.IsClayVcl
                  colr = colrs(1, :);
              else
                  colr = colrs(2, :);
              end
              fill(x, z, colr); 
           end
           for n=1:numel(obj.HW.Thickness)
              x =  [limx(2), xHWf(n), xHWf(n+1), limx(2), limx(2)];
              z = [zHW(n), zHW(n), zHW(n+1), zHW(n+1), zHW(n)];
              if obj.HW.Vcl(n) >= obj.HW.IsClayVcl
                  colr = colrs(1, :);
              else
                  colr = colrs(2, :);
              end
              fill(x, z, colr); 
           end
           hold off
           axis equal
           h = gca; h.XAxis.Visible = 'off';
           xlim(limx)
           ylim([0 throw])
           ylabel('$z$ [m]', 'fontSize', sz(2), latx{:})
           title(['$z_\mathrm{f} =$ ' num2str(unique(obj.DepthFaulting)) ...
                  ' m $\mid$ $z_\mathrm{max} =$ ' num2str(unique(obj.FW.DepthBurial)) ...
                  ' m'], latx{:}, 'fontSize', sz(1))
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
        
    end
end

