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
        
    end
end

