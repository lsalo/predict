classdef Stratigraphy       % Value class to define a data structure.
    %
    % SUMMARY:
    %   Define a stratigraphy object with appropriate fields, such
    %   as number of layers, thickness of each layer, clay content of each
    %   layer, and so on. 
    %   For each input array, the first index corresponds to the 
    %   deepest layer, and the last one to the shallowest.
    %
    % 
    % REQUIRED PARAMETERS:
    %   thickness:  Thickness of each layer on that side of the 
    %               fault (hangingwall or footwall). Must be a double array 
    %               of size 1xN, where N is the number of layers. In meters.
    %
    %   vcl:        Clay volume fraction of each layer. Must be a double   
    %               array of size 1xN, where N is the number of layers.
    %
    %   dip:        Dip angle of each layer, in degrees. The dip angle is 
    %               the angle between the horizontal and the line of max
    %               slope on a given layer. It can take an array of size
    %               1xN, where N is the number of layers, or a scalar value
    %               (same dip angle for all layers).
    %
    %
    % RECOMMENDED PARAMETERS:
    %   'DepthFaulting':  Faulting depth, in meters.
    %   'DepthBurial':    Maximum burial depth, in meters.
    %
    %
    % OPTIONAL PARAMETERS:
    %   'property' - Set property to the specified value.
    %
    %
    % RETURNS:
    %   Class instance.
    %
    %
    % EXAMPLES:
    %   footwall = Stratigraphy(thickness, vcl, dip)
    %   footwall = Stratigraphy(thickness, vcl, dip, 'DepthFaulting', 50)
    %   footwall = Stratigraphy(thickness, vcl, dip, ...
    %                           'DepthBurial', 2000, 'ResFric', [15, 45])
    %   hangingwall = Stratigraphy(thickness, vcl, dip, 'IsHW', 1, ...
    %                              'NumLayersFW', footwall.NumLayers)
    % 
    % ____________________________________________________________________
    
    properties
        % Required as inputs to instantiate the class.
        Thickness       % Thickness of each layer [m]
        Vcl             % Clay volume fraction [-]
        Dip             % Dip of each layer (w/ respect to hzntal) [deg]
        
        % Optional
        DepthFaulting   % Depth of faulting [m]
        DepthBurial     % Maximum burial depth of the layer [m]
        ClayMine        % Predominant clay mineral in clay beds
        IsHW            % pass 1 if this is the HangingWall object.
        Perm            % Permeability perpendicular to bedding [mD]
        ResFric         % Residual Friction Angle [deg]
    end
    properties (Hidden) % Required if HW is being defined
        NumLayersFW     % if IsHW = 1, pass obj.NumLayers of FW object.
        IsClayVcl       % clay smears for Vcl >= this value.
    end  
    properties (Dependent, Hidden)    % Computed based on other inputs
        NumLayers       % Number of layers
        Id              % 1: NumLayers (FW), N(FW)+1:N(FW)+NumLayers(HW) (HW)
        IsClay          % true for Vcl >= IsClayValue   
        
    end
    
    methods
        function obj = Stratigraphy(thickness, vcl, dip, varargin)
            %  Construct an instance of this class with properties

            % Required inputs
            obj.Thickness = thickness;
            obj.Vcl       = vcl;
            if numel(dip) == 1
                obj.Dip = repelem(dip, numel(vcl));
            else
                obj.Dip = dip;
            end
            
            % Optional inputs
            obj.DepthFaulting = 50;    % default (used if not provided)
            obj.DepthBurial   = repelem(2000, numel(vcl));   % "
            obj.ClayMine      = 'kao';  % (kaolinite) "    
            obj.IsHW          = 0;      % "
            obj = merge_options_relaxed(obj, varargin{:});
            
            % Predetermined properties
            % IsClayVcl: consider smear potential at Vcl >= 0.4 
            %            (Fisher & Knipe, 1998)
            obj.IsClayVcl = 0.4;      
            
            % Assert sizes
            len = numel(thickness);
            assert(numel(vcl) == len, 'vcl must be same size as thickness.')
            assert(numel(obj.DepthBurial) == len, ...
                   'DepthBurial must be same size as thickness')  
            assert(numel(dip) == 1 || numel(dip) == len, ...
                   strcat('dip must be an array of size 1xN, where N', ...
                          ' is the number of layers, or a scalar value', ...
                          ' (same dip angle for all layers).'))
            assert(numel(obj.DepthFaulting) == 1, ...
                   'DepthFaulting must be a scalar.')
            assert(isempty(obj.Perm) || numel(obj.Perm) == len, ...
                   strcat('if Perm is provided, it must be an array', ...
                          ' of size 1xN, where N is the number of', ...
                          ' layers. NaN entries are allowed.'))
            assert(isempty(obj.ResFric) || numel(obj.ResFric) == len, ...
                   strcat('if ResFric is provided, it must be an array', ...
                          ' of size 1xN, where N is the number of', ...
                          ' layers. NaN entries are allowed.'))
               
        end
        
        
        % Calculate values of dependent properties
        function id = get.Id(obj)
            id = 1:numel(obj.Thickness);
            if obj.IsHW == 1
                id = id + obj.NumLayersFW;                
            end   
        end       
        
        function numLayers = get.NumLayers(obj)
            numLayers = numel(obj.Thickness);
        end
        
        function isClay = get.IsClay(obj)
            isClay = false(1, obj.NumLayers);
            isClay(obj.Vcl >= obj.IsClayVcl) = true;
        end
        
    end
end

