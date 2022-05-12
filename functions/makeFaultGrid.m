function G = makeFaultGrid(T, D, L, targetCellDim, makeplot)
% Generate MRST grid (an MRST installation is required).
%
% INPUTS
%   thick: fault thickness (grid x dimension) [m]
%   disp: fault displacement (grid y (z) dimension) [m]
%   targetCellDim: target cell dimensions (thickness and length). An array
%                  of 1xn elements [m], where d is number of dimensions.
%                  Note that the output grid must be square in [x, z],
%                  so only the targetCellDim leading to higher
%                  resolution (more cells) is respected.
%   makeplot: [optional] pass 1 to generate grid plot coloring cells by
%             index.
%
% OUTPUT
%   G: MRST grid structure (refer to MRST documentation for details). A
%      structured, quadrilateral 2D grid is created. The grid has the same
%      number of cells in x and y dimensions, corresponding to the maximum
%      of {thick/targetCellDim(1), disp/targetCellDim(2)}. Hence, cell size
%      in one of the dimensions may be << than what was passed in
%      targetCellDim.
%
% EXAMPLE
%   thick = 0.5;
%   disp = 100;
%   targetCellDim = [0.1, 1]
%   G = makeFaultGrid(thick, disp, targetCellDim); 
%

% User can pass desired grid resolution
if nargin < 4 || isempty(targetCellDim)
    if nargin < 3 || isempty(L)
        dim = 2;
        targetCellDim = [D/1000, D/100];             % [m], [thick., disp.]
    else
        dim = 3;
        targetCellDim = [D/1000, D/100, D/100];   % [thick, length, disp.]
    end
else
    assert(isa(targetCellDim, 'double') && ...
        numel(targetCellDim) == 2 || numel(targetCellDim) == 3, ...
        ['If grid resolution is passed, it must be a ', ...
        'double array with 2 or 3 elements: ', ...
        '[res. in thick. dir. (x), (res in len. dir. (y)), res. in displ. dir (y/z)]'])
end



nelem = max([round(T / targetCellDim(1)), ...
             round(D / targetCellDim(2))]);
G = cartGrid([nelem, nelem], [T, D]);
G.cellDim = [T/nelem, D/nelem];
if dim == 2
    G = computeGeometry(G);
elseif dim == 3     % extrude 2D cartesian grid
    nelem = max([round(T / targetCellDim(1)), ...
                 round(L / targetCellDim(2)), ...
                 round(D / targetCellDim(3))]);
    G.cellDim = [L/nelem, T/nelem, D/nelem];  
    nLayers = nelem;   % Number of layers
    strikeLayerThickness = G.cellDim(1);
    G = makeLayeredGrid(G, ones(nLayers, 1)*strikeLayerThickness);
    G.nodes.coords(G.nodes.coords(:, 3) == 1, 3) = strikeLayerThickness;
    G.nodes.coords = G.nodes.coords(:, [3 1 2]); % x -> y, y -> z, along-strike = x
    G = computeGeometry(G);
    G.xySwap = true;
end

if nargin > 4 && makeplot == 1
    f1 = figure(1);
    colormap(turbo)
    plotCellData(G, (1:G.cells.num)', 'EdgeColor', [0.2 0.2 0.2], 'EdgeAlpha', 0.1);
    %plotCellData(G, (1:G.layerSize)', (1:G.layerSize)', 'EdgeColor', [0.2 0.2 0.2], ...
    %             'EdgeAlpha', 0.1);
    ax = gca; 
    colorbar;
    %xlim([0 T]); ylim([0 D]);
    title(['Number of cells = ' num2str(G.cells.num)])
    set(f1, 'position', [400, 100, 500, 600]); 
    if dim == 2
        ax.DataAspectRatio = [0.1 1];
        xlabel('x [m]'), ylabel('y [m]');
        view([-30, 75])
    elseif dim == 3
        ax.DataAspectRatio = [1 0.1 1];
        xlabel('y [m]'), ylabel('x [m]'); zlabel('z [m]');
        view([-60 20])
        ax.ZDir = 'normal';
    end
end
end