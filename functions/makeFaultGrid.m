function G = makeFaultGrid(thick, disp, targetCellDim, makeplot)
% Generate MRST grid (an MRST installation is required).
%
% INPUTS
%   thick: fault thickness (grid x dimension) [m]
%   disp: fault displacement (grid y (z) dimension) [m]
%   targetCellDim: target cell dimensions (thickness and length). An array
%                  of 1x2 elements [m]. Note that the output grid must be 
%                  square, so only the targetCellDim leading to higher
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
nelem = max([round(thick / targetCellDim(1)), ...
    round(disp / targetCellDim(2))]);
G = computeGeometry(cartGrid([nelem, nelem], [thick, disp]));
G.cellDim = [thick/nelem, disp/nelem];

if nargin > 3 && makeplot == 1
    f1 = figure(1);
    plotToolbar(G, G, 'EdgeColor', [0.2 0.2 0.2], ...
        'EdgeAlpha', 0.1);
    axis equal; colorbar;
    xlim([0 thick]); ylim([0 disp]);
    xlabel('x [m]'), ylabel('y [m]')
    title(['Number of cells = ' num2str(G.cells.num)])
    set(f1, 'position', [400, 0, 350, 900]); view([-30, 75])
end
end