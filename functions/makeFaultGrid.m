function G = makeFaultGrid(thick, disp, targetCellDim, makeplot)
%
%
%
nelem = max([round(thick / targetCellDim(1)), ...
    round(disp / targetCellDim(2))]);
G = computeGeometry(cartGrid([nelem, nelem], [thick, disp]));
G.CellDim = [thick/nelem, disp/nelem];

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