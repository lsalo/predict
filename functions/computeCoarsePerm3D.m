function Perm = computeCoarsePerm3D(G, permG, U, CG)
%
% SUMMARY
% Obtain the permeability of a single-cell grid equivalent to that of the
% fine grid, i.e. upscale the fine-grid permeability. This is accomplished
% using flow-based upscaling. Flow-based upscaling consists in imposing an
% axial p-drop in each direction, while maintaining the other two
% boundaries sealed.
% 

% Initial variables
L = max(G.faces.centroids) - min(G.faces.centroids);
fluid = initSingleFluid('mu' , 1, 'rho', 1);
rock.perm = permG;

% Compute equivalent/upscaled perm according to input method
gravity reset off
if strcmp(U.method, 'tpfa')
    %CG = generateCoarseGrid(G, p);
    Perm = myUpscalePerm(G, CG, rock, 'method', U.method);
    
elseif strcmp(U.method, 'mpfa')
    assert(all(U.coarseDims==1), ...
           "mpfa not implemented for upscaled grids with ncells > 1")
    Dp{1} = 5*barsa;
    if U.useAcceleration == 1
        inB = 'mex';
    else
        inB = [];
    end
    hTmp = computeMultiPointTrans(G, rock, 'invertBlocks', inB);
    psolver = @(state0, G, fluid, bc) incompMPFA(state0, G, hTmp, ...
                                                 fluid, 'bc', bc);
    Perm = myupscalePermeabilityFixed(G, Dp{1}, psolver, fluid, L);
%     psolver = @(state0, G, fluid, bc, rock) incompMPFA(state0, G, hTmp, ...
%                                                  fluid, 'bc', bc);
%     K = diag(upscalePermeabilityFixed(G, Dp{1}, psolver, fluid, rock, L));

else
    error("U.method not supported. Choose 'tpfa' or 'mpfa'.")
end

% Plot
% cmap = copper;
% f2 = figure(2);
% colormap(cmap);
% plotCellData(G, log10(extrudedPerm(:,1)/(milli*darcy)), 'EdgeColor', [0.2 0.2 0.2], 'EdgeAlpha', 0.1)
% ax = gca;
% colorbar;
% %xlim([0 T]); ylim([0 D]);
% title(['Number of cells = ' num2str(G.cells.num)])
% set(f2, 'position', [400, 100, 600, 600]);
% ax.DataAspectRatio = [0.2 1 1];
% xlabel('x [m]'), ylabel('y [m]'); zlabel('z [m]');
% view([30 20])
% ax.ZDir = 'normal';
% 
% f3 = figure(3);
% subplot(1,2,1)
% plotCellData(G, log10(extrudedPerm(isSmear,1)/(milli*darcy)), isSmear, ...
%              'EdgeColor', [0.8 0.8 0.8], 'EdgeAlpha', 0.1)
% ax = gca;
% colormap(ax, cmap(1:100, :));
% colorbar;
% %xlim([0 T]); ylim([0 D]);
% title(['Clay smear perm'])
% set(f3, 'position', [100, 100, 500, 600]);
% ax.DataAspectRatio = [0.2 1 1];
% xlabel('x [m]'), ylabel('y [m]'); zlabel('z [m]');
% view([30 20])
% ax.ZDir = 'normal';
% grid on
% 
% subplot(1,2,2)
% plotCellData(G, log10(extrudedPerm(~isSmear,1)/(milli*darcy)), ~isSmear, ...
%              'EdgeColor', [0.2 0.2 0.2], 'EdgeAlpha', 0.1)
% ax = gca;
% colormap(ax, cmap(156:end,:));
% colorbar;
% %xlim([0 T]); ylim([0 D]);
% title(['Sand smear perm'])
% set(f3, 'position', [100, 100, 900, 500]);
% ax.DataAspectRatio = [0.2 1 1];
% %xlabel('x [m]'), ylabel('y [m]'); zlabel('z [m]');
% view([30 20])
% ax.ZDir = 'normal';
% grid on

end