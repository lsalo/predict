function [Perm, Perm2] = comparePermZDir(G, permG, U, CG)
% Check that reversing the Z axis does not affect result. This can be run
% from within computeCoarsePerm3D to check.

physDim = max(G.faces.centroids);
G2 = cartGrid(G.cartDims, physDim);
G2.nodes.coords(:,end) = -1*(G2.nodes.coords(:,end) - physDim(end));
G2 = computeGeometry(G2);

figure(randi(1000, 1)); subplot(1,2,1)
plotCellData(G, log10(permG(:,1)/(milli*darcy)), 'edgecolor', 'none'); 
view([30 20]); ax=gca; ax.ZDir = 'normal'; ax.DataAspectRatio = [0.05 1 1];
colormap(copper);
subplot(1,2,2)
plotCellData(G2, log10(permG(:,1)/(milli*darcy)), 'edgecolor', 'none'); 
view([30 20]); ax=gca; ax.DataAspectRatio = [0.05 1 1];
colormap(copper);

Perm = myUpscalePerm(G, CG, rock, 'method', U.method);
Perm2 = myUpscalePerm(G2, CG, rock, 'method', U.method);

end