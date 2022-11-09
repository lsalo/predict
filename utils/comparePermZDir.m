function [Perm, Perm2, res] = comparePermZDir(G, permG, U, CG)
% Check that reversing the Z axis does not affect result. This can be run
% from within computeCoarsePerm3D to check.

physDim = max(G.faces.centroids);
G2 = cartGrid(G.cartDims, physDim);
G2 = computeGeometry(G2);

idG = reshape(1:G.cells.num, G.cartDims(1)*G.cartDims(2), G.cartDims(3));
idG = reshape(fliplr(idG), G.cells.num, 1);
rock.perm = permG;
rock2.perm = permG(idG, :);

figure(randi(1000, 1)); subplot(1,2,1)
plotCellData(G, log10(rock.perm(:,1)/(milli*darcy)), 'edgecolor', 'none'); 
view([30 20]); ax=gca; ax.ZDir = 'normal'; ax.DataAspectRatio = [0.1 1 1];
colormap(copper);
subplot(1,2,2)
plotCellData(G2, log10(rock2.perm(:,1)/(milli*darcy)), 'edgecolor', 'none'); 
view([30 20]); ax=gca; ax.DataAspectRatio = [0.1 1 1];
colormap(copper);

Perm = myUpscalePerm(G, CG, rock, 'method', U.method);
Perm2 = myUpscalePerm(G2, CG, rock2, 'method', U.method);
res = abs(Perm2 - Perm);


end