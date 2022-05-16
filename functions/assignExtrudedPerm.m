function [extrudedPerm, isSmear] = assignExtrudedPerm(G, extrudedPerm, isSmear, ...
                                                      faultSection, ...
                                                      segLen, cellDim)
%
%
%

% Initial checks
if ~isnan(G.cellDim(2))
    assert(mod(segLen/cellDim, 1)==0)
    nklayers = segLen/cellDim;              % repeat 2D section for nklayers
else
    nklayers = 1;
end
if isempty(extrudedPerm)
    extrudedPerm = zeros(G.cells.num, 6);   %[kxx, kxy, kxz, kyy, kyz, kzz]
    isSmear = false(G.cells.num, 1);
end
idFirst0 = find(all(extrudedPerm==0, 2), 1); 

% cartGrid
kxx = faultSection.Grid.perm(:,1);
kxz = faultSection.Grid.perm(:,2);
kyy = faultSection.Grid.permy;
kzz = faultSection.Grid.perm(:,3);
idPerm3D = [1 3 4 6];                       % rotated about the y axis (along-strike)
[nx, ny, nz] = deal(G.cartDims(1), G.cartDims(2), G.cartDims(3));
idCellsLayer1 = repmat((1:nx)', nz, 1) + idFirst0-1;
idCellsLayer1 = idCellsLayer1 + repelem(nx*ny*(0:nz-1)', nx, 1);
idCellsLayers = repmat(idCellsLayer1, 1, nklayers);
if nklayers > 1
    idCellsLayers(:,2:end) = idCellsLayers(:,2:end) + ...
                  repelem(repelem(nx, numel(idCellsLayer1), 1), 1, nklayers-1) + ...
                  ((nx:nx:(nklayers-1)*nx) - nx);
    idCellsLayers = reshape(idCellsLayers, nx*nz*nklayers, 1);
end
extrudedPerm(idCellsLayers, idPerm3D) = repmat([kxx, kxz, kyy, kzz], nklayers, 1); 

% extruded
% kyy = faultSection.Grid.perm(:,1);
% kyz = faultSection.Grid.perm(:,2);
% kzz = faultSection.Grid.perm(:,3);
% kxx = faultSection.Grid.permy;
% idPerm3D = [1 4 5 6];                       % rotated about the x axis (along-strike)
% layerSize = size(kyy, 1);
% extrudedPerm(idFirst0:(idFirst0-1)+nklayers*layerSize, idPerm3D) = ...
%             repmat([kxx, kyy, kyz, kzz], nklayers, 1); 

% isSmear
M = faultSection.MatMap;
nCellsLayer = G.cartDims(1)*G.cartDims(3);
isSmear(idCellsLayers) = repmat(reshape(transpose(flipud(M.vals)), ...
                                        nCellsLayer, 1), nklayers, 1);

end
     