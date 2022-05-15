function extrudedPerm = assignExtrudedPerm(G, extrudedPerm, faultSection, ...
                                           segLen, cellDim)
%
%
%

% Initial checks
assert(mod(segLen/cellDim, 1)==0)
nklayers = segLen/cellDim;              % repeat 2D section for nklayers

% cartGrid
kxx = faultSection.Grid.perm(:,1);
kxz = faultSection.Grid.perm(:,2);
kyy = faultSection.Grid.permy;
kzz = faultSection.Grid.perm(:,3);
idPerm3D = [1 3 4 6];                       % rotated about the y axis (along-strike)
[nx, ny, nz] = deal(G.cartDims(1), G.cartDims(2), G.cartDims(3));
idLayer1 = repmat((1:nx)', nz, 1);
idLayer1 = idLayer1 + repelem(nx*ny*(0:nz-1)', nx, 1);
idLayers = repmat(idLayer1, 1, nklayers);
idLayers(:,2:end) = idLayers(:,2:end) + ...
              repelem(repelem(nx, numel(idLayer1), 1), 1, nklayers-1) + ...
              ((nx:nx:(nklayers-1)*nx) - nx);
idLayers = reshape(idLayers, nx*nz*nklayers, 1);
extrudedPerm(idLayers, idPerm3D) = repmat([kxx, kyy, kyz, kzz], nklayers, 1); 
% extruded
% idFirst0 = find(extrudedPerm==0, 1); 
% kyy = faultSection.Grid.perm(:,1);
% kyz = faultSection.Grid.perm(:,2);
% kzz = faultSection.Grid.perm(:,3);
% kxx = faultSection.Grid.permy;
% idPerm3D = [1 4 5 6];                       % rotated about the x axis (along-strike)
% layerSize = size(kyy, 1);
% extrudedPerm(idFirst0:(idFirst0-1)+nklayers*layerSize, idPerm3D) = ...
%             repmat([kxx, kyy, kyz, kzz], nklayers, 1); 

end
     