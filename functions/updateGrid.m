function G = updateGrid(G, T)
% update grid geometry (volumes, centroids, areas, etc) once the thickness
% has been sampled. This is straightforward because this is a cartesian
% grid and only the thickness is udpdated.
%
% NOTE: The initial G thickness must be 1 for the scaling to be correct.
%

% Checks
id = 1;
assert(G.cartDims(id)*G.cellDim(id)==1)
% if G.griddim == 3
%     % Must be layeredGrid (where fault thickness is 2nd dimension!)
%     assert(strcmp(G.type{3}, 'makeLayeredGrid'))
%     assert(isfield(G, 'xySwap') && G.xySwap)
%     id = 2;
%     assert(sqrt(G.layerSize)*G.cellDim(id)==1)   
%     assert(~isfield(G, 'cartDims'))
%     assert(mod(sqrt(G.layerSize), 1) == 0)
%     G.cartDims = [G.numLayers sqrt(G.layerSize) sqrt(G.layerSize)];
% else
%     assert(G.cartDims(id)*G.cellDim(id)==1)
% end

% cellDim
G.cellDim(id) = T*G.cellDim(id);

% Nodes
G.nodes.coords(:,id) = T*G.nodes.coords(:,id);

% Cells
G.cells.volumes = T*G.cells.volumes;
G.cells.centroids(:,id) = T*G.cells.centroids(:,id);

% Faces
G.faces.centroids(:,id) = T*G.faces.centroids(:,id);
if G.griddim == 3   % only faces inolving the T dimension
    facId = G.faces.areas < 0.99*G.cellDim(2)*G.cellDim(3);
    %plotFaces(G,~facId,'faceColor', 'none', 'edgeColor', 'r');
    %xlabel('x [m]'), ylabel('y [m]'); zlabel('z [m]'); view([30 20])
    G.faces.areas(facId) = T*G.faces.areas(facId);
    G.faces.normals(facId,:) = T*G.faces.normals(facId,:); 
else
    G.faces.areas = T*G.faces.areas;
    G.faces.normals = T*G.faces.normals;
end

end