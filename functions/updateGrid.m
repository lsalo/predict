function G = updateGrid(G, T)
% update grid geometry (volumes, centroids, areas, etc) once the thickness
% has been sampled. This is straightforward because this is a cartesian
% grid and only the thickness is udpdated.
%
% NOTE: The initial G thickness must be 1 for the scaling to be correct.
%

% Checks
id = 1;
if G.griddim == 3
    % Must be layeredGrid (where fault thickness is 2nd dimension!)
    assert(strcmp(G.type{3}, 'makeLayeredGrid'))
    assert(isfield(G, 'xySwap') && G.xySwap)
    id = 2;
    assert(sqrt(G.layerSize)*G.cellDim(id)==1)      
else
    assert(G.cartDims(id)*G.cellDim(id)==1)
end


% cellDim
G.cellDim(id) = T*G.cellDim(id);

% Nodes
G.nodes.coords(:,id) = T*G.nodes.coords(:,id);

% Cells
G.cells.volumes = T*G.cells.volumes;
G.cells.centroids(:,id) = T*G.cells.centroids(:,id);

% Faces
G.faces.centroids(:,id) = T*G.faces.centroids(:,id);
if G.griddim == 3   % only faces inolving the T dimension (y in swap grid)
    facId = G.faces.areas < 0.99*G.cellDim(1)*G.cellDim(3);
    %plotFaces(G,~facId,'faceColor', 'none', 'edgeColor', 'r');
    %xlabel('y [m]'), ylabel('x [m]'); zlabel('z [m]'); view([-60 20])
    G.faces.areas(facId) = T*G.faces.areas(facId);
    G.faces.normals(facId,:) = T*G.faces.normals(facId,:); 
else
    G.faces.areas = T*G.faces.areas;
    G.faces.normals = T*G.faces.normals;
end

end