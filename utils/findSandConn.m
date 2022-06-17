function [c, M] = findSandConn(M, method, dim, G)
% From a binary matrix with 0 for sand and 1 for clay, find whether any
% object (group of connected sand pixels) connects the left boundary to the
% right boundary (x dimension) and/or the top boundary with the bottom
% boundary (z dimension).
%
% Note that here, "object" is pixel-defined in a binary image. Any single
% sand object here can be the combination of multiple sand layers and types
% in the fault.
%
% INPUTS:
%   M: binary matrix where each entry is 0 (sand) or 1 (clay)
%   method: method to approximate fluxes on the cell interfaces during
%           permeability upscaling. Can be either 'tpfa' or 'mpfa'. If
%           'tpfa' is passed, then the pixel connectivity is 4-connected,
%           otherwise is 8-connected.
%           If no method is passed, the default (8-connected) is used.
%
% OUTPUT:
%   c: structure with fields:
%       x: bounding box x dimension (in pixel units = number of matrix
%          entries) for each object.
%       z: same as x in z dimension.
%       bc: 1x2 logical array, first entry true if L and R boundaries are
%       connected by a single object, second entry true if top and bottom
%       boundaries are connected by a single object.
%
% REQUIREMENTS:
%   Image Processing Toolbox
%
if nargin < 3 || dim == 2
    assert(ismatrix(M));        % check input is matrix
elseif dim == 3                 % reshape to 3D array for connectivity
    M_grid = M;
    [nx, ny, nz] = deal(G.cartDims(1), G.cartDims(2), G.cartDims(3));
    M = zeros(nx,nz,ny);
    layerSize = nx*nz;
    id_layer1 = repmat(1:nx, nz, 1);
    id_layer1(2:end, :) = id_layer1(2:end, :) + (1:nz-1)'*(nx*ny);
    id_layer1 = reshape(id_layer1', layerSize, 1);
    %spy(flipud(transpose(reshape(M_grid(id_layer1), nx, nz))));
    id_layers = [id_layer1 repmat(id_layer1, 1, ny-1) + (1:ny-1)*nx];
    id_layers = reshape(id_layers, nx*ny*nz, 1);
    M(:) = M_grid(id_layers);
    M = flipud(pagetranspose(M));
    if ismatrix(M)  % one layer only, we duplicate layer to make 3D
        M(:,:,2) = M;
    end
end
M = 1 - M;                  % we want 0 = clay, 1 = sand

% Set pixel connectivity
conn = 8;
if nargin == 2 && strcmp(method, 'tpfa') || ...
   dim == 2 && strcmp(method, 'tpfa')
    conn = 4;     
elseif dim == 3 && strcmp(method, 'tpfa')
    conn = 6;
end

% Get bounding box
L = labelmatrix(bwconncomp(M, conn));
B = regionprops(L, 'BoundingBox');
%imshow(label2rgb(repelem(L, 10, 1),'jet'), 'InitialMagnification','fit');
%imshow(L, 'InitialMagnification', 'fit');
%axis on

% Find extension of each object
sz = size(M);
if nargin < 3 || dim == 2
    c.x = arrayfun(@(s) s.BoundingBox(3), B)';
    c.z = arrayfun(@(s) s.BoundingBox(4), B)';
    c.bc = false(1, 2);
    if any(c.x == sz(1)), c.bc(1) = true; end
    if any(c.z == sz(2)), c.bc(2) = true; end
elseif dim == 3
    c.x = arrayfun(@(s) s.BoundingBox(4), B)';
    c.z = arrayfun(@(s) s.BoundingBox(5), B)'; 
    c.y = arrayfun(@(s) s.BoundingBox(6), B)';
    c.bc = false(1, 3);
    if any(c.x == sz(1)), c.bc(1) = true; end
    if any(c.z == sz(2)), c.bc(2) = true; end
    if any(c.y == sz(3)), c.bc(3) = true; end
end


end