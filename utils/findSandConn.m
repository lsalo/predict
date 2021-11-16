function c = findSandConn(M, method)
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
assert(ismatrix(M));        % check input is matrix
M = 1 - M;                  % we want 0 = clay, 1 = sand

% Set pixel connectivity
conn2D = 8;
if nargin > 1 && strcmp(method, 'tpfa')
    conn2D = 4;     
end

% Get bounding box
L = labelmatrix(bwconncomp(M, conn2D));
B = regionprops(L, 'BoundingBox');
imshow(label2rgb(L,'jet'), 'InitialMagnification','fit');

% Find extension of each object
dim = size(M);
c.x = arrayfun(@(s) s.BoundingBox(3), B)';
c.z = arrayfun(@(s) s.BoundingBox(4), B)';
c.bc = false(1, 2);
if any(c.x == dim(1)), c.bc(1) = true; end
if any(c.z == dim(2)), c.bc(2) = true; end


end