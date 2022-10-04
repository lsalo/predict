function Perm = computeCoarsePerm3D(G, permG, U, CG)
%
% SUMMARY
% Obtain the permeability of the coarse grid equivalent to that of the
% fine grid, i.e. upscale the fine-grid permeability. This is accomplished
% using flow-based upscaling. Flow-based upscaling consists in imposing an
% axial p-drop in each direction, while maintaining the other two
% boundaries sealed (similar to laboratory experiments).
%
% INPUTS: 
%   G:      fault grid (MRST grid, fine-scale)
%   permG:  fine-grid anisotropic permeability, typically ncells x 6
%   U:      upscaling options (U.method 'tpfa' or 'mpfa', however 'mpfa' is
%                              restricted to certain conditions, see below)
%   CG:     coarse grid (MRST grid, based on partition of G)
%   * see examples for workflow, and Fault3D for input definition
%
%  OUTPUT: 
%   Perm =  upscaled permeability (n coarse cells x 3). First column is
%           kxx, second is kyy and third is kzz.
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
    assert(~isfield(U, 'flexible') || U.flexible, ...
           "mpfa too slow for fine-scale grids with full along-strike resolution")
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

end