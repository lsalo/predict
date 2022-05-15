function [Perm] = computeCoarsePerm3D(G, permG, partDims, U)
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
if strcmp(U.method, 'tpfa')
    p2 = partitionCartGrid(G.cartDims, partDims);
    CG2 = generateCoarseGrid(G, p2);
    K = diag(upscalePerm(G, CG2, rock, 'method', U.method));
    
elseif strcmp(U.method, 'mpfa')
    Dp{1} = 5*barsa;
    if U.useAcceleration == 1
        inB = 'mex';
    else
        inB = [];
    end
    hTmp = computeMultiPointTrans(G, rock, 'invertBlocks', inB);
    psolver = @(state0, G, fluid, bc) incompMPFA(state0, G, hTmp, ...
                                                 fluid, 'bc', bc);
    tic
    K = diag(myupscalePermeabilityFixed(G, Dp{1}, psolver, fluid, L));
    toc
%     psolver = @(state0, G, fluid, bc, rock) incompMPFA(state0, G, hTmp, ...
%                                                  fluid, 'bc', bc);
%     K = diag(upscalePermeabilityFixed(G, Dp{1}, psolver, fluid, rock, L));

else
    error("U.method not supported. Choose 'tpfa' or 'mpfa'.")
end

% Pass the 3 independent components   
Perm = [K(1), K(5), K(9)];

end