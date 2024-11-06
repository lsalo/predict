function [Perm, Perm_ha, Perm_per] = computeCoarsePerm3D(G, permG, U, CG, compare)
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


if nargin > 4 && compare % based on upscaling examples from MRST book module
    
    if G.griddim == 3
        assert(size(rock.perm, 2) == 6)
        id_ha = [1 4 6];
        id_per = [1 5 9];
    elseif G.griddim == 2
        assert(size(rock.perm, 2) == 3)
        id_ha = [1 3];
        id_per = [1 4];
    end

    % 1. Harmonic-arithmetic
    q = partitionUI(G, U.coarseDims);
    vol = G.cells.volumes;
    Perm_ha = zeros(1, G.griddim);
    for i = 1:3
        dims = G.cartDims; 
        dims(i)=U.coarseDims(i);
        qq = partitionUI(G, dims);
        K = accumarray(qq,vol)./accumarray(qq,vol./rock.perm(:,id_ha(i)));
        Perm_ha(i) = accumarray(q,K(qq).*vol)./accumarray(q,vol);
    end
    
    % 2. Flow-based with periodic boundary condition
    bcr{1}=pside([],G,'RIGHT',0);   bcl{1}=pside([],G,'LEFT',0);
    bcr{2}=pside([],G,'FRONT',0);   bcl{2}=pside([],G,'BACK',0);
    bcr{3}=pside([],G,'BOTTOM',0);  bcl{3}=pside([],G,'TOP',0);
    [Gp, bcp]=makePeriodicGridMulti3d(G, bcl, bcr, {0, 0, 0});
    
    % Do a single-phase periodic upscale
    psolver = @(state, Grid, Fluid, BCP, Rock)...
               incompTPFA(state, Grid, computeTransGp(G, Grid, Rock), Fluid, 'bcp', BCP);
    L = max(G.faces.centroids) - min(G.faces.centroids);
    fluid_pure = initSingleFluid('mu', 1, 'rho', 1);
    warning('off', 'mrst:periodic_bc')
    Perm_per = upscalePermeabilityPeriodic(Gp, bcp, 1, psolver, fluid_pure, rock, L);
    Perm_per = Perm_per(id_per);
    warning('on', 'mrst:periodic_bc')

    
%     % Manual implementation of flow-based with sealing boundaries - to check
%     % 1. Half-trans and fluid
%     hT    = computeTrans(G, rock);
%     fluid = initSingleFluid('mu' ,1, 'rho', 1);
%     
%     % 2. Structures with BC
%     d = G.griddim;
%     [bcl, bcr, Dp]=deal(cell(d,1));
%     bcsides = {'XMin', 'XMax'; 'YMin', 'YMax'; 'ZMin', 'ZMax'};
%     for j = 1:d
%         bcl{j} = pside([], G, bcsides{j, 1}, 0);
%         bcr{j} = pside([], G, bcsides{j, 2}, 0);
%         Dp{j}  = 0;
%     end
%     Dp{1} = 4*barsa;
%     L  = max(G.faces.centroids)-min(G.faces.centroids);
%     
%     % 3. Upscale
%     [v,dp] = deal(zeros(d, 1));
%     for i=1:d
%         bc = addBC([], bcl{i}.face, 'pressure', Dp{1});
%         bc = addBC(bc, bcr{i}.face, 'pressure', Dp{2});
%         
%         xr = initResSol(G, 100*barsa, 1);
%         xr = incompTPFA(xr, G, hT, fluid, 'bc', bc);
%         
%         v(i)  = sum(xr.flux(bcr{i}.face)) / sum(G.faces.areas(bcr{i}.face));
%         dp(i) = Dp{1}/L(i);
%     end
%     K = convertTo(v./dp, milli*darcy);  

end

end