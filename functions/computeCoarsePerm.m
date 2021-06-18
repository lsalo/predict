function [Perm, axialDropVals, K2] = computeCoarsePerm(G, permG, ...
                                                       kz_loc, U, ...
                                                       faultDisp, ...
                                                       faultThick)
%
% SUMMARY
% Obtain the permeability of a single-cell grid equivalent to that of the
% fine grid, i.e. upscale the fine-grid permeability. This is accomplished
% using flow-based upscaling for k_xx and k_zz, while an area-weighted
% average is used for k_yy. Flow-based upscaling consists in imposing an
% axial p-drop in each direction, while maintaining the other two
% boundaries sealed.
% Additionally, this function allows (1) comparing the outflux obtained 
% between the fine grid with detailed permeability and the coarse grid with 
% upscaled permeability, and (2) checking if the grid aspect ratio is too
% high for accurate flow-based upscaling.
%
% INPUTS
%   G: MRST grid structure obtained using makeFaultGrid
%   permG: nx3 array of full (symmetric) tensor permeability. n is the
%          number of grid cells in G, columns are xx, xz, and zz.
%   kz_loc: nx1 array with permeability along the fault materials (before
%           tensor transformation applied in setGridPoroPerm.m)
%   U: structure with options for flow-based upscaling.
%   faultDisp: fault displacement
%   faultThick: fault thickness
%
% OUTPUTS
%   Perm: 1x3 array with upscaled permeabilities [k_xx, k_yy, k_zz]
%   axialDropVals: fluxes at the outflow boundaries 
%                  [fine grid, coarse (upscaling) grid] 
%                  (U.outflux must be set to 1).
%   K2: Upscaled permeability obtained with a grid with low cell aspect
%       ratio. U.ARcheck must be set to 1.
% 

% Initial variables
L = max(G.faces.centroids) - min(G.faces.centroids);
fluid = initSingleFluid('mu' , 1, 'rho', 1);
rock.perm = permG;

% Compute equivalent/upscaled perm according to input method
if strcmp(U.method, 'tpfa')
    p2 = partitionCartGrid(G.cartDims, [1 1]);
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
    K = diag(myupscalePermeabilityFixed(G, Dp{1}, psolver, fluid, L));
%     psolver = @(state0, G, fluid, bc, rock) incompMPFA(state0, G, hTmp, ...
%                                                  fluid, 'bc', bc);
%     K = diag(upscalePermeabilityFixed(G, Dp{1}, psolver, fluid, rock, L));

else
    error("U.method not supported. Choose 'tpfa' or 'mpfa'.")
end

% Compute perm in third dimension and pass the 3 independent components
kyy = mean(kz_loc);     
Perm = [K(1), kyy, K(4)];

% Compare outflux if requested
if U.outflux == 1
    % Pressure drop along the axial (Y) dimension, with sealed lateral
    % boundaries. Remember that MRST uses SI units (perm m^2).
    
    % Coarse grid
    upscaled = [1 10];              % cellnum in each direction
    disp('________________________________________________________________')
    disp(' Pressure drop with sealed boundaries: ')
    disp(['Coarse grid cells along axial direction is ' ...
          num2str(upscaled(2))]);
    G_ups = computeGeometry(cartGrid(upscaled, [faultThick faultDisp]));
    p = partitionCartGrid(G.cartDims, upscaled);
    CG = generateCoarseGrid(G, p);
    
    % Coarse grid permeability is the result of the upscaling
    crock.perm = repmat([K(1) K(4)], CG.cells.num, 1);
    
    % Define fluid
    fluid = initSingleFluid('mu', 1*centi*poise, ...
                            'rho', 1014*kilogram/meter^3);       % water
    
    % Define representative initial conditions
    watCol = 50*meter;
    topy   = 1000*meter;
    g      = 9.861;
    p_r    = (watCol + topy + min(G.cells.centroids(:,2)))*fluid.rhoWS*g;
    [z_0, z_max] = deal(min(G.cells.centroids(:,2)), ...
        max(G.cells.centroids(:,2)));
    equil  = ode23(@(z,p) g .* fluid.rhoWS, [z_0, z_max], p_r);
    p0     = flipud(reshape(deval(equil, G.cells.centroids(:,2)), [], 1));
    state0 = initResSol(G, p0, 1); clear equil
    
    p_s = (watCol + topy + min(G_ups.cells.centroids(:,2)))*fluid.rhoWS*g;
    [z_0u, z_maxu] = deal(min(G_ups.cells.centroids(:,2)), ...
                          max(G_ups.cells.centroids(:,2)));
    if z_0u ~= z_maxu
        equil  = ode23(@(z,p) g .* fluid.rhoWS, [z_0u, z_maxu], p_s);
        p0u = flipud(reshape(deval(equil, G_ups.cells.centroids(:,2)), [], 1));
        clear equil
    elseif G_ups.cells.num == 1
        p0u = (watCol + topy + G_ups.cells.centroids(:,2))*fluid.rhoWS*g;
    end
    state0c = initResSol(G_ups, p0u, 1);
    
    % Fine-scale problem
    bc        = pside([], G, 'North', 0.75*min(p0));
    faces     = bc.face;
    bc        = pside(bc, G, 'South',  1.25*max(p0));
    if strcmp(U.method, 'mpfa')
        x   = incompMPFA(state0, G, hTmp, fluid, 'bc', bc);
        str  = 'Sum outflux on fine scale (MPFA): ';
    else
        hT   = computeTrans(G, rock);
        x    = incompTPFA(state0, G, hT, fluid, 'bc', bc);
        str  = 'Sum outflux on fine scale (TPFA): ';
    end
    
    % Coarse-scale problem
    bc_ups    = pside([], G_ups, 'North', 0.75*min(p0));
    faces_ups = bc_ups.face;
    bc_ups    = pside(bc_ups, G_ups, 'South', 1.25*max(p0));
    if strcmp(U.method, 'mpfa')
        hTmp_ups = computeMultiPointTrans(G_ups, crock, 'invertBlocks', inB);
        x_ups = incompMPFA(state0, G, hTmp_ups, fluid, 'bc', bc);
        str_ups  = 'Sum outflux on coarse scale (MPFA): ';
    else
        T_ups = computeTrans(G_ups, crock);
        x_ups = incompTPFA(state0c, G_ups, T_ups, fluid, 'bc', bc_ups);
        str_ups  = 'Sum outflux on coarse scale (TPFA): ';
    end
    
    
    % Results
    flux = sum(x.flux(faces));
    flux_ups = sum(x_ups.flux(faces_ups));
    disp([str, num2str(flux)]);
    disp([str_ups, num2str(flux_ups)]);
    disp('______________________________________________________________');
    axialDropVals = [flux, flux_ups];
    
end

if U.ARcheck == 1 % check if same perm with smaller aspect ratio. Slow.
    assert(G.griddim == 2)               % way too many elements in 3D
    if ~U.useAcceleration
        error('U.useAcceleration must be set to true. Too slow otherwise')
    end
    ar = G.xzFaceDim(2)/G.xzFaceDim(1);
    ar2 = 5;
    if ar > ar2
        disp(['Checking if cell aspect ratio is too high for', ...
              ' flow-based upscaling']);
        Lz = G.cartDims(2)*G.xzFaceDim(2);
        nels = [G.cartDims(1) round(Lz/(ar2*G.xzFaceDim(1)))];
        G2 = computeGeometry(cartGrid([nels(1), nels(2)], ...
                                      [faultThick, faultDisp]));
        disp(['Finer grid with cell aspect ratio = ' num2str(ar2) ...
              ' created.'])
        disp(['Finer grid has a total of ' num2str(G2.cells.num) ...
              ' cells. Now computing cell distances...'])
        
        % assign perm based on cell centroid distance to closest cell
        distx = pdist2(G2.cells.centroids(:,1), G.cells.centroids(:,1));
        distz = pdist2(G2.cells.centroids(:,end), G.cells.centroids(:,end));
        dist = sqrt(distx.^2 + distz.^2);
        disp(['Distances computed. Now calculating MPFA perm with ' ...
              'finer z mesh...'])
        
        % Compute 2nd perm
        [~, mapG2toG] = min(dist,[],2);
        rock2.perm = rock.perm(mapG2toG, :);
        hTmp2 = computeMultiPointTrans(G2, rock2, 'invertBlocks', 'mex');
        psolver = @(state0, G, fluid, bc) incompMPFA(state0, G, hTmp2, ...
                                                     fluid, 'bc', bc);
        fluid = initSingleFluid('mu' , 1, 'rho', 1);
        K2 = diag(myupscalePermeabilityFixed(G2, Dp{1}, psolver, fluid, L));
        disp(['Finer z mesh perm [x, z], in m^2, is [' ...
              num2str(K2(1)) ' ' num2str(K2(end)) ']'])
        disp(['Ratio Perm(finerGrid) / Perm(originalGrid) [x, z] is [' ...
              num2str(K2(1)/K(1)) ' ' num2str(K2(end)/K(end)) ']'])
        
        % Plot
        %figure(15); plotCellData(G2, log10(rock2.perm(:,end)./(milli*darcy)), ...
        %                         'facea', 1, 'edgea', 0);
        %colormap(copper); title('Finer Z Grid vertical Perm')
        %xlim([0 f.T]); ylim([0 f.D]); c = colorbar;
    else
        disp(['aspect ratio < ' num2str(ar2) '. Check not needed'])
    end
end

end