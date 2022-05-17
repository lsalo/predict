function [poroG, permG, kz_loc, vcl, units] = setGridPoroPerm(fault, G, FS)
% Assign permeability and porosity to each grid cell, based on the material
% mapping outcome.
%
% INPUT:
%   fault: an instance of Fault with MatMap property.
%   G: MRST grid structure.
%   FS: an instance of FaultedSection.
%
%
% DESCRIPTION:
%   This function assigns porosity and permeability to each grid cell.
%   The process is done one domain at a time. First, we determine whether 
%   1 unit (continuous clay smear with window filling the whole domain, or 
%   just sand) or 2 units (clay smear + sand) are present in the domain.
%   Then, sand or clay porosity/permeability is assigned based on the
%   parent units, and a small range of variation is introduced to account 
%   for their higher variation frequency (e.g. Grant, 2019). Permeability 
%   is allowed to vary a $\log_{10}$ value of $\pm 0.2$, while porosity is 
%   allowed $\pm 0.015$. The permeability along each
%   material is given based on the computed permeability anisotropy ratios.
%   Finally, these permeabilities (along and across each material, not
%   parallel to the fault zone) are transformed to the fault local axes, so
%   that a full 2D perm tensor is obtained for each cell (permG = kxx, kxz, 
%   kzz). These transformed values are the cell permeabilities used for
%   flow-based upscaling, while Kyy (upscaled) is the average of the
%   k-along prior to tensor rotation (kz_loc).
%   
%
% NOTE:
%   Sand-material properties in domains with smear but not fully covered by
%   smear are assigned based on the closest sand unit in the stratigraphy
%   (closest parent). If no sand is present in the stratigraphy, but the
%   clay smears are discontinuous, sand-based fault material with generic
%   properties is assigned. This case (no sand in stratigraphy) should be
%   better modeled future code versions, but it requires some changes for 
%   handling discontinuous smears.
%
% OUTPUT:
%  poroG: nx1 array with cell porosity, where n is the number of cells in G
%  permG: nx3 array with kxx, kxz, kzz (aligned with fault axes, ready for
%         upscaling).
%  kz_loc: nx1 array with permeability along each fault material, for
%          calculating upscaled permeability along strike (Kyy).
%  vcl:   nx1 array with clay content of each fault cell, based on vcl of
%         parent contributing material to each cell.
%--------------------------------------------------------------

% Rename
M = fault.MatMap;
alpha = fault.Alpha;
unitPoro = fault.MatProps.poro;
unitPerm = fault.MatProps.perm;
permAnisoRatio = fault.MatProps.permAnisoRatio;

% Griddim
id_dim = 1;
if G.griddim == 3
    %ncell = G.layerSize;    % extruded
    %id_dim = 2;
    ncell = G.cartDims(1)*G.cartDims(3); 
elseif G.griddim == 2
    ncell = G.cells.num;
end


% 1) Map diagonals to Grid indexing: Grid indexing starts at bottom left,
%    columns (x) move faster. Matrix starts counting at top left and rows
%    (y) move faster. This is for 2D grid. For 3D grid, the x-z plane
%    starts at top left, so we just need to transpose.

% figure(9); spy(sparse(M.vals))       % check that smears are correctly placed in M.
% 2D grid (xy) mapping
isSmear = reshape(transpose(flipud(M.vals)), ncell, 1);
%whichUnit = reshape(transpose(flipud(M.units)), ncell, 1);

% 2) Assign permeability and porosity based on parent unit(s) in each 
%    domain. This is assumed to be the permeability across the clay/sand 
%    smears (principal directions). 
%    We then compute the permeability along the sand and clay smears based 
%    on anisotropy ratios.
%    This diagonal permeability tensor needs to be rotated to the fault 
%    local axes, since clay/sand smears are not parallel to the fault.

% Poro
% Assign based on endMemberPoro to corresponding cells. Random sampling for
% clay.     
poroG = nan(ncell, 1);                            % sand poro default
permG = nan(ncell, 3);                            % [kxx, kxz, kzz]
vcl   = nan(ncell, 1);
kx_loc = nan(ncell, 1);
kz_loc = nan(ncell, 1);
units = nan(ncell, 1);
T = [cosd(alpha) sind(alpha); ...
     -sind(alpha) cosd(alpha)];                         % Transform. mat
fn = 0.01;                                              % porosity var factor
fk = 0.2;                                               % log perm var factor

for n=1:numel(M.unit)
    % To get cellIds, we use idsUnitBlock instead of whichUnit == M.unit(n)
    % since same unit can appear more than once, in separate blocks.
    idsUnitBlock = ~full(spdiags(zeros(G.cartDims(id_dim),M.nDiag(n)), ...
                                 M.DiagBot(n):M.DiagTop(n), M.units));
    idsUnitBlock = reshape(transpose(flipud(idsUnitBlock)), ncell, 1);
    cellIds = [all([isSmear, idsUnitBlock], 2), ...
               all([~isSmear, idsUnitBlock], 2)];
    cellNum = [sum(cellIds(:, 1)), sum(cellIds(:, 2))];
    
    if M.isclay(n) && cellNum(2) == 0
        cellIds = cellIds(:, 1);
        cellNum = cellNum(1);
    elseif ~M.isclay(n)
        assert(cellNum(1) == 0);
        cellIds = cellIds(:, 2);
        cellNum = cellNum(2);
    end
    
    if numel(cellNum) == 1              % clay or sand domain, continuous
        % Give a factor of about 3 max variation in permeability (already
        % correlated to fault thickness through shear strain), and 0.03 to
        % porosity.
        units(cellIds) = M.unit(n);
        unitPoroCorrRange = [unitPoro(M.unit(n))-fn unitPoro(M.unit(n))+fn];
        assert(all(unitPoroCorrRange > 0))
        unitPermCorrRangeLog = log10(unitPerm(M.unit(n))) + [-fk, fk];
        % Variable value each cell
        poroG(cellIds) = unitPoroCorrRange(1) + ...
                         rand(cellNum, 1)*diff(unitPoroCorrRange);    
        kx_loc(cellIds) = 10.^(unitPermCorrRangeLog(1) + ...
                               rand(cellNum, 1)*diff(unitPermCorrRangeLog));
        vcl(cellIds) = FS.Vcl(M.unit(n));
        % Permeability
        %kx_loc = unitPerm{n}(cellNum);
        kz_loc(cellIds) = kx_loc(cellIds) * permAnisoRatio(M.unit(n));
        sz = [2, 2];
        nels = sum(sz);
        S = zeros(sz(1), sz(2), numel(kx_loc(cellIds)));
        idDiag = [1, 4];               % [kxx, kzz] of smear (diag perm)
        idS = [idDiag; idDiag + ...
               cumsum(repmat(repelem(nels, 2), cellNum-1, 1), 1)]';
        S(reshape(idS, numel(idS), 1)) = [kx_loc(cellIds), kz_loc(cellIds)]';
        kmat = transformTensor(S, T, 'posDef');
        id = [1,2,4];                  % [kxx, kxz, kzz] of fault
        idk = [id; id + cumsum(repmat(repelem(nels, 3), cellNum-1, 1), 1)]';
        idk = reshape(idk, numel(idk), 1);
        permG(cellIds, :) = reshape(kmat(idk), 3, cellNum)';
    
    else                                % clay domain, disc. smear
         units(cellIds(:, 1)) = M.unit(n);
         assert(~isnan(M.unitInClayGaps(n)))
         units(cellIds(:, 2)) = M.unitInClayGaps(n);
         
        % Get sand/clay values        
        % Clay porosity
        unitPoroCorrRange = [unitPoro(M.unit(n))-fn unitPoro(M.unit(n))+fn];
        assert(all(unitPoroCorrRange > 0))
        poroG(cellIds(:, 1)) = unitPoroCorrRange(1) + ...
                               rand(cellNum(1), 1)*diff(unitPoroCorrRange);
        % Clay permeability
        unitPermCorrRangeLog = log10(unitPerm(M.unit(n))) + [-fk, fk];
        kx_loc(cellIds(:, 1)) = 10.^(unitPermCorrRangeLog(1) + ...
                                rand(cellNum(1), 1)*diff(unitPermCorrRangeLog));
        kz_loc(cellIds(:, 1)) = kx_loc(cellIds(:, 1))*permAnisoRatio(M.unit(n));
        
        % Clay Vcl
        vcl(cellIds(:,1)) = FS.Vcl(M.unit(n));
        
        if any(~M.isclayIn)     % if sand in stratigraphy
            phi_s = unitPoro(M.unitInClayGaps(n));
            permx_s = unitPerm(M.unitInClayGaps(n));
            kprime_s = permAnisoRatio(M.unitInClayGaps(n));
            vcl_s = FS.Vcl(M.unitInClayGaps(n));
            
        else                    % no sand in stratigraphy
            disp('_______________________________________________________')
            warning('No sand in stratigraphy.')
            disp(['Properties for sand material between discontinuous ' ...
                  'smears calculated assuming the following:'])
            disp('     Clay content = 0.15')
            disp('     Faulting Depth = Footwall Faulting Depth')
            disp('     Burial Depth = Footwall Burial Depth');
            disp('     Permeability anisotropy ratio = 1');
            disp('_______________________________________________________')
            
            % FS data
            zf     = FS.FW.DepthFaulting;
            zmax   = FS.FW.DepthBurial;
            isuc   = FS.IsUndercompacted;
            
            % Compute sand poro and perm
            cap = 1000;          % [mD]
            vcl_s = 0.15;
            phi_s = getPorosity(vcl_s, FS.IsClayVcl, zf, zmax, 'zmax', isuc);
            permx_s = getPermeability(vcl_s, FS.IsClayVcl, zf, zmax, cap);
            kprime_s = 1;
        end
        
        % Sand porosity
        poroSandRange = [phi_s-fn phi_s+fn];
        assert(all(poroSandRange > 0))
        poroG(cellIds(:, 2)) = poroSandRange(1) + ...
                               rand(cellNum(2), 1)*diff(poroSandRange); 
        
        % Sand Permeability
        permxSandRangeLog = log10(permx_s) + [-fk, fk];
        kx_loc(cellIds(:, 2)) = 10.^(permxSandRangeLog(1) + ...
                                rand(cellNum(2), 1)*diff(permxSandRangeLog));
        kz_loc(cellIds(:, 2)) = kx_loc(cellIds(:, 2)) * kprime_s;
        
        % Sand Vcl
        vcl(cellIds(:, 2)) = vcl_s;
        
        % Transform
        cellIds = sort([find(cellIds(:, 1)); find(cellIds(:, 2))]);        
        sz = [2, 2];
        nels = sum(sz);
        S = zeros(sz(1), sz(2), numel(kx_loc(cellIds)));
        idDiag = [1, 4];               % [kxx, kzz] of smear (diag perm)
        id = [1,2,4];                  % [kxx, kxz, kzz] of fault
        cellNum = sum(cellNum);
        if cellNum == 2
            idS = [idDiag; idDiag + repelem(nels, 2)]';
            idk = [id; id + repelem(nels, 3)]';
        else
            idS = [idDiag; idDiag + ...
                          cumsum(repmat(repelem(nels, 2), cellNum-1, 1))]';
            idk = [id; id + cumsum(repmat(repelem(nels, 3), cellNum-1, 1))]';
        end
        S(reshape(idS, numel(idS), 1)) = [kx_loc(cellIds), kz_loc(cellIds)]';
        kmat = transformTensor(S, T, 'posDef');
        idk = reshape(idk, numel(idk), 1);
        permG(cellIds, :) = reshape(kmat(idk), 3, cellNum)';  
    end
end

end