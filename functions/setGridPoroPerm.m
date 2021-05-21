function [poroG, permG, kz_loc, M] = setGridPoroPerm(obj, G, FS)
%
%
%

% Rename
M = obj.MatMap;
alpha = obj.Alpha;
unitPoro = obj.MatProps.poro;
unitPerm = obj.MatProps.perm;
permAnisoRatio = obj.MatProps.permAnisoRatio;
M.unitInClayGaps = nan(1, numel(M.unit));


% 1) Map diagonals to Grid indexing: Grid indexing starts at bottom left,
%    columns (x) move faster. Matrix starts counting at top left and rows
%    (y) move faster. This is for 2D grid. For 3D grid, the x-z plane
%    starts at top left, so we just need to transpose.

% figure(9); spy(sparse(M.vals))       % check that smears are correctly placed in M.
% 2D grid (xy) mapping
isSmear = reshape(transpose(flipud(M.vals)), G.cells.num, 1);
%whichUnit = reshape(transpose(flipud(M.units)), G.cells.num, 1);

% 2) Assign permeability and porosity. We have end-member porosity, which 
%    is assigned to the corresponding unit (clay or sand). For the
%    permeability, we generate N random samples according to the
%    distribution for each domain, and number of cells for each domain in
%    the fault. This is then directly assumed to be the permeability
%    across the clay/sand smears (principal directions). 
%    We then compute the permeability along the sand and clay smears based 
%    on anisotropy ratios.
%    This diagonal permeability tensor needs to be rotated to the fault 
%    local axes, since clay/sand smears are not parallel to the fault.

% Poro
% Assign based on endMemberPoro to corresponding cells. Random sampling for
% clay.     
poroG = nan(G.cells.num, 1);                            % sand poro default
permG = nan(G.cells.num, 3);                            % [kxx, kxz, kzz]
kx_loc = nan(G.cells.num, 1);
kz_loc = nan(G.cells.num, 1);
T = [cosd(alpha) sind(alpha); ...
     -sind(alpha) cosd(alpha)];                         % Transform. mat
fn = 0.01;                                              % porosity var factor
fk = 1.75;                                              % perm var factor

for n=1:numel(M.unit)
    % To get cellIds, we use idsUnitBlock instead of whichUnit == M.unit(n)
    % since same unit can appear more than once, in separate blocks.
    idsUnitBlock = ~full(spdiags(zeros(G.cartDims(1),M.nDiag(n)), ...
                                 M.DiagBot(n):M.DiagTop(n), M.units));
    idsUnitBlock = reshape(transpose(flipud(idsUnitBlock)), G.cells.num, 1);
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
    
    if numel(cellNum) == 1        
        % Give a factor of about 3 max variation in permeability (already
        % correlated to fault thickness through shear strain), and 0.03 to
        % porosity.
        unitPoroCorrRange = [unitPoro(M.unit(n))-fn unitPoro(M.unit(n))+fn];
        assert(all(unitPoroCorrRange > 0))
        unitPermCorrRangeLog = log10([unitPerm(M.unit(n))/fk unitPerm(M.unit(n))*fk]);
        % Variable value each cell
        poroG(cellIds) = unitPoroCorrRange(1) + ...
                         rand(cellNum, 1)*diff(unitPoroCorrRange);    
        kx_loc(cellIds) = 10.^(unitPermCorrRangeLog(1) + ...
                               rand(cellNum, 1)*diff(unitPermCorrRangeLog));
        
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
        % Get sand/clay values        
        % Clay porosity
        unitPoroCorrRange = [unitPoro(M.unit(n))-fn unitPoro(M.unit(n))+fn];
        assert(all(unitPoroCorrRange > 0))
        poroG(cellIds(:, 1)) = unitPoroCorrRange(1) + ...
                               rand(cellNum(1), 1)*diff(unitPoroCorrRange);
        % Clay permeability
        unitPermCorrRangeLog = log10([unitPerm(M.unit(n))/fk unitPerm(M.unit(n))*fk]);
        kx_loc(cellIds(:, 1)) = 10.^(unitPermCorrRangeLog(1) + ...
                                rand(cellNum(1), 1)*diff(unitPermCorrRangeLog));
        kz_loc(cellIds(:, 1)) = kx_loc(cellIds(:, 1))*permAnisoRatio(M.unit(n));
        
        if any(~M.isclayIn)     % if sand in stratigraphy
            cCenter = mean([M.DiagBot(n); M.DiagTop(n)]);
            sCenters = M.layerDiagCenter(~M.isclayIn);
            [~, sid] = min(abs(cCenter - sCenters));
            sandIds = M.unitIn(~M.isclayIn);
%             if any(~M.isclay)
%                 [~, sid] = min(abs(M.unit(~M.isclay) - M.unit(n)));
%                 sandIds = M.unit(~M.isclay);
%             else
%                 [~, sid] = min(abs(M.unitIn(~M.isclayIn) - M.unit(n)));
%                 sandIds = M.unitIn(~M.isclayIn);
%             end
            closestSandId = sandIds(sid);
            phi_s = unitPoro(closestSandId);
            permx_s = unitPerm(closestSandId);
            kprime_s = permAnisoRatio(closestSandId);
            M.unitInClayGaps(n) = closestSandId;
            
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
            phi_s = getPorosity(0.15, 0.4, zf, zmax, 'zmax', isuc);
            permx_s = getPermeability(0.15, 0.4, zf, zmax, cap);
            kprime_s = 1;
        end
        
        % Sand porosity
        poroSandRange = [phi_s-fn phi_s+fn];
        assert(all(poroSandRange > 0))
        poroG(cellIds(:, 2)) = poroSandRange(1) + ...
                               rand(cellNum(2), 1)*diff(poroSandRange); 
        
        % Sand Permeability
        permxSandRangeLog = log10([permx_s/fk permx_s*fk]);
        kx_loc(cellIds(:, 2)) = 10.^(permxSandRangeLog(1) + ...
                                rand(cellNum(2), 1)*diff(permxSandRangeLog));
        kz_loc(cellIds(:, 2)) = kx_loc(cellIds(:, 2)) * kprime_s;
        
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