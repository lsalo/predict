function [poroG, permG, kz_loc] = setGridPoroPerm(obj, G, FS)
%
%
%

% Rename
M = obj.MatMap;
alpha = obj.Alpha;
unitPoro = obj.MatProps.Poro;
unitPerm = obj.MatProps.Perm;
permAnisoRatio = obj.MatProps.PermAnisoRatio;


% 1) Map diagonals to Grid indexing: Grid indexing starts at bottom left,
%    columns (x) move faster. Matrix starts counting at top left and rows
%    (y) move faster. This is for 2D grid. For 3D grid, the x-z plane
%    starts at top left, so we just need to transpose.

% figure(9); spy(sparse(M.vals))       % check that smears are correctly placed in M.
% 2D grid (xy) mapping
isSmear = reshape(transpose(flipud(M.vals)), G.cells.num, 1);
whichUnit = reshape(transpose(flipud(M.units)), G.cells.num, 1);

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
T = [cosd(alpha) sind(alpha); -sind(alpha) cosd(alpha)];% Transform. mat

for n=1:numel(M.unit)
    cellIds = [all([isSmear, whichUnit == M.unit(n)], 2), ...
               all([~isSmear, whichUnit == M.unit(n)], 2)];
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
        % Get MatProps values/fcns
        phi = unitPoro(M.unit(n), :);     % M.unit can differ from M.unitIn
        kx_loc(cellIds) = unitPerm{M.unit(n)}(cellNum);

        % Porosity
        if abs(diff(phi)) > 0.01
            poroG(cellIds) = phi(1) + rand(cellNum, 1) .* abs(diff(phi));
        else
            poroG(cellIds) = phi(1);
        end
        
        % Permeability
        %kx_loc = unitPerm{n}(cellNum);
        kz_loc(cellIds) = kx_loc(cellIds) * permAnisoRatio(M.unit(n));
        sz = [2, 2];
        nels = sum(sz);
        S = zeros(sz(1), sz(2), numel(kx_loc(cellIds)));
        idDiag = [1, 4];               % [kxx, kzz] of smear (diag perm)
        idS = [idDiag; idDiag + ...
            cumsum(repmat(repelem(nels, 2), cellNum-1, 1))]';
        S(reshape(idS, numel(idS), 1)) = [kx_loc(cellIds), kz_loc(cellIds)]';
        kmat = transformTensor(S, T, 'posDef');
        id = [1,2,4];                  % [kxx, kxz, kzz] of fault
        idk = [id; id + cumsum(repmat(repelem(nels, 3), cellNum-1, 1))]';
        idk = reshape(idk, numel(idk), 1);
        permG(cellIds, :) = reshape(kmat(idk), 3, cellNum)';
    
    else                                % clay domain, disc. smear
        % Get sand/clay values
        phi_c = unitPoro(M.unit(n), :);     
        kx_loc_c = unitPerm{M.unit(n)}(cellNum(1));
        kz_loc_c = kx_loc_c * permAnisoRatio(M.unit(n));
        
        if any(~M.isclayIn)
            if any(~M.isclay)
                [~, sid] = min(abs(M.unit(~M.isclay) - M.unit(n)));
                sandIds = M.unit(~M.isclay);
            else
                [~, sid] = min(abs(M.unitIn(~M.isclayIn) - M.unit(n)));
                sandIds = M.unitIn(~M.isclayIn);
            end
            closestSandId = sandIds(sid);
            phi_s = unitPoro(closestSandId, :);
            kx_loc_s = unitPerm{closestSandId}(cellNum(2));
            kz_loc_s = kx_loc_s * permAnisoRatio(closestSandId);
            
        else
            disp('_______________________________________________________')
            warning('No sand in stratigraphy.')
            disp(['Properties for sand material between discontinuous ' ...
                  'smears calculated assuming the following:'])
            disp('     Clay content = 0.15')
            disp('     Faulting Depth = Footwall Faulting Depth')
            disp('     Burial Depth = Footwall Burial Depth');
            disp('     Permeability anisotropy ratio = 1');
            disp('_______________________________________________________')
            
            % create dummy section, we will just use the FW;
            thicks = [25, 25];
            vcls   = [0.15, 0.5];
            dips   = [0, 0];
            zf     = FS.FW.DepthFaulting;
            zmax   = FS.FW.DepthBurial;
            
            % FW and HW
            fw = Stratigraphy(thicks(1), vcls(1), dips(1), ...
                              'DepthFaulting', zf(1), ...
                              'DepthBurial', zmax(1));
            hw = Stratigraphy(thicks(2), vcls(2), dips(2), 'IsHW', 1, ...
                              'NumLayersFW', fw.NumLayers);
            
            % Strati in Faulted Section
            mysect = FaultedSection(fw, hw);
            
            % Compute sand poro and perm
            sandPermUpperBound = 1000;          % [mD]
            [perm, poro] = getPermeability(mysect, sandPermUpperBound);
            
            % Assign;
            phi_s = poro(1, :);
            kx_loc_s = perm{1}(cellNum(2));
            kz_loc_s = kx_loc_s;
        end
        
        
        % Porosity
        if abs(diff(phi_c)) > 0.01
            poroG(cellIds(:, 1)) = phi_c(1) + rand(cellNum(1), 1) .* abs(diff(phi_c));
        else
            poroG(cellIds(:, 1)) = phi_c(1);
        end
        if abs(diff(phi_s)) > 0.01
            poroG(cellIds(:, 2)) = phi_s(1) + rand(cellNum(2), 1) .* abs(diff(phi_s));
        else
            poroG(cellIds(:, 2)) = phi_s(1);
        end
        
        % Permeability
        cellNum = sum(cellNum);
        %kx_loc = zeros(10,1); kz_loc = zeros(10, 1);
        kx_loc(cellIds(:, 1)) = kx_loc_c; kx_loc(cellIds(:, 2)) = kx_loc_s;
        kz_loc(cellIds(:, 1)) = kz_loc_c; kz_loc(cellIds(:, 2)) = kz_loc_s;
        %kx_loc(kx_loc == 0) = []; kz_loc(kz_loc == 0) = [];
        cellIds = sort([find(cellIds(:, 1)); find(cellIds(:, 2))]);
        
        sz = [2, 2];
        nels = sum(sz);
        S = zeros(sz(1), sz(2), numel(kx_loc(cellIds)));
        idDiag = [1, 4];               % [kxx, kzz] of smear (diag perm)
        idS = [idDiag; idDiag + ...
                       cumsum(repmat(repelem(nels, 2), cellNum-1, 1))]';
        S(reshape(idS, numel(idS), 1)) = [kx_loc(cellIds), kz_loc(cellIds)]';
        kmat = transformTensor(S, T, 'posDef');
        id = [1,2,4];                  % [kxx, kxz, kzz] of fault
        idk = [id; id + cumsum(repmat(repelem(nels, 3), cellNum-1, 1))]';
        idk = reshape(idk, numel(idk), 1);
        permG(cellIds, :) = reshape(kmat(idk), 3, cellNum)';
        
    end
end

end