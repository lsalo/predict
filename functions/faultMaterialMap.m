function M = faultMaterialMap(G, FS, smear)
%
% -----------------------------SUMMARY------------------------------------
% This function takes as inputs the Grid structure (G), faulted section (FS)
% and smear dimensions, and creates a structure, M, which contains
% the information regarding the clay smears and sand distribution in
% the modeled domain (the fault). The mapping matrix itself, stored in
% M.vals (0 = sand, 1 = clay smear) and M.units (parent unit in each 
% domain), provides a direct map to the materials in the simulation grid.
%
% The matrix M.vals is initialized with all potential 1s (all cells in the 
% fault that would contain smear if all smears were continuous), and 0s 
% (all cells in the fault that will surely contain sand).
%
% Each smear (which occupies a given number of diagonals in M.vals as well
% as the simulation grid) is initialized around the middle of the
% corresponding unit in the FW or HW. The number of diagonals are decided
% based on the corresponding thickness. In case of clay smear overlaps,
% a single parent is selected according to a uniform PMF (random).
%
% **Reminder**: Directly superposed sand layers (consecutive in FW or HW)
% should be collapsed in FW and HW variables. That is, FW, for example,
% should be passed as thicknesses [20 60 20] for vcl [A B A], where A is 
% above smearing threshold and B is below, and not [20 20 20 20 20] with 
% vcl [A B B B A]. Conversely, consecutive clay layers above the smearing 
% threshold can be passed, e.g. thickness [20 20 20 40 20] for vcl 
% [A A A B A].
% 
% Note: MRST grid indexing (G) starts at bottom left, columns (x) move 
% faster. Standard MATLAB matrices start counting at top left, and rows (y) 
% move faster.
%
% ------------------------------INPUTS------------------------------------
% G   = MRST Grid structure
% FS  = FaultedSection object with valid fields.
% smear = Smear object with valid fields.
%
% ------------------------------OUTPUT------------------------------------
% M = matrix structure. Contains the following fields:
%   nDiagTot        = number of total diagonals in grid of size n*n, i.e. 
%                     2*n - 1.
%   vals            = matrix of n*n that contains the mapping of 1s (smear)
%                     and 0s (sand) to the grid. Note that the domains with
%                     parent materials with vcl > smearing threshold are
%                     just initialized with all 1s here, and the final
%                     configuration is defined in another function, based 
%                     on the chosen geostatistical method.
%   units           = matrix of n*n that indicates which parent unit is
%                     present in each cell.
%   unit            = 1xk array where k is the number of different domains
%                     in the fault. Note that, in some cases, this number
%                     can be different from the input number of units given
%                     overlap handling.
%   isclay          = logical array indicating whether the material in each
%                     domain is sand-based (0) or clay smear (1).
%   nDiag           = number of diagonals that each domain has in the
%                     grid.
%   layerCenter     = center height of each layer in the FW and HW, 
%                     relative to the bottom of the modeled fault portion. 
%   layerTop        = top height of each layer in the FW and HW.
%   layerBot        = base height of each layer in the FW and HW.
%   layerDiagCenter = Fault grid diag number corresponding to the center of 
%                     each layer in the FW and HW. Follows the same
%                     convention as MATLAB's spdiags, i.e. 0 is the main
%                     diagonal, lower diagonals are negative and upper are
%                     positive. 
%   layerDiagTop    = Grid diagonal number corresponding to the top of 
%                     each layer in the FW and HW.
%   layerDiagBot    = Grid diagonal number corresponding to the base of 
%                     each layer in the FW and HW.
%   nDiagLayer      = Number of diagonals between layerDiagTop and
%                     layerDiagBot, i.e. the number of diagonals in contact
%                     with each stratigraphic layer in the FW and HW.
%   DiagTop         = Grid diagonal number corresponding to the top of 
%                     each material domain in the fault.
%   DiagBot         = Grid diagonal number corresponding to the bottom of 
%                     each material domain in the fault.
%   clayDiagBot     = Grid diagonal number corresponding to the bottom of 
%                     each clay source domain in the fault.
%   unitIn          = array with input units (parent materials), i.e. 1 
%                     (bottom of fw) to j (top of hw). Same as FS.ParentId.
%   isclayIn        = logical array indicating whether each input unit is
%                     sand (0) or clay (1).
%   divLayerDiag    = Number of diagonals in the lower and upper triangle
%                     of the grid, when there is a layer whose 
%                     corresponding subdomain in the fault crosses the 
%                     main diagonal. Otherwise [0, 0];
%   idSmearInRemoved = index of M.unitIn(M.isclayIn) corresponding to
%                      smearing sources not appearing (i.e. not 
%                      contributing material) to the fault, as a result of 
%                      overlap handling.
%   Psmear          = fraction of domains with smearing sources as parent
%                     material that is actually occupied by smear 
%                     (continuous smear = 1). Result of calculation based
%                     on SSFc performed by the Smear class.
%   unitInClayGaps  = sand unit to be placed in each smear domain in case
%                     of discontinuous smear. The selected sand is the
%                     closest to the domain location in the fault. Only
%                     defined if at least one sand unit is present in the
%                     stratigraphy.
%   windowBot       = In the fault, a smear contributed by a given parent 
%                     can only appear at an elevation (z) <= top of parent
%                     in FW, and >= bottom of parent in HW (window in which
%                     it is sheared). windowBot corresponds to the bottom
%                     id of each smear window in a diagonal with m 
%                     elements, where m is the grid size.
%   windowTop       = Top id of each smear window in a diagonal with m
%                     elements, where m is the grid size.
%   P               = Psmear (1st row) and result of object simulation in
%                     placeSmearObjects.m (performed later).
%   
%__________________________________________________________________________

% Initial values to Matrix structure
if G.griddim == 3
    assert(~isfield(G, 'cartDims'))
    assert(mod(sqrt(G.layerSize), 1) == 0)
    G.cartDims = [G.numLayers sqrt(G.layerSize) sqrt(G.layerSize)];
    id_dim = 2;
elseif G.griddim == 2
    id_dim = 1;
end
M.nDiagTot    = 2*G.cartDims(id_dim) - 1;            % total number of diagonals
M.vals        = false(G.cartDims(id_dim));           % actual matrix of 0s and 1s
M.units       = zeros(size(M.vals));            % Unit domain of each cell  (parent Id)
M.unit        = FS.ParentId;                    % Unit domain of each group (parent Id)
M.isclay      = [FS.FW.IsClay, FS.HW.IsClay];   % total units and clay or not  
M.unitIn      = M.unit;                         % For reference (unchanged)
M.isclayIn    = M.isclay;                       % "
idc           = find(M.isclay);
        

%% 1. Place clay smears

% 1.1 Calculate number of diagonals with potential smear
layerLs = zeros(1, max(FS.ParentId));
layerLs(M.isclay) = smear.ThickInFault;
M.nDiag = round((layerLs./smear.DomainLength)*M.nDiagTot);


% 1.2 Position of stratigraphic layers with respect to diags in M
faultDisp = sum(FS.Tap(FS.FW.Id));
M.layerCenter = [cumsum(FS.Tap(FS.FW.Id))-FS.Tap(FS.FW.Id)/2 ...
                 cumsum(FS.Tap(FS.HW.Id))-FS.Tap(FS.HW.Id)/2];
M.layerTop    = [cumsum(FS.Tap(FS.FW.Id)) cumsum(FS.Tap(FS.HW.Id))];
M.layerBot    = [0 M.layerTop(FS.FW.Id(1:end-1)) ...
                 0 M.layerTop(FS.HW.Id(1:end-1))];
M.layerDiagCenter = round(M.layerCenter.*(G.cartDims(id_dim)/faultDisp)) - ...
                          G.cartDims(id_dim);
M.layerDiagCenter(FS.HW.Id) = M.layerDiagCenter(FS.HW.Id) + G.cartDims(id_dim);

% Assign initial DiagTop and DiagBot
M.DiagTop = M.layerDiagCenter + fix((M.nDiag-1)/2);
M.DiagBot = M.layerDiagCenter - round((M.nDiag-1)/2);
M.DiagTop(~M.isclay) = 0;
M.DiagBot(~M.isclay) = 0;

% Just for convenience in plotting, etc.
M.layerDiagTop = round(M.layerTop.*(G.cartDims(id_dim)/faultDisp)) - ...
                       G.cartDims(id_dim);
M.layerDiagTop(FS.HW.Id) = M.layerDiagTop(FS.HW.Id) + G.cartDims(id_dim);
M.layerDiagTop(M.layerDiagTop>0) = M.layerDiagTop(M.layerDiagTop>0)-1;
M.layerDiagBot = round(M.layerBot.*(G.cartDims(id_dim)/faultDisp)) - ...
                       G.cartDims(id_dim);
M.layerDiagBot(FS.HW.Id) = M.layerDiagBot(FS.HW.Id) + G.cartDims(id_dim)+1;
M.layerDiagBot(M.layerDiagBot==0) = 1;
M.layerDiagBot(1) = -(G.cartDims(id_dim)-1);
M.nDiagLayer = (M.layerDiagTop - M.layerDiagBot)+1;


% 1.3. Adjust nDiag to within -G.cartDims and +G.cartDims
diagEnd = G.cartDims(id_dim) - 1;
M.DiagBot(M.DiagBot < -diagEnd) = -diagEnd;
M.DiagTop(M.DiagTop > diagEnd) = diagEnd;
M.nDiag(idc) = abs(M.DiagTop(idc) - M.DiagBot(idc)) + 1;


% 1.4 Determine whether the total smear area occupies the full fault (do
%     not add overlapping areas)
smearThickAsFault = 0;
if sum(M.nDiag) >= sum(M.nDiagTot) % smear may occupy the full fault area
    diagIds = -G.cartDims(id_dim)+1:G.cartDims(id_dim)-1;
    clayDiag = cell2mat(arrayfun(@(x,y) x:y, M.DiagBot(idc), ...
                                 M.DiagTop(idc),'uniformoutput',false));
    clayDiag = unique(clayDiag);
    diagIds(ismember(diagIds, clayDiag)) = [];
    
    if isempty(diagIds)            % smear occupies the full fault area
        smearThickAsFault = 1;
    end
end


% 1.5 Unit selection in each diagonal group, accounting for potential 
%     overlaps: Randomly select unit in overlapping areas. Note that, in 
%     case of overlaps, this may lead to the same unit appearing more than 
%     once and non-consecutively. Moreover, units may no longer be centered 
%     with respect to source layer in HW or FW (only in case of overlaps).

% Find all units potentially present in each diagonal
nDiag = sum(G.cartDims(1:2)) - 1;
nc = sum(M.isclay);
Omap = zeros(nDiag, nc);
for n = 1:nc
   id0 = [M.DiagBot(idc(n)) + G.cartDims(id_dim), ...
          M.DiagTop(idc(n)) + G.cartDims(id_dim)]; 
   Omap(id0(1):id0(2), n) = idc(n);
end
%Omap(~any(Omap, 2), : ) = [];

% Find diagonal groups with same potentially present units
idChange = find(any(diff(Omap) ~= 0, 2));
diagsGroup = [[1; idChange + 1], [idChange; size(Omap, 1)]];
diagsGroup(~any(Omap(diagsGroup(:, 1), :), 2), :) = []; % remove diags w/o clay


% Select final unit randomly, and assign unit to each diagonal group 
unitGroup = zeros(size(diagsGroup, 1), 1);
DiagBot = nan(1, numel(M.DiagBot));
DiagTop = nan(1, numel(M.DiagBot));
len = numel(M.DiagBot);
for n = 1:size(diagsGroup, 1)
    repeatedUnitNonConsec = false;
    vals = unique(Omap(diagsGroup(n, 1), :));  
    vals(vals == 0) = [];
    idSelectedUnit = randi(numel(vals), 1);
    unitGroup(n) = vals(idSelectedUnit);
    if isnan(DiagBot(unitGroup(n)))                     % new unit
        DiagBot(unitGroup(n)) = diagsGroup(n, 1) - G.cartDims(id_dim);
        DiagTop(unitGroup(n)) = diagsGroup(n, 2) - G.cartDims(id_dim);
        
    else  % Unit was already assigned to a group of diags, so we need to
          % check what limits (DiagBot and/or DiagTop) we need to extend.
        itBot = diagsGroup(n, 1) - G.cartDims(id_dim);
        itTop = diagsGroup(n, 2) - G.cartDims(id_dim);
        idsThisUnit = unitGroup == unitGroup(n);
        assert(itBot > DiagBot(unitGroup(n))) 
            %error('Check what is going on.')
            %DiagBot(unitGroup(n)) = itBot;     
        if idsThisUnit(n-1) == false    % non consecutive unit appearance
            repeatedUnitNonConsec = true;
            len = len + 1;
            DiagBot(len) = itBot; 
            M.unit(len) = unitGroup(n);
            M.isclay(len) = true;
        end
        lastIdThisUnit = find(M.unit == unitGroup(n), 1, 'last');
        if idsThisUnit(n-1) == true && itTop > DiagTop(lastIdThisUnit)
            DiagTop(lastIdThisUnit) = itTop;
        elseif idsThisUnit(n-1) == false
            assert(repeatedUnitNonConsec == true)
            DiagTop(len) = itTop;
        end
    end
end
M.nDiag = (DiagTop - DiagBot) + 1;
M.nDiag(~M.isclay) = 0;                     % 0 diags for sand intervals
clayUnitsAssigned = unique(unitGroup);      
idNotPresent = ~ismember(idc, clayUnitsAssigned);
M.nDiag(idc(idNotPresent)) = 0;             % 0 diags for non-appearing c
assert(sum(M.nDiag) <= M.nDiagTot);
M.DiagBot = DiagBot;    M.DiagBot(isnan(M.DiagBot)) = 0;
M.DiagTop = DiagTop;    M.DiagTop(isnan(M.DiagTop)) = 0;

% 1.6  If any M.nDiag is 1 we neglect it (For now, we don't)
if any(M.nDiag == 1) 
    nid = find(M.nDiag == 1);
    assert(all(M.DiagTop(nid) == M.DiagBot(nid)))
    %M.DiagTop(nid) = 0;
    %M.DiagBot(nid) = 0;
    %M.nDiag(nid) = 0;
    %for n = 1:numel(nid)
    %    
    %    if smearThickAsFault % add to other unit randomly
    %        
    %    end
    %end
end


% 1.7 Check smearThickAsFault
if smearThickAsFault
    try
        assert(sum(M.nDiag) == M.nDiagTot)
    catch
        error('smear not as thick as fault!')
    end
end


% 1.8 Track which smear domains were removed (if 0 diag) or repeated (for
%     Psmear)
M.Psmear = zeros(1, max(M.unitIn));
M.Psmear(M.isclayIn) = smear.Psmear;
M.Psmear = M.Psmear(M.unit);       % repeat Psmear for repeated units
M.Psmear(M.Psmear == 0) = [];      % remove 0s (sands)
if any(M.nDiag(idc) == 0)
    M.idSmearInRemoved = find(M.nDiag(idc) == 0);  
    M.Psmear(M.idSmearInRemoved) = [];      % remove non-appearing smear domains
    idSelectedUnit = idc(M.idSmearInRemoved);
    M.unit(idSelectedUnit) = [];
    M.isclay(idSelectedUnit) = [];
    M.nDiag(idSelectedUnit) = [];
    M.DiagBot(idSelectedUnit) = []; M.DiagTop(idSelectedUnit) = [];
end


% 1.9 Sort based on appearing diagonal group, and finalize this stage with
%     only clays in the following fields.
idToSort = find(M.nDiag > 0);
[~, pos] = sort(M.DiagBot(idToSort));
M.DiagBot = M.DiagBot(idToSort(pos));
M.DiagTop = M.DiagTop(idToSort(pos));
M.nDiag   = M.nDiag(idToSort(pos));
M.unit    = M.unit(idToSort(pos));
M.isclay  = M.isclay(idToSort(pos));
M.Psmear = M.Psmear(pos);



%% 2. Diagonals with sand
idc     = find(M.isclay);
diagIds = -G.cartDims(id_dim)+1:G.cartDims(id_dim)-1;
clayDiag = cell2mat(arrayfun(@(x,y) x:y, M.DiagBot(idc), ...
                             M.DiagTop(idc),'uniformoutput',false));
clayDiag = unique(clayDiag);
diagIds(ismember(diagIds, clayDiag)) = [];
%flag = 0;
if numel(diagIds) > 0
        if numel(diagIds) > 1
            diffsf = diff(diagIds)>1;
            diffsi = [false diffsf(1:end-1)];
            if diffsf(end) == 1                 
                diffsi(end+1) = 1;
            end
            sandIds = [diagIds(1) diagIds(diffsi); diagIds(diffsf) diagIds(end)];
        else
            sandIds = [diagIds; diagIds];
        end
        M.clayDiagBot = M.DiagBot;
        M.DiagBot = sort([M.DiagBot sandIds(1,:)]);
        M.DiagTop = sort([M.DiagTop sandIds(2,:)]);
        
        idSandAll = M.unitIn(~M.isclayIn);
        sandCenterDiag = mean(sandIds)';
        stratiSandCenterDiag = M.layerDiagCenter(~M.isclayIn)';
        stratiCenterRep = repmat(stratiSandCenterDiag, ...
                                  [1 length(sandCenterDiag)]);
        [~, idClosestSand] = min(abs(stratiCenterRep-sandCenterDiag'), [], 1);
        sandParent = idSandAll(idClosestSand);
        [~, ids] = intersect(M.DiagBot, sandIds(1,:));
%         idSandAll = M.unitIn(~M.isclayIn);
%         sandParent = zeros(1, numel(ids));
%         for n=1:numel(ids)
%             [~, idClosestSand] = min(abs(ids(n) - idSandAll));
%             sandParent(n) = idSandAll(idClosestSand);
%         end

        nunits = sum([sum(M.isclay), numel(ids)]);
        unitAll = zeros(1, nunits);
        unitAll(ids) = sandParent;
        idclay = unitAll == 0;
        unitAll(idclay) = M.unit(M.isclay);
        M.unit = unitAll;
        M.isclay = false(1, numel(M.unit));
        M.isclay(idclay) = true;
%         if any(M.DiagTop(M.DiagBot == 0) == 0)      % remove clay smear domain with 0 diags.
%             %error('To code.')
%             idTopZ = find(M.DiagTop == 0);
%             idBotZ = find(M.DiagBot == 0);
%             idTopBotZ = idTopZ(idTopZ == idBotZ);
%             M.DiagTop(idTopBotZ) = [];
%             M.DiagBot(idTopBotZ) = [];
%             M.unit(idTopBotZ) = [];
%             M.isclay(idTopBotZ) = [];
%         end
%         if any(M.DiagBot == 0) && any(M.DiagTop == 0) %&& ~all(M.DiagTop(M.DiagBot==0)==0)
%             error('To code.')
%             %M.DiagTop(M.DiagTop == 0) = [];
%             %M.DiagBot(M.DiagBot == 0) = [];
%         end
        M.nDiag = abs(M.DiagBot - M.DiagTop)+1;
end



%% 3. Add field divLayerDiag for layer with lower and upper diags
if any(all([M.DiagBot<0; M.DiagTop>0]))
    idLay = all([M.DiagBot<0; M.DiagTop>0]);
    %M.DiagBot(idLay) = M.DiagBot(idLay)+1;
    M.divLayerDiag = [abs(M.DiagBot(idLay))+1 M.DiagTop(idLay)];
else
    if M.DiagBot(M.DiagTop == 0) == 0
        idLay = M.DiagBot == 0;
        assert(M.nDiag(idLay) == 1)
        M.divLayerDiag = [1 0];
    elseif any(M.DiagBot == 0)
        idLay = M.DiagBot == 0;
        M.divLayerDiag = [0 M.nDiag(idLay)];
    elseif any(M.DiagTop == 0)
        idLay = M.DiagTop == 0;
        M.divLayerDiag = [M.nDiag(idLay) 0];
    end
end



%% Check that we have the correct number of diagonals for each layer.
assert(all(M.nDiag > 0))
try
    assert(sum(M.nDiag) == M.nDiagTot)
catch
    error('Material mapping was unsuccessful.') 
end



%% 4. Populate mapping matrix with all potential 1s and sure 0s
M.unitInClayGaps = nan(1, numel(M.unit));
M.windowTop = nan(1, numel(M.unit));
M.windowBot = nan(1, numel(M.unit));
for n = 1:numel(M.nDiag)    
    assert(M.DiagTop(n) >= M.DiagBot(n))
    if M.isclay(n) == 1
        % Select sand unit to be placed in smear gap, if any (used later)
        if any(~M.isclayIn)
            cCenter = mean([M.DiagBot(n); M.DiagTop(n)]);
            sCenters = M.layerDiagCenter(~M.isclayIn);
            [~, sid] = min(abs(cCenter - sCenters));
            sandIds = M.unitIn(~M.isclayIn);
            closestSandId = sandIds(sid);
            M.unitInClayGaps(n) = closestSandId;
        end
        
        % 4.1 Select elevation window in which each clay smear can be found 
        % (never at higher z than top of source in FW, and never at lower z 
        % than bottom of source in HW).
        % Top and bottom ids of smear window
        if ismember(M.unit(n), FS.FW.Id)
            idTop = round(( M.layerTop(M.unit(n))/...
                            (G.cartDims(end)*G.cellDim(end)) ) * G.cartDims(end));
            idTop(idTop > G.cartDims(end)) = G.cartDims(end);
            idBot = 1;
        else
            idTop = G.cartDims(end);
            idBot = round(( M.layerBot(M.unit(n))/...
                            (G.cartDims(end)*G.cellDim(end)) ) * G.cartDims(end)); 
            idBot(idBot == 0) = 1;
        end
        M.windowTop(n) = idTop;
        M.windowBot(n) = idBot;     
        diagVals = false(G.cartDims(end), M.nDiag(n));        % unitInClayGaps (sand)
        diagVals(idBot:idTop, :) = true;                    % clay
        diagVals = flipud(diagVals);                        % for spdiags
        
        % 4.2 Populate mapping matrix
        M.vals = full(spdiags(diagVals, -M.DiagTop(n):-M.DiagBot(n), ...
                              M.vals));
        
    else
        % Populate mapping matrix directly
        M.vals = full(spdiags(false(G.cartDims(end), M.nDiag(n)), ...
                              -M.DiagTop(n):-M.DiagBot(n), M.vals));
    end
    M.units = full(spdiags(M.unit(n)*ones(G.cartDims(end), M.nDiag(n)), ...
                   -M.DiagTop(n):-M.DiagBot(n), M.units));
end
M.units = transpose(M.units);
% M.vals also needs to be transposed. We do it later in placeSmearObjects.

end