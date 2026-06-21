function M = faultMaterialMap(G, FS, smear, varargin)
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
    % -------------------------OPTIONAL INPUTS--------------------------------
    % 'SmearOverlapRule' = how clay-smear overlaps are handled. One of:
    %     'random' (default) : legacy behavior. In each overlapping diagonal
    %                          group a single parent clay source is selected
    %                          (uniform PMF), then smears are realized by
    %                          object-based simulation (placeSmearObjects).
    %     'cell_union_psmear': each clay source is kept as a separate,
    %                          untrimmed domain here; the stochastic
    %                          realization, cell-wise union of active source
    %                          masks, and parent labeling are deferred to
    %                          placeSmearObjects. In this mode the function
    %                          returns per-source band descriptors and sets
    %                          M.cellUnion = true (see below); it does NOT
    %                          resolve overlaps or place sand diagonals.
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
    
    % Method to handle clay smear overlaps and sand units
    sand_method = 'closest';
    overlap_method = 'random';                    % 'random' or 'pffr' or 'vcl'

    % Optional input: clay-smear overlap rule. 'random' keeps the legacy
    % single-parent path; 'cell_union_psmear' switches to the deferred
    % per-source realization + union handled in placeSmearObjects.
    opt.SmearOverlapRule = 'random';
    opt = merge_options_relaxed(opt, varargin{:});
    smearOverlapRule = normalizeSmearOverlapRule(opt.SmearOverlapRule);
    
    % Initial values to Matrix structure
    nx            = G.cartDims(1);
    nz            = G.cartDims(end);
    lz            = G.cellDim(end);
    M.nDiagTot    = 2*nx - 1;                       % total number of diagonals
    M.vals        = false(nx);                      % actual matrix of 0s and 1s
    M.units       = zeros(size(M.vals));            % Unit domain of each cell  (parent Id)
    M.unit        = FS.ParentId;                    % Unit domain of each group (parent Id)
    M.isclay      = [FS.FW.IsClay, FS.HW.IsClay];   % total units and clay or not  
    M.unitIn      = M.unit;                         % For reference (unchanged)
    M.isclayIn    = M.isclay;                       % "
    idc           = find(M.isclay);


    %% 1. Place clay smears
    % 1.1: Compute physical and grid geometry for each layer
    M = computeLayerGeometry(FS, M, nx);
    
    % 1.2 Calculate number of diagonals with potential smear and assign
    %     initial domains
    M = assignInitialDiags(FS, M, smear, nx, idc);

    % 1.2b cell_union_psmear: hand off untrimmed per-source band descriptors
    %      to the placement stage and return. We deliberately skip overlap
    %      resolution (1.4), the full-fault check (1.3/1.5), and sand-diagonal
    %      assignment (2.x): in this mode each clay source is realized
    %      independently and overlaps are resolved by cell-wise union (with
    %      nearest-source / nearest-sand labeling) inside placeSmearObjects.
    if strcmp(smearOverlapRule, 'cell_union_psmear')
        M = buildCellUnionDescriptors(M, FS, smear, idc, nz, lz);
        return
    end

    % 1.3: Check if smear covers the whole fault area
    smearThickAsFault = checkFullFaultSmear(M, nx, idc);
    
    % 1.4: Assign parent units to diagonal groups accounting for overlaps
    M = assignUnitsToDiags(M, nx, idc, overlap_method);

    % 1.5 Check smearThickAsFault
    if smearThickAsFault
        assert(sum(M.nDiag) == M.nDiagTot)
    end
    
    % 1.6 Post-processing of clay domains for clarity
    M = finalizeClayDomains(FS, M, smear, idc);


    %% 2. Diagonals with sand
    % 2.1 Assign sand domains between clay smear domains
    idc = find(M.isclay); % needs to be updated
    M = assignSandDomains(M, nx, idc, sand_method);  

    % 2.2 Add field divLayerDiag for layer with lower and upper diags
    M = assignDivLayerDiag(M);

    % 2.3 Check that we have the correct number of diagonals for each layer.
    assert(all(M.nDiag > 0))
    assert(sum(M.nDiag) == M.nDiagTot)


    %% 3. Populate mapping matrix with all potential 1s and sure 0s
    M = populateMappingMatrices(FS, M, nz, lz, sand_method);


    %% Helper functions
    function M = computeLayerGeometry(FS, M, nx)
        %
        % Position of stratigraphic layers with respect to diags in M
        % 

        TapFW = FS.Tap(FS.FW.Id);
        TapHW = FS.Tap(FS.HW.Id);
        faultDisp = sum(TapFW);

        M.layerCenter = [cumsum(TapFW) - TapFW/2, cumsum(TapHW) - TapHW/2];
        M.layerTop = [cumsum(TapFW), cumsum(TapHW)];
        M.layerBot = [0, M.layerTop(FS.FW.Id(1:end-1)), 0, M.layerTop(FS.HW.Id(1:end-1))];


        c = nx / faultDisp;
        M.layerDiagCenter = round(M.layerCenter * c) - nx;
        M.layerDiagCenter(FS.HW.Id) = M.layerDiagCenter(FS.HW.Id) + nx;

        % Just for plotting/etc
        M.layerDiagTop = round(M.layerTop * c) - nx;
        M.layerDiagTop(FS.HW.Id) = M.layerDiagTop(FS.HW.Id) + nx;
        M.layerDiagTop(M.layerDiagTop > 0) = M.layerDiagTop(M.layerDiagTop > 0) - 1;

        M.layerDiagBot = round(M.layerBot * c) - nx;
        M.layerDiagBot(FS.HW.Id) = M.layerDiagBot(FS.HW.Id) + nx + 1;
        M.layerDiagBot(M.layerDiagBot == 0) = 1;
        M.layerDiagBot(1) = -(nx - 1);

        M.nDiagLayer = (M.layerDiagTop - M.layerDiagBot) + 1;
    end


    function M = assignInitialDiags(FS, M, smear, nx, idc)
        %
        %
        %
        assert(M.nDiagTot == 2*nx-1);
        
        layerLs = zeros(1, max(FS.ParentId));
        layerLs(M.isclay) = smear.ThickInFault;
        M.nDiag = round((layerLs ./ smear.DomainLength) * M.nDiagTot);
        
        % Assign initial DiagTop and DiagBot
        M.DiagTop = M.layerDiagCenter + fix((M.nDiag-1)/2);
        M.DiagBot = M.layerDiagCenter - round((M.nDiag-1)/2);
        M.DiagTop(~M.isclay) = 0;
        M.DiagBot(~M.isclay) = 0;


        % 1.3. Adjust nDiag to within -G.cartDims and +G.cartDims
        diagEnd = nx - 1;
        M.DiagBot(M.DiagBot < -diagEnd) = -diagEnd;
        M.DiagTop(M.DiagTop > diagEnd) = diagEnd;
        M.nDiag(idc) = abs(M.DiagTop(idc) - M.DiagBot(idc)) + 1;
    end


    function smearThickAsFault = checkFullFaultSmear(M, nx, idc)
        
        smearThickAsFault = 0;
        if sum(M.nDiag) >= sum(M.nDiagTot)  % smear may occupy the full fault area
            diag_ids = -nx + 1 : nx - 1;
            clay_diag = cell2mat(arrayfun(@(x, y) x:y, M.DiagBot(idc), M.DiagTop(idc), 'uniformoutput', false));
            clay_diag = unique(clay_diag);
            diag_ids(ismember(diag_ids, clay_diag)) = [];
            if isempty(diag_ids)
                smearThickAsFault = 1;
            end
        end
    end


    function M = assignUnitsToDiags(M, nx, idc, overlap_method)
        %
        % Unit selection in each diagonal group, accounting for potential 
        % overlaps. Methods: see selectUnitIndex
        %    
        % Note that, in case of overlaps, this may lead to the same unit 
        % appearing more than once and non-consecutively, and units 
        % may no longer be centered with respect to source layer.
        %
        
        % Find all potential units present in each diagonal
        Omap = buildOverlapMap(M, nx, idc);
        
        % Find diagonal groups with same potentially present units
        diagsGroup = groupDiagsByOverlap(Omap);
        
        % For each diagonal group, select a parent unit and assign
        [M, unitGroup] = selectAndAssignUnits(M, nx, diagsGroup, Omap, overlap_method);
        
        % Finalize: update nDiag, check for non-appearing units, assign DiagBot/DiagTop
        [M, ~] = finalizeDiagUnitAssignments(M, idc, unitGroup);
        
        % -----------------
        % Internal helpers
        % -----------------
        function Omap = buildOverlapMap(M, nx, idc)
            nDiag = nx - 1;
            nc = sum(M.isclay);
            Omap = zeros(nDiag, nc);
            for k = 1:nc
                id0 = [M.DiagBot(idc(k)) + nx, M.DiagTop(idc(k)) + nx];
                Omap(id0(1):id0(2), k) = idc(k);
            end
        end
        
        function diagsGroup = groupDiagsByOverlap(Omap)
            idChange = find(any(diff(Omap) ~= 0, 2));
            diagsGroup = [[1; idChange + 1], [idChange; size(Omap, 1)]];
            diagsGroup(~any(Omap(diagsGroup(:, 1), :), 2), :) = []; % remove diags w/o clay
        end
        
        function [M, unitGroup] = selectAndAssignUnits(M, nx, diagsGroup, Omap, overlap_method)
            %
            %   Select final unit randomly, and assign unit to each diagonal group 
            %
            unitGroup = zeros(size(diagsGroup, 1), 1);
            DiagBot = nan(1, numel(M.DiagBot));
            DiagTop = nan(1, numel(M.DiagBot));
            len = numel(M.DiagBot);
            
            for k = 1:size(diagsGroup, 1)
                repeatedUnitNonConsec = false;
                id_source = unique(Omap(diagsGroup(k, 1), :));
                id_source(id_source == 0) = [];
                
                % Choose using the selection method 
                selectedIdx = selectUnitIndex(id_source, overlap_method);
                unitGroup(k) = id_source(selectedIdx);
                
                if isnan(DiagBot(unitGroup(k)))     % new unit
                    DiagBot(unitGroup(k)) = diagsGroup(k, 1) - nx;
                    DiagTop(unitGroup(k)) = diagsGroup(k, 2) - nx;
                else
                    % Unit was already assigned to a group of diags, so we need to
                    % check what limits (DiagBot and/or DiagTop) we need to extend.
                    itBot = diagsGroup(k, 1) - nx;
                    itTop = diagsGroup(k, 2) - nx;
                    idsThisUnit = unitGroup == unitGroup(k);
                    assert(itBot > DiagBot(unitGroup(k)));
                    if ~idsThisUnit(k-1)            % non consecutive unit appearance
                        repeatedUnitNonConsec = true;
                        len = len + 1;
                        DiagBot(len) = itBot;
                        M.unit(len) = unitGroup(k);
                        M.isclay(len) = true;
                    end
                    lastIdThisUnit = find(M.unit == unitGroup(k), 1, 'last');
                    if idsThisUnit(k-1) && itTop > DiagTop(lastIdThisUnit)
                        DiagTop(lastIdThisUnit) = itTop;
                    elseif ~idsThisUnit(k-1)
                        assert(repeatedUnitNonConsec == true)
                        DiagTop(len) = itTop;
                    end
                end
            end
            M.DiagBot = DiagBot;
            M.DiagTop = DiagTop;
        end
        
        function idx = selectUnitIndex(id_source, method)
            % 'random': all units overlapping have same probability of
            %           being selected.
            switch lower(method)
                case 'random'
                    idx = randi(numel(id_source), 1);   
                    
                otherwise
                    error('Selection method not implemented');
            end
        end
        
        function [M, clayUnitsAssigned] = finalizeDiagUnitAssignments(M, idc, unitGroup)
            % Update nDiag, remove non-appearing units, fill DiagBot/DiagTop, etc
            M.nDiag = (M.DiagTop - M.DiagBot) + 1;
            M.nDiag(~M.isclay) = 0; % 0 for sand intervals
            clayUnitsAssigned = unique(unitGroup);
            idNotPresent = ~ismember(idc, clayUnitsAssigned);
            M.nDiag(idc(idNotPresent)) = 0; % 0 for non-appearing clays
            assert(sum(M.nDiag) <= M.nDiagTot);
            M.DiagBot(isnan(M.DiagBot)) = 0;
            M.DiagTop(isnan(M.DiagTop)) = 0;
        end  
    end


    function M = finalizeClayDomains(FS, M, smear, idc)
        % Track which smear domains were removed or repeated; finalize sorting.
        %
        % M.Psmear: Vector of Psmear values for each domain (repeats for repeated unit)
        % M.idSmearInRemoved: Indices of removed (zero-diagonal) domains
        % M.unit, .isclay, .nDiag, etc: Reduced to only real clay domains
        % All fields sorted by DiagBot for consistent order
        
        % Step 1: Assign Psmear for all clay units in input, repeat for repeated units,
        %         and remove zero-Psmear (i.e. removed) domains.
        M.Psmear = zeros(1, max(M.unitIn));
        M.Psmear(M.isclayIn) = smear.Psmear;
        M.Psmear = M.Psmear(M.unit);             % repeat Psmear for repeated units
        M.Psmear(M.Psmear == 0) = [];            % remove zeros (non-appearing or sand)
        
        if any(M.nDiag(idc) == 0)
            M.idSmearInRemoved = find(M.nDiag(idc) == 0);   % zero-diag units
            M.Psmear(M.idSmearInRemoved) = [];              % remove non-appearing smear domains from Psmear
            idSelectedUnit = idc(M.idSmearInRemoved);
            % Clean all arrays for these units
            M.unit(idSelectedUnit)    = [];
            M.isclay(idSelectedUnit)  = [];
            M.nDiag(idSelectedUnit)   = [];
            M.DiagBot(idSelectedUnit) = [];
            M.DiagTop(idSelectedUnit) = [];
        end
        
        % Step 2: Sort by DiagBot and keep only actual clays (nDiag > 0)
        idToSort = find(M.nDiag > 0);
        [~, pos] = sort(M.DiagBot(idToSort));
        
        M.DiagBot = M.DiagBot(idToSort(pos));
        M.DiagTop = M.DiagTop(idToSort(pos));
        M.nDiag   = M.nDiag(idToSort(pos));
        M.unit    = M.unit(idToSort(pos));
        M.isclay  = M.isclay(idToSort(pos));
        M.Psmear  = M.Psmear(pos);  % Only the clays that survived
        M.isPFFRIn = all([FS.Vcl >= 0.3; FS.Vcl < 0.4]);
        M.TapIn   = FS.Tap;
        M.VclIn   = FS.Vcl;
        M.disp    = sum(FS.Tap(1:FS.FW.Id(end)));
    end


    function M = assignSandDomains(M, nx, idc, method_sand)
        %
        % Place sand domains in diagonals not yet filled by clay, using the chosen method.
        %
        diagIds = -nx + 1 : nx - 1;
        clayDiag = unique(cell2mat(arrayfun(@(x, y) x:y, M.DiagBot(idc), M.DiagTop(idc), 'uniformoutput', false)));
        sandDiag = diagIds(~ismember(diagIds, clayDiag));  % These need to be sand
        
        if isempty(sandDiag)
            return;
        end
        
        % 1. Find contiguous sand diagonal groups
        sandGroups = findContiguousGroups(sandDiag);
        
        % 2. Decide sand parent unit for each group (using method_sand)
        sandParents = chooseSandParents(M, sandGroups, method_sand);
        
        % 3. Add sand domains to M arrays, merge/clay as needed
        M = insertSandDomains(M, sandGroups, sandParents);
        
        % -----------------
        % Internal helpers
        % -----------------
        function sandGroups = findContiguousGroups(sandDiag)
            %
            % Returns a 2 x nGroups array where each column is [start_diag; end_diag]
            %
            breakIdx = [find(diff(sandDiag) > 1), numel(sandDiag)];
            startIdx = [1, breakIdx(1:end-1)+1];
            sandGroups = [sandDiag(startIdx); sandDiag(breakIdx)];
        end
        
        function sandParents = chooseSandParents(M, sandGroups, method_sand)
            nGroups = size(sandGroups, 2);
            sandParents = zeros(1, nGroups);
            sand_unit_ids = M.unitIn(~M.isclayIn);
            sand_centers_grid = M.layerDiagCenter(~M.isclayIn);
            
            for k = 1:nGroups
                center_diag = mean(sandGroups(:, k));
                sandParents(k) = getSandUnit(M, center_diag, sand_centers_grid, sand_unit_ids, method_sand);
            end
        end
        
        function M = insertSandDomains(M, sandGroups, sandParents)
            % Add new sand domains, re-sort all as needed
            
            % Concatenate with existing
            M.DiagBot = [M.DiagBot, sandGroups(1,:)];
            M.DiagTop = [M.DiagTop, sandGroups(2,:)];
            M.unit    = [M.unit,   sandParents];
            M.isclay  = [M.isclay, false(1, numel(sandParents))];
            
            % Sort all by DiagBot
            [M.DiagBot, sortIdx] = sort(M.DiagBot);
            M.DiagTop = M.DiagTop(sortIdx);
            M.unit    = M.unit(sortIdx);
            M.isclay  = M.isclay(sortIdx);
            
            % Recompute nDiag
            M.nDiag = abs(M.DiagBot - M.DiagTop) + 1;
        end

    end


    function M = assignDivLayerDiag(M)
        %
        % Adds or updates the field divLayerDiag in M.
        %
        if any(all([M.DiagBot < 0; M.DiagTop > 0]))
            idLay = all([M.DiagBot < 0; M.DiagTop > 0]);
            M.divLayerDiag = [abs(M.DiagBot(idLay)) + 1, M.DiagTop(idLay)];
            
        elseif any(M.DiagBot(M.DiagTop == 0) == 0)
            idLay = M.DiagBot == 0;
            assert(M.nDiag(idLay) == 1);
            M.divLayerDiag = [1, 0];
            
        elseif any(M.DiagBot == 0)
            idLay = M.DiagBot == 0;
            M.divLayerDiag = [0, M.nDiag(idLay)];
            
        elseif any(M.DiagTop == 0)
            idLay = M.DiagTop == 0;
            M.divLayerDiag = [M.nDiag(idLay), 0];
            
        else
            M.divLayerDiag = []; 
        end
        
    end

    function M = populateMappingMatrices(FS, M, nz, lz, sand_method)
        %
        % Populates mapping matrices M.vals (clay/sand) and M.units (parent material),
        % and fills M.unitInClayGaps and .windowTop/Bot for each domain.
        %
        nUnits = numel(M.nDiag);
        M.unitInClayGaps = nan(1, nUnits);
        M.windowTop = nan(1, nUnits);
        M.windowBot = nan(1, nUnits);
        
        
        for n = 1:nUnits
            assert(M.DiagTop(n) >= M.DiagBot(n))
            
            if M.isclay(n)
                % 1. If this is a clay unit, find sand unit for gaps
                if any(~M.isclayIn)
                    M.unitInClayGaps(n) = selectSandUnitForGap(M, n, sand_method);
                end
                
                % 2. Smear window for this clay domain
                [idTop, idBot] = getSmearWindow(M, FS, nz, lz, M.unit(n));
                M.windowTop(n) = idTop;
                M.windowBot(n) = idBot;
                
                % 3. Populate diag values (M.vals) for this clay domain
                diagVals = false(nz, M.nDiag(n));
                diagVals(idBot:idTop, :) = true;
                diagVals = flipud(diagVals);  % for spdiags
                
            else
                diagVals = false(nz, M.nDiag(n));
            end
            diagUnits = M.unit(n)*ones(nz, M.nDiag(n));
            
            % 4. Write values and parent units to mapping matrices
            M.vals = full(spdiags(diagVals, -M.DiagTop(n):-M.DiagBot(n), M.vals));
            M.units = full(spdiags(diagUnits, -M.DiagTop(n):-M.DiagBot(n), M.units));
        end        
        M.units = transpose(M.units);  % M.vals is transposed in placeSmearObjects (as before)
        
        % -----------------
        % Internal helpers
        % -----------------
        function sandId = selectSandUnitForGap(M, k, method)
            
            % Selects sand unit (parent id) for gap in smear domain n, given desired method.           
            sand_unit_ids = M.unitIn(~M.isclayIn);
            sand_centers = M.layerDiagCenter(~M.isclayIn);
            clay_center = mean([M.DiagBot(k); M.DiagTop(k)]);
            sandId = getSandUnit(M, clay_center, sand_centers, sand_unit_ids, method);
        end
        
    end


    function sandId = getSandUnit(M, target, sand_values, sand_unit_ids, method)
        %
        % target: scalar, typically central diagonal or elevation
        % sand_values: array (same units as target), represents sand domain "centers"
        % sand_unit_ids: array of parent IDs, corresponds to sand_values
        % method: string, e.g. 'closest', ...
        %
        
        % Calculate absolute distance, consider juxtaposition proximity as well
        n = round((M.nDiagTot + 1) / 2);            
        if target < 0
            target = n + target;
        end
        ids = sand_values < 0;
        sand_values(ids) = n + sand_values(ids);
        dist = abs(target - sand_values);       
        
        switch lower(method)
            case 'closest'
                [~, sid] = min(dist);
                sandId = sand_unit_ids(sid);
                
            otherwise
                error('Unknown method for sand gap filling: %s', method);
        end
        
    end

    function [idTop, idBot] = getSmearWindow(M, FS, nz, lz, unitId)
        % Top and bottom cell indices (rows of "diagVals" used in spdiags)
        % of the elevation window in which a clay source can smear: a FW
        % source only at z <= its layer top, a HW source only at z >= its
        % layer bottom. Shared by populateMappingMatrices (legacy path, called
        % with M.unit(n)) and buildCellUnionDescriptors (cell-union path,
        % called per clay source). Takes the parent unit id directly.
        zTot = nz * lz;     % total height

        if ismember(unitId, FS.FW.Id)     % Smear source is FW
            idTop = round((M.layerTop(unitId) / zTot) * nz);
            idTop = min(idTop, nz);
            idBot = 1;
        else                              % Smear source is HW
            idTop = nz;
            idBot = round((M.layerBot(unitId) / zTot) * nz);
            if idBot == 0, idBot = 1; end
        end
    end

    function M = buildCellUnionDescriptors(M, FS, smear, idc, nz, lz)
        %
        % cell_union_psmear handoff. Reduce M to the list of clay SOURCES,
        % each kept untrimmed (no overlap resolution), described by its
        % diagonal band (DiagBot/DiagTop/nDiag), elevation window
        % (windowTop/windowBot) and Psmear. The stochastic realization,
        % cell-wise union of active source masks, and parent-unit labeling
        % are deferred to placeSmearObjects (which detects M.cellUnion).
        %
        % Full-length input/geometry fields are intentionally retained:
        %   - M.unitIn, M.isclayIn : needed to identify sand parent units and
        %                            to map per-source smear lengths.
        %   - M.layerDiagCenter, M.layerTop, M.layerBot, ... : needed for the
        %                            nearest-source and nearest-sand labeling
        %                            in the placement stage (indexed by the
        %                            actual parent unit id, not by source #).
        %
        % NOTE on indexing: after this function the reduced descriptor arrays
        % (M.unit, M.DiagBot/Top, M.nDiag, M.windowTop/Bot, M.Psmear) are
        % aligned to the clay-source order (idc), i.e. indexed by source #.
        % The retained layer-geometry arrays remain indexed by parent unit id.

        nSrc = numel(idc);
        assert(numel(smear.Psmear) == nSrc, ...
               'smear.Psmear length must equal the number of clay sources.')

        % Per-source elevation window (rows of the diagonal reachable by the
        % smear). Mirrors getSmearWindow in the legacy populate step: a FW
        % source can only smear at z <= its layer top; a HW source only at
        % z >= its layer bottom.
        windowTop = nan(1, nSrc);
        windowBot = nan(1, nSrc);
        for s = 1:nSrc
            [windowTop(s), windowBot(s)] = getSmearWindow(M, FS, nz, lz, idc(s));
        end

        % Reduce to clay-source descriptors (aligned to idc)
        M.unit      = idc;                       % clay source unit ids
        M.isclay    = true(1, nSrc);
        M.DiagBot   = M.DiagBot(idc);
        M.DiagTop   = M.DiagTop(idc);
        M.nDiag     = M.nDiag(idc);
        M.windowTop = windowTop;
        M.windowBot = windowBot;
        M.Psmear    = reshape(smear.Psmear, 1, []);   % aligned to idc

        % Mode flag (consumed by placeSmearObjects and Fault2D.placeMaterials)
        M.cellUnion        = true;
        M.smearOverlapRule = 'cell_union_psmear';
    end


    function rule = normalizeSmearOverlapRule(ruleIn)
        %
        % Canonical name for the clay-smear overlap rule. Only 'random'
        % (legacy) and 'cell_union_psmear' are supported here; the 'geologic'
        % rule from the prototype is intentionally omitted.
        %
        if isstring(ruleIn)
            ruleIn = char(ruleIn);
        end
        assert(ischar(ruleIn), ...
               'SmearOverlapRule must be a character vector or string scalar.')
        rule = lower(strtrim(ruleIn));
        switch rule
            case {'random', 'legacy', 'uniform_random'}
                rule = 'random';
            case {'cell_union_psmear', 'cell-union-psmear', ...
                  'union_psmear', 'psmear_cell_union'}
                rule = 'cell_union_psmear';
            otherwise
                error(['Unknown SmearOverlapRule "%s". Use "random" or ' ...
                       '"cell_union_psmear".'], ruleIn)
        end
    end


end