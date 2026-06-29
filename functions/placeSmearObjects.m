function M = placeSmearObjects(M, smear, FS, G, tolerance, verbose)
%
% -----------------------------SUMMARY------------------------------------
% This function randomly places smear objects using object-based simulation, 
% so that each realization provides a random case consistent with the
% inputs. The final goal is to assign correct number of 0s and 1s in M.vals 
% such that probabilities in M.Psmear are matched to the indicated 
% tolerance. These probabilities are matched in 2D, ie accounting for the 
% total number of cells in each domain. Smears, however, are treated as 
% 1D objects, so that their thickness in a given domain stays constant.
%
% For each clay layer in the FW and HW, the iterative loop takes into
% account the maximum length of a smear segment, indicated as an input 
% (smear.SegLenMax), as well as M.Psmear. Then, it places smear objects of
% maximum length SegLenMax, and successively adds or removes cells until P
% within each subdomain is matched to P(smear) within a tolerance. The 
% smear objects/segments are placed one at a time within an inner for loop;
% the algorithm randomly selects a location within the longest sand segment 
% at the beginning of each iteration. This initial process is 1D, i.e. 
% using a single array of 0s and 1s, since the smear thickness in a given 
% subdomain stays constant.
% 
% Visual example of smear segment placement:
% - Initially, there is 0 smear along the path (diagonal):
%   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
% - Iteration 1
%   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 
%   (For this case, we see that each smear segment has a maximum length
%   corresponding to 6 cells. This depends on SegLenMax, but may be 
%   overwritten to match Psmear in each subdomain).
%   Suppose P calculated (in 2D, with all corresponding diags) + tol <
%   Psmear. We need another iteration.
% - Iteration 2
%   0 0 1 1 1 1 1 1 0 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0
%   Suppose P calculated - tol > Psmear, we need another iteration.
% - Now remove smear cells one (smear row) at a time
%   0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0
%   until abs(Psmear - P calc.) < tol. 
%
%
% -----------------------------INPUTS-------------------------------------%
% M     = mapping matrix structure. For arrays, the FW goes first (left) 
%         and then the HW. See faultMaterialMap.m for full documentation.
% smear = smear object with the corresponding fields (see Smear.m)
% FS    = faulted section object (see FaultedSection.m)
% G     = MRST grid structure (see MRST documentation).
% tolerance = tolerance for matching P calculated to Psmear, typically
%             around 0.01-0.05.
% verbose = 1 for displaying info on smear placement (verbose), 0 to avoid 
%           it. Optional argument, default value is set to 0 (preferred if
%           many simulations are going to be run).
%
%
% -----------------------------OUTPUT-------------------------------------
% M     = mapping matrix structure. For arrays, the FW goes first (left) 
%         and then the HW. See faultMaterialMap.m for full documentation.
%         This function updates the field M.vals with the correct number
%         and placement of 0s (sand) and 1s (smear) after object-based
%         simulation. M.vals can be directly mapped to the grid cell 
%         indexes in G by using (this is for 2D grid):
%         >>  gridIndices = reshape(transpose(flipud(M.vals)), ...
%                                   G.cells.num, 1);
%         We need to flip up-down and transpose since G indexing starts at 
%         bottom left, columns (x) move faster. MATLAB matrix indexing 
%         starts counting at top left and rows (y) move faster.
%_________________________________________________________________________
 
% verbose
if nargin < 6
    verbose = 0;
end

% cell_union_psmear: realize each clay source independently (reusing the legacy
% single-band placement below via a recursive call), union the active masks,
% and label cells. Returns a final per-cell map (M.cellBasedFinal = true). 
if isfield(M, 'cellUnion') && M.cellUnion
    M = placeSmearObjectsCellUnion(M, smear, FS, G, tolerance, verbose);
    return
end
 
% Options
opt_merge_thres = 0.4; % Merge smear segments shorter than this fraction of the maximum segment length, scaled.
 
% Get smear lengths and maximum segment lengths for *discontinuous* smear 
% domains in fault, in the order in which they appear in the fault (both 
% the order and number of smears in the fault may be different from the 
% strati due to overlaps)
smear_lens = zeros(1, max(M.unitIn));    
smear_lens(M.isclayIn) = smear.Length;   
smear_lens = smear_lens(M.unit);              
smear_lens(smear_lens == 0) = [];                                           % all in fault
smear_lens = smear_lens(M.Psmear<1);                                        % in fault & discontinuous
 
smear_seg_max_lens = zeros(1, max(M.unitIn)); 
smear_seg_max_lens(M.isclayIn) = smear.SegLenMax;
smear_seg_max_lens = smear_seg_max_lens(M.unit);
smear_seg_max_lens(smear_seg_max_lens == 0) = [];
smear_seg_max_lens = smear_seg_max_lens(M.Psmear<1);                        % in fault & discontinuous 
 
% Define position-smear-variables from M for easy access
p = M.Psmear(M.Psmear<1);
p_fin = zeros(1, numel(M.Psmear)); 
p_fin(M.Psmear==1) = 1;
id_smear = find(M.isclay); id_smear = id_smear(M.Psmear<1);                 % Smears with Psmear < 1
 
smear_diag_bound = [M.DiagBot(id_smear);M.DiagTop(id_smear)];
smear_diag_main  = min(abs(smear_diag_bound));                              % Main domain diagonals: Ndiag to subtract from G.cartDims(1)
if any(all([smear_diag_bound(1,:)<0; smear_diag_bound(2,:)>0]))
    id_layer = all([smear_diag_bound(1,:)<0; smear_diag_bound(2,:)>0]);
    smear_diag_main(id_layer) = 0;
end
smear_ndiag = M.nDiag(id_smear);                                            % N of diags of each clay unit (Psmear < 1)
smear_win_bot = M.windowBot(id_smear);     
smear_win_top = M.windowTop(id_smear);
 
% Check & revise tolerance value
if tolerance < 0.01
    tolerance = 0.01;
    disp('Tolerance for probability matching is too low. Re-set to 0.01')
elseif tolerance > 0.2
    tolerance = 0.2;
    disp('Tolerance for probability matching is too high. Re-set to 0.2')
end
 
% We need to do the process for each smear with Psmear < 1.
p_obj = zeros(1,numel(id_smear));
for j = 1:numel(id_smear)
    
    % Verbose output
    if verbose == 1 
        disp(['Material subdomain ' num2str(id_smear(j)) ' is discontinuous smear.']) 
    end
    
    % Get max number of diag entries (cells) within smear window, and
    % number of available cells for smear for the domain j in M.
    smear_parent = M.unit(id_smear(j));
    smear_diag_bounds = [smear_diag_bound(1,j) smear_diag_bound(2,j)];
    smear_win = [smear_win_bot(j) smear_win_top(j)];
    [diag_cells_num, n_cells_domain] = getCellCounts(G, FS, smear_parent, smear_diag_bounds, smear_win, ...
                                                     smear_diag_main(j), smear_ndiag(j));
    
    % Adjust maximum smear segment length if longer than grid
    cell_diag_len  = sqrt(sum(G.cellDim(~isnan(G.cellDim).^2)));
    smear_seg_max_lens(j) = min([smear_seg_max_lens(j), smear_lens(j)]);
    if smear_seg_max_lens(j) > G.cartDims(end)*cell_diag_len
        smear_seg_max_lens(j) = G.cartDims(end)*cell_diag_len;
    end
    
    % Estimate number of cells corresponding to smearSegMaxLen
    smear_seg_ncell_max = min(round(smear_seg_max_lens(j)*G.cartDims(end)/smear.DomainLength), diag_cells_num);
    smear_seg_ncell_max(smear_seg_ncell_max == 0) = 1;
    
    % Iterate to end up with number of cells with smear = P(smear).
    % This outer iterative loop matches the smear fraction to a value that 
    % is either within the indicated tolerance, or to a value that is 
    % higher than P(smear)
    ncell_merge_thres   = round(opt_merge_thres*(smear_seg_max_lens(j)*diag_cells_num/smear.DomainLength));
    smear_cells         = [];
    p_obj_it            = 0;
    max_its             = 10;
    i                   = 0;
    sand_cells          = (1:diag_cells_num)';
    while i < max_its
        
        % Updates and initial checks
        i = i + 1;
        smear_cells_toadd = (p(j)-p_obj_it)*diag_cells_num;
        smear_num = round(smear_cells_toadd/smear_seg_ncell_max);
        
        if smear_cells_toadd < 1
            smear_cells_toadd = round(smear_cells_toadd);
            if smear_cells_toadd == 0
                if verbose == 1
                    disp(['Tolerance of ' num2str(tolerance) ' cannot be met: not enough cell resolution.'])
                end
                break
            end
        end
        
        if smear_num == 0
            smear_num = 1;
            smear_seg_ncell_max = ceil(smear_cells_toadd);  % ceil to minimize its.
        end
        
        % Place smear segments if sand cells available
        if numel(smear_cells) < diag_cells_num
            for k = 1:smear_num              
                [sand_cells, smear_cells] = placeSmearSegment(sand_cells, smear_cells, smear_seg_ncell_max, ...
                                                              ncell_merge_thres, diag_cells_num);
            end
        end
        
        % Initialize logical vector representing the main diagonal. This will be
        % repmat to the total number of diagonals in the domain later.
        diag_is_smear = false(diag_cells_num, 1);
        diag_is_smear(smear_cells) = true;
        
        % Force split of segments of length > segLMax, which may occur from the
        % random placing of smears.
        if numel(smear_cells) < diag_cells_num
            [diag_is_smear, smear_cells] = splitSmearSegments(smear, G, diag_cells_num, diag_is_smear, ...
                                                              smear_seg_max_lens(j));
        end
        
        % Expand diag_is_smear to domain (vector to matrix)
        [diag_is_smear_domain, diag_ids] = getDiagIsSmearDomain(G, FS, smear_diag_bound(:, j), smear_parent, ...
                                                                smear_win, smear_ndiag(j), diag_is_smear);
        
        % Compute iteration probability
        M_test      = transpose(full(spdiags(diag_is_smear_domain, diag_ids, false(G.cartDims(1)))));
        p_obj_it    = sum(sum(M_test)) / sum(n_cells_domain);
        %figure(9); spy(Mtest); pause           % to debug
        
        % Compare with desired probability p(j)
        if abs(p(j) - p_obj_it) < tolerance 
            break   % already within tolerance, so we are done for this domain.
   
        elseif p_obj_it > p(j)
            % If it is larger than true value, we remove cells one row at 
            % a time (iteratively), until P is within tolerance. Alternatively, 
            % if P becomes smaller than true value again, it means that true 
            % P cannot be matched with this cell resolution, so we break anyway.
            % Cells are removed randomly from either the top or the bottom 
            % of the longest segment.
            
            if verbose == 1
                disp('2nd iterative loop (removal) required to match Psmear in the 2D subdomain:')
            end
            
            [diag_is_smear_domain, p_obj_it] = adjustSmearCells(G, FS, smear_cells, diag_cells_num, ...
                                                                smear_diag_bound(:, j), smear_win, smear_parent, ...
                                                                smear_ndiag(j), diag_ids, n_cells_domain, p(j), ...
                                                                tolerance, verbose, 'remove');
            break
        
        elseif i == max_its && abs(p(j) - p_obj_it) > tolerance 
            % Less smear than required. Add one cell at a time, picking a
            % random cell with 0 that is attached to a smear (avoid placing
            % extra short smear segments). Shortest smear is chosen.
            
            if numel(smear_cells) == diag_cells_num
                if verbose == 1
                    disp('Psmear cannot be matched, not enough cell resolution')
                end
                
            else
                if verbose == 1
                    disp('2nd iterative loop (addition) required to match Psmear in the 2D subdomain:')
                end
                
                [diag_is_smear_domain, p_obj_it] = adjustSmearCells(G, FS, smear_cells, diag_cells_num, ...
                                                                 smear_diag_bound(:, j), smear_win, smear_parent, ...
                                                                 smear_ndiag(j), diag_ids, n_cells_domain, p(j), ...
                                                                 tolerance, verbose, 'add');
                
                break
            end
        end
    end
 
    % Final messages
    if verbose == 1
        disp(['SSFc-derived P = ' num2str(p(j)) '. Final smear ' num2str(j) ' P = ' ...
              num2str(p_obj_it) ' (' num2str(i) ' main iterations)'])
    end
    
    if p_obj_it > 0
        if verbose == 1 
            numSmearSegFin = ceil(numel(find(diff(diag_is_smear_domain(:,1))))*0.5);
            disp(['Input smear length [m] = ' num2str(smear_seg_max_lens(j)) ...
                  '. Number of final smear segments = ' num2str(numSmearSegFin)]);
        end
        
    else
        if verbose == 1
            disp(['Input smear length [m] = ' num2str(smear_seg_max_lens(j)) ...
                  '. Number of final smear segments = 0']);
        end
        
        diag_is_smear = false(diag_cells_num, 1);
        [diag_is_smear_domain, diag_ids] = getDiagIsSmearDomain(G, FS, smear_diag_bound(:, j), smear_parent, ...
                                                                smear_win, smear_ndiag(j), diag_is_smear, true); 
    end
    
    if verbose == 1
        disp('___________________________________________________')
    end
    
    %  Assign values to mapping matrix & save final smear probability
    M.vals = full(spdiags(diag_is_smear_domain, diag_ids, M.vals));
    p_obj(j) = p_obj_it;
    
end
 
% To compromise spdiags truncation rules with smears chopped parallel
% to f.T (what we want) and not to f.L, the last step is to transpose. 
% Note that this needs to be done since the diagonal indices (diagIds) used 
% in spdiags during M.vals and Mtest assignment are flipped. Mtest is 
% directly transposed above, so that P is computed and matched with the 
% same exact configuration as in the definitive mapping matrix (M.vals).
M.vals = transpose(M.vals);
 
% Pass desired and final probabilities to output
p_fin(M.Psmear < 1) = p_obj;             
M.P = [M.Psmear; p_fin];               % [desired; obtained] (incl. continuous smears)
   
   
%% Helper functions
    function [diag_cells_num, n_cells_domain] = getCellCounts(G, FS, smear_parent, smear_diag_bounds, ...
                                                              smear_win, smear_diag_main, smear_ndiag)
        %
        % Get max number of diag entries (cells) within smear window (diag_cells_num),
        % and number of available cells for smear for the domain j in M (n_cells_domain)
        %
        if ismember(smear_parent, FS.HW.Id) && smear_diag_bounds(1) < 0 || ...
           ismember(smear_parent, FS.FW.Id) && smear_diag_bounds(2) > 0
            diag_cells_num = numel(smear_win(1):smear_win(2)) - smear_diag_main;
            
            if all([smear_diag_bounds(1)<0 smear_diag_bounds(2)>0]) && ismember(smear_parent, FS.HW.Id)
                n_cells_domain = [diag_cells_num+smear_diag_bounds(1):diag_cells_num ...
                    fliplr(G.cartDims(1)-smear_diag_bounds(2):G.cartDims(1)-1)];
                
            elseif all([smear_diag_bounds(1)<0 smear_diag_bounds(2)>0]) && ismember(smear_parent, FS.FW.Id)
                n_cells_domain = [G.cartDims(1)+smear_diag_bounds(1):G.cartDims(1) ...
                    fliplr(diag_cells_num-smear_diag_bounds(2):diag_cells_num-1)];
                
            elseif abs(smear_diag_bounds(1)) >= abs(smear_diag_bounds(2))
                n_cells_domain = diag_cells_num-(smear_ndiag - 1):diag_cells_num;
            else
                n_cells_domain = fliplr(diag_cells_num-(smear_ndiag - 1):diag_cells_num);
            end
            
        else
            diag_cells_num = numel(smear_win(1):smear_win(2));
            if abs(smear_diag_bounds(1)) > abs(smear_diag_bounds(2))
                n_cells_domain = G.cartDims(1)-abs(smear_diag_bounds(1)):G.cartDims(1)-abs(smear_diag_bounds(2));
            else
                n_cells_domain = fliplr(G.cartDims(1)-abs(smear_diag_bounds(2)):G.cartDims(1)-abs(smear_diag_bounds(1)));
            end
        end
        n_cells_domain(n_cells_domain > diag_cells_num) = diag_cells_num;
        assert(~isempty(n_cells_domain))
        assert(all(n_cells_domain >= 0)) 
    end
 
 
    function [id, seg_cells] = findSegment(cell_ids, type)
        %
        % Find minimum or maximum number of consecutive entries (segment)
        % in a vector containing cell_ids
        %
        id_breaks    = find(diff(cell_ids) > 1)';
        id_start_end = [unique([1, id_breaks+1]); 
                        unique([id_breaks, numel(cell_ids)])];
        seg_cell_num = diff(id_start_end) + 1;
        
        if strcmp(type, "longest")
            [~, id] = max(seg_cell_num);
        elseif strcmp(type, "shortest")
            [~, id] = min(seg_cell_num);
        end
        
        seg_cells = cell_ids(id_start_end(1,id):id_start_end(2,id));
        
    end
 
 
    function [sand_cells, smear_cells] = placeSmearSegment(sand_cells, smear_cells, smear_seg_ncell_max, ...
                                                           ncell_merge_thres, diag_cells_num)
        %
        % Places a smear segment in a random position within the longest
        % sand segment.
        %
        
        % Cells of longest sand segment (where smear segment will be placed)
        [~, sand_seg_cells] = findSegment(sand_cells, "longest");
        
        % Determine where to add smear cells
        short_segment = smear_seg_ncell_max == 1 && ~isempty(smear_cells) || ...
                        smear_seg_ncell_max < ncell_merge_thres;
                    
        if short_segment
            % Attach them to existing smear segments
            if numel(sand_seg_cells) <= smear_seg_ncell_max
                possible_cell_starts = sand_seg_cells(1);
            else
                possible_cell_starts = sand_seg_cells(1):(sand_seg_cells(end) - smear_seg_ncell_max-1);
            end
            
            min_start = min(possible_cell_starts);
            max_start = max(possible_cell_starts);
            if min_start == 1
                % Attach at bottom
                cells_to_add = max_start:(max_start + smear_seg_ncell_max-1);
            elseif (max_start + (smear_seg_ncell_max-1)) == max(sand_cells)
                % Attach at top
                cells_to_add = min_start:(min_start + smear_seg_ncell_max-1);
            else
                % Random
                if rand(1) > 0.5
                    cells_to_add = min_start:(min_start + smear_seg_ncell_max-1);
                else
                    cells_to_add = max_start:(max_start + smear_seg_ncell_max-1);
                end
            end
            
        else
            % Randomly placed, each cell has equal probability
            rand_cell_start = sand_seg_cells(randi(numel(sand_seg_cells), 1));
            num_cells_available = (max(sand_seg_cells) - rand_cell_start) + 1;
            
            if num_cells_available < smear_seg_ncell_max % wrap to top
                cells_bot = rand_cell_start:sand_seg_cells(end);
                cells_top = sand_seg_cells(1):(sand_seg_cells(1) + smear_seg_ncell_max - (numel(cells_bot)+1));
                cells_to_add = [cells_top cells_bot];
            else
                cells_to_add = rand_cell_start:(rand_cell_start + smear_seg_ncell_max-1);
            end
        end
        
        % Updates
        smear_cells = unique([cells_to_add'; smear_cells]);
        sand_cells = (1:diag_cells_num)';
        if max(sand_cells) < max(smear_cells)
           smear_cells(smear_cells > max(sand_cells)) = [];
        end
        sand_cells(smear_cells) = [];
        
    end
 
 
    function [diag_is_smear, smear_cells] = splitSmearSegments(smear, G, diag_cells_num, diag_is_smear, smear_seg_max_len)
        %
        % Splits smear segments if their length is larger than
        % smear_seg_max_len. Currently splits all segments in the center.
        %
        length_avail = smear.DomainLength*(diag_cells_num/G.cartDims(end));
        [id_start, seg_ncells, smear_seg_lens, long_smear_segments] = ...
            findLongSmearSegments(diag_is_smear, length_avail, diag_cells_num, smear_seg_max_len);
        
        max_its_div = 5;
        it_div = 1;
        while any(long_smear_segments)
            ids = smear_seg_lens > smear_seg_max_len;
            smear_seg_centers = id_start(ids) + fix(seg_ncells(ids)/2);
            diag_is_smear(smear_seg_centers) = false;
            
            [id_start, seg_ncells, smear_seg_lens, long_smear_segments] = ...
                findLongSmearSegments(diag_is_smear, length_avail, diag_cells_num, smear_seg_max_len);
            
            if all(~long_smear_segments)
                break
            elseif it_div == max_its_div
                error('max num it. reached, check smear segment length.')
            else
                it_div = it_div + 1;
            end
            
        end
        
        smear_cells = find(diag_is_smear);  % update smearCells
 
        % Small helper to avoid code repetition
        % Threshold set to 1.05
        function [id_start, ncells, seg_lens, long_seg] = findLongSmearSegments(diag_is_smear, length_avail, ...
                                                                                diag_cells_num, smear_seg_max_len)
            threshold_factor = 1.05;
            id_start = strfind([0, diag_is_smear'],[0 1]);
            id_end = strfind([diag_is_smear', 0],[1 0]);
            ncells = id_end - id_start + 1;
            seg_lens = ncells*(length_avail/diag_cells_num);
            long_seg = seg_lens > threshold_factor * smear_seg_max_len;
        end
        
    end
 
 
    function [diag_is_smear_domain, diag_ids] = getDiagIsSmearDomain(G, FS, smear_diag_bound, smear_parent, ...
                                                                     smear_win, smear_ndiag, diag_is_smear, is0)
        %
        % Expands diag_is_smear to the whole domain. 
        % If called with 7 args, expands diag_is_smear.
        % If called with 8th arg is0==true, fills with zeros.
        %
        if nargin == 7 
            if smear_diag_bound(1) < 0 && ismember(smear_parent, FS.HW.Id)
                add_zeros_bot = zeros(smear_win(1)-1, 1);
                add_zeros_top = zeros(G.cartDims(end)-numel([diag_is_smear; add_zeros_bot]), 1);
                diag_is_smear_domain = repmat([add_zeros_top; diag_is_smear; add_zeros_bot], 1, smear_ndiag);
 
            elseif smear_diag_bound(2) > 0 && ismember(smear_parent, FS.FW.Id)
                add_zeros_top = zeros(G.cartDims(end)-smear_win(2), 1);
                add_zeros_bot = zeros(G.cartDims(end)-numel([add_zeros_top; diag_is_smear]), 1);
                diag_is_smear_domain = repmat([add_zeros_top; diag_is_smear; add_zeros_bot], 1, smear_ndiag);
 
            elseif -smear_diag_bound(2) >= 0
                add_zeros_h = zeros(G.cartDims(1)-numel(diag_is_smear), 1);
                diag_is_smear_domain = repmat([add_zeros_h; diag_is_smear], 1, smear_ndiag);
 
            else
                add_zeros_h = zeros(G.cartDims(1)-numel(diag_is_smear), 1);
                diag_is_smear_domain = repmat([diag_is_smear; add_zeros_h], 1, smear_ndiag);
            end
 
            num_low_h = -smear_diag_bound(2); 
            num_up_h = -smear_diag_bound(1);
            diag_ids = num_low_h:num_up_h;
            
        elseif nargin == 8 && is0 % fill with 0s
            num_low = -smear_diag_bound(2); 
            num_up = -smear_diag_bound(1);
            add_zeros = zeros(G.cartDims(1) - numel(diag_is_smear), 1);
            diag_is_smear_domain = repmat([add_zeros; diag_is_smear], 1, smear_ndiag);
            diag_ids = num_low:num_up;
        end
        
    end
 
    function [diag_is_smear_domain, p_obj_it, smear_cells] = adjustSmearCells(G, FS, smear_cells, diag_cells_num, ...
                                                           smear_diag_bound, smear_win, smear_parent, smear_ndiag, ...
                                                           diag_ids, n_cells_domain, p, tolerance, verbose, mode)
        %
        % Remove smear cells, one row at a time, to match smear probability.
        % Selects the longest smear segment to remove cells.
        %
        
        iter  = 0;
        MAXITS = 25;
        rval = round(rand(1));        % only once per smear for consistency
        
        while iter < MAXITS
            % Update iteration number
            iter = iter + 1;
            
            % Update smear cells
            if strcmp(mode, 'remove')
                [~, long_seg_cells] = findSegment(smear_cells, 'longest');
                if rval == 1
                    smear_cells(smear_cells == long_seg_cells(1))= []; 
                else
                    smear_cells(smear_cells == long_seg_cells(end))= [];
                end      
            
            elseif strcmp(mode, 'add')
                [~, short_seg_cells] = findSegment(smear_cells, 'shortest');
                r = rand(1);
                if short_seg_cells(end) == diag_cells_num || short_seg_cells(1) > 1 && r >= 0.5
                    cell_to_add = short_seg_cells(1) - 1;
                elseif short_seg_cells(1) == 1 || short_seg_cells(end) < diag_cells_num && r < 0.5
                    cell_to_add = short_seg_cells(end) + 1;
                end
                smear_cells = unique([smear_cells; cell_to_add]);
            
            else
                error('Unknown mode: %s', mode);  
            end
            
            % Update logical vector of main diagonal
            diag_is_smear_var = false(diag_cells_num,1);
            diag_is_smear_var(smear_cells)  = true;
            
            % Expand diag_is_smear to the full domain
            diag_is_smear_domain = getDiagIsSmearDomain(G, FS, smear_diag_bound, smear_parent, ...
                                                        smear_win, smear_ndiag, diag_is_smear_var);
            
            % Compute updated P
            M_test_var = transpose(full(spdiags(diag_is_smear_domain, diag_ids, false(G.cartDims(1)))));
            p_obj_it = sum(M_test_var(M_test_var==true)) / sum(n_cells_domain);
            
            % Break if finished and verbose messages
            if abs(p_obj_it - p) < tolerance
                if verbose == 1
                    disp(['Tolerance of ' num2str(tolerance) ' met in ' num2str(iter) ' iterations.'])
                end
                break
                
            elseif iter == MAXITS
                if verbose == 1
                    disp('Tolerance cannot be met: max number of iterations reached.')
                end
                
            elseif (strcmp(mode,'remove') && p_obj_it < p) || (strcmp(mode,'add') && p_obj_it > p)
                if verbose == 1
                    disp(['Tolerance of ' num2str(tolerance) ' cannot be met: not enough cell resolution.'])
                    disp(['Number of iterations in loop= ' num2str(iter)])
                end
                break
            end
        end
        
    end
    
end

% ======================================================================
%  cell_union_psmear placement (deferred from faultMaterialMap)
% ======================================================================
% These local functions implement the per-source realization, cell-wise
% union, and parent-unit labeling for the 'cell_union_psmear' rule. The
% legacy object-placement core above is reused unchanged: each clay source
% is realized by a recursive placeSmearObjects call on a single-domain M
% (which has no cellUnion field, so it runs the legacy path).

function M = placeSmearObjectsCellUnion(M, smear, FS, G, tolerance, verbose)
    % Realize each clay source's discontinuous smear independently, take the
    % cell-wise union of the active masks, then label every cell with a parent
    % unit: overlapping active clay cells get the nearest clay source; non-clay
    % cells get the nearest sand source. Produces a final per-cell map and sets
    % M.cellBasedFinal = true.
    
    n  = G.cartDims(1);          % mapping-matrix dimension (square grid)
    nz = G.cartDims(end);
    
    % Diagonal index of each cell in the final (transposed) orientation, i.e. the
    % same orientation as the legacy post-placement M.vals.
    [rowId, colId] = ndgrid(1:n, 1:n);
    cellDiag = rowId - colId;
    
    nSrc     = numel(M.unit);    % clay-source descriptor list from faultMaterialMap
    clayMask = false(n, n);
    unitMap  = zeros(n, n);
    bestDist = inf(n, n);
    actualP  = zeros(1, nSrc);   % realized coverage per source (aligned to M.unit)
    
    for s = 1:nSrc
        unitId = M.unit(s);
    
        % Skip empty / zero-probability sources
        if M.nDiag(s) <= 0 || ~(M.Psmear(s) > 0)
            continue
        end
    
        % Cells this source could occupy: diagonal band intersected with its
        % elevation window (window precomputed in faultMaterialMap).
        candidateMask = sourceCandidateMask(n, nz, M.DiagBot(s), M.DiagTop(s), ...
            M.windowBot(s), M.windowTop(s));
        if ~any(candidateMask(:))
            continue
        end
    
        % Realize this source's smear. Continuous sources fill the candidate;
        % discontinuous sources are realized by the legacy object placement.
        if M.Psmear(s) >= 1
            sourceMask = candidateMask;
            actualP(s) = 1;
        else
            [sourceMask, actualP(s)] = realizeSourcePsmear(M, smear, FS, G, s, ...
                candidateMask, ...
                tolerance, verbose);
        end
        if ~any(sourceMask(:))
            continue
        end
    
        % Union, and assign parent = nearest clay source (ties: lower unit id).
        % NOTE (suspected sign error, from prototype). This only affects PARENT LABELING in OVERLAP cells
        % (which source's poro/perm/Vcl a shared cell gets); Left as-is to match Shaowen's
        % prototype; verify and switch to the line below during the MATLAB pass.
        srcDist = abs(cellDiag - M.layerDiagCenter(unitId));
        % srcDist = abs(cellDiag + M.layerDiagCenter(unitId));   % suspected fix
    
        tie     = abs(srcDist - bestDist) <= 1e-12 & (unitMap == 0 | unitId < unitMap);
        upd     = sourceMask & (srcDist < bestDist | tie);
    
        clayMask      = clayMask | sourceMask;
        unitMap(upd)  = unitId;
        bestDist(upd) = srcDist(upd);
    end
    
    % Fill non-clay cells with the nearest sand parent.
    sandMap = nearestSandUnitMap(M, cellDiag);
    if any(~clayMask(:)) && any(sandMap(~clayMask) == 0)
        error(['cell_union_psmear requires at least one sand parent unit when ' ...
            'some cells are not occupied by clay smear.'])
    end
    unitMap(~clayMask) = sandMap(~clayMask);
    if any(clayMask(:) & unitMap(:) == 0)
        error('cell_union_psmear failed to assign a parent unit to some clay cells.')
    end
    
    % Finalize per-cell map. Orientation matches the legacy post-placement M.vals
    % (transposed once), so setGridPoroPerm's transpose(flipud(.)) maps correctly.
    M.vals  = clayMask;
    M.units = unitMap;
    M.cellBasedFinal = true;
    
    % Rebuild domain descriptors from the realized map (band fields NaN; M.unit
    % lists every unit actually present in the fault).
    M = rebuildCellUnionDescriptors(M, actualP);

end


function mask = sourceCandidateMask(n, nz, diagBot, diagTop, idBot, idTop)
    % Cells reachable by one clay source: its diagonal band [diagBot, diagTop]
    % intersected with its elevation-window rows [idBot, idTop]. Built in the same
    % transposed-once orientation as the legacy M.vals.
    
    nDiag = diagTop - diagBot + 1;
    if nDiag <= 0
        mask = false(n, n);
        return
    end
    diagVals = false(nz, nDiag);
    diagVals(idBot:idTop, :) = true;
    diagVals = flipud(diagVals);                       % for spdiags (as in legacy)
    maskPreT = full(spdiags(diagVals, -diagTop:-diagBot, false(n)));
    mask = transpose(maskPreT);
end


function [activeMask, actualP] = realizeSourcePsmear(M, smear, FS, G, s, candidateMask, tol, verbose)
    % Realize one discontinuous clay source using the legacy single-band object
    % placement: build a one-domain M and call placeSmearObjects recursively. The
    % sub-call has no cellUnion field, so it runs the legacy path. The realized
    % mask is clipped to the source's candidate cells.
    
    Msrc           = struct();
    Msrc.vals      = false(size(M.vals));
    Msrc.units     = zeros(size(M.units));
    Msrc.unit      = M.unit(s);
    Msrc.isclay    = true;
    Msrc.unitIn    = M.unitIn;
    Msrc.isclayIn  = M.isclayIn;
    Msrc.nDiag     = M.nDiag(s);
    Msrc.DiagBot   = M.DiagBot(s);
    Msrc.DiagTop   = M.DiagTop(s);
    Msrc.Psmear    = M.Psmear(s);
    Msrc.windowBot = M.windowBot(s);
    Msrc.windowTop = M.windowTop(s);
    
    Msrc = placeSmearObjects(Msrc, smear, FS, G, tol, verbose);
    
    activeMask = logical(Msrc.vals) & candidateMask;
    actualP    = sum(activeMask(:)) / max(1, sum(candidateMask(:)));
end


function unitMap = nearestSandUnitMap(M, cellDiag)
    % Assign each cell to the closest sand source layer (ties: lower unit id).
    % NOTE (possible sign error, from prototype): same convention as the clay
    % labeling, the commented line is the suspected fix to verify during the MATLAB pass.
    
    sandUnits = M.unitIn(~M.isclayIn);
    unitMap = zeros(size(cellDiag));
    if isempty(sandUnits)
        return
    end
    bestDistance = inf(size(cellDiag));
    for i = 1:numel(sandUnits)
        unitId  = sandUnits(i);
        srcDist = abs(cellDiag - M.layerDiagCenter(unitId));
        % srcDist = abs(cellDiag + M.layerDiagCenter(unitId));   % fix?
        tie     = abs(srcDist - bestDistance) <= 1e-12 & ...
                  (unitMap == 0 | unitId < unitMap);
        update  = srcDist < bestDistance | tie;
        unitMap(update)      = unitId;
        bestDistance(update) = srcDist(update);
    end
end


function M = rebuildCellUnionDescriptors(M, actualP)
    % Rebuild M.unit / M.isclay / Psmear etc. from the realized per-cell M.units,
    % enumerating every unit actually present in the fault (clay sources that
    % labeled >= 1 cell, plus every sand unit that labeled a cell). Band fields are
    % set to NaN: cell-union domains do not form non-overlapping diagonal bands.
    
    % Map per-source quantities (aligned to the clay-source list M.unit) to
    % per-unit-id arrays BEFORE M.unit is overwritten.
    pSmearByUnit = zeros(1, max(M.unitIn));
    actualByUnit = zeros(1, max(M.unitIn));
    pSmearByUnit(M.unit) = M.Psmear;
    actualByUnit(M.unit) = actualP;
    inputClayUnits = M.unit;                 % clay sources we attempted
    
    % Units present in the realized map
    appearing = unique(M.units(:))';
    appearing(appearing == 0) = [];
    appearingClay = unique(M.units(logical(M.vals)))';
    appearingClay(appearingClay == 0) = [];
    
    % Descriptor arrays for all present units (bands undefined in this mode)
    M.unit           = appearing;
    M.isclay         = ismember(M.unit, M.unitIn(M.isclayIn));
    M.nDiag          = nan(1, numel(M.unit));
    M.DiagBot        = nan(1, numel(M.unit));
    M.DiagTop        = nan(1, numel(M.unit));
    M.windowTop      = nan(1, numel(M.unit));
    M.windowBot      = nan(1, numel(M.unit));
    M.unitInClayGaps = nan(1, numel(M.unit));
    M.clayDiagBot    = [];
    M.divLayerDiag   = [];
    
    % Probabilities for clay units that actually label cells
    M.Psmear             = pSmearByUnit(appearingClay);
    M.P                  = [M.Psmear; actualByUnit(appearingClay)];
    M.PsmearByUnit       = pSmearByUnit;
    M.actualPsmearByUnit = actualByUnit;
    
    % Clay sources that produced no labeled clay cells (removed by overlap/union)
    M.idSmearInRemoved = find(~ismember(inputClayUnits, appearingClay));
end