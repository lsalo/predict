function M = placeSmearObjects(M, G, tolerance, smearLength, ...
                               smearSegLenMax, smearDomainLen, ...
                               faultDisp, faultThick, verbose)
%
% -----------------------------SUMMARY------------------------------------
% This function randomly places smear objects using object-based simulation, 
% so that each realization provides a random case consistent with the
% inputs. The final goal is to assign correct number of 0s and 1s in M.vals 
% such that probabilities in M.Psmear are matched to the indicated 
% tolerance. These probabilities are matched in 2D, ie accounting for the 
% total number of cells in each subdomain. Smears, however, are treated as 
% 1D objects, so that their thickness in a given subdomain stays constant.
%
% For each clay layer in the FW and HW, the iterative loop takes into
% account the maximum length of a smear segment, indicated as an input 
% (smear.inL), as well as P(smear). Then, it places smear objects of
% maximum length inL, and successively adds or removes cells until P(1)
% within each subdomain is matched to P(smear) to the given tolerance. The 
% smear objects/segments are placed one at a time within an inner for loop;
% the algorithm randomly selects a location within the longest sand segment 
% at the beginning of each iteration. This process is done in 1D, i.e. 
% using a single array of 0s and 1s, since the smear thickness in a given 
% subdomain stays constant.
% 
% Visual example of smear segment placement:
% - Initially, there is 0 smear along the path (diagonal):
%   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
% - Iteration 1
%   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 
%   (For this case, we see that each smear segment has a maximum length
%   corresponding to 6 cells. This depends on the user input (desired
%   maximum length) but is limited by the fraction of smear that can be
%   present in each subdomain)
%   P calculated (in 2D, with all corresponding diags) + tol < P(smear)
% - Iteration 2
%   0 0 1 1 1 1 1 1 0 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0
%   P calculated > P (smear)
% - Inner loop to remove smear cells one (smear row) at a time
%   0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0
%   abs(P(smear) - P calc.) < tol 
%
% -----------------------------INPUTS-------------------------------------
% Length units are in meters (m), and angles are in degrees (º).
%
% M     = mapping matrix structure. For arrays, the FW goes first (left) 
%         and then the HW. See smearMap.m for full documentation of fields.
%
% smear = smear structure with the following fields:
%       inL       = desired maximum smear length. Either 'Lsmear' or double
%                   value.
%       modeledAs = algorithm to use, in this case only 'objects' is 
%                   accepted.
%       tolerance = tolerance to match simulated P with P(smear). Scalar
%                   double.
%       Ts        = smear thickness for each smear (horizontal array, FW 
%                   first). 
%       Ls        = smear apparent length in a diagonal crossing the fault 
%                   from a lower extreme (eg bottom left) to the opposite 
%                   upper extreme (eg top right). Not used here.
%       L         = smear length at which it becomes discontinuous, i.e. 
%                   the upper bound for smear length as determined from 
%                   Tap*SSFc. Horizontal array.
%       Ds        = path length along which the smear is distributed.
%                   Scalar double.
%
% f     = fault structure with the following fields:
%       T     = thickness (scalar double)
%       dip   = dip angle     (")
%       t     = throw         (")
%       D     = displacement  (")
%       L     = fault length, along the fault dip, for each FW and HW layer, 
%               from top of that layer in the FW to bottom of the layer 
%               in the HW.
%       delta = delta angle, ie the angle of the shear zone within the
%               fault, along which smears are distributed, with the
%               horizontal. This would be equal to the fault dip if smears
%               were parallel to the fault.
%       alpha = 90 - gamma - delta, where gamma = 90 - f.dip.
%
% G     = MRST grid structure.
%
% verbose = 1 for displaying info on smear placement (verbose), 0 to avoid it. 
%           Optional argument, default value is set to 1 (verbose).
%
% -----------------------------OUTPUT-------------------------------------
% M     = mapping matrix structure. For arrays, the FW goes first (left) 
%         and then the HW. See smearMap.m for full documentation of fields.
%         This function updates the field M.vals with the correct number
%         and placement of 0s (sand) and 1s (smear) according to inputs.
%         M.vals can be directly mapped to the grid cell indexes in G by
%         using (this is for 2D grid):
%         >>  gridIndices = reshape(transpose(flipud(M.vals)), ...
%                                   G.cells.num, 1);
%         We need to flip up-down and transpose since G indexing starts at 
%         bottom left, columns (x) move faster. MATLAB matrix indexing 
%         starts counting at top left and rows (y) move faster.
%_________________________________________________________________________

% verbose
if nargin < 9
    verbose = 1;
end

% Add Psmear and check if  removed smears
if isfield(M, 'idSmearInRemoved')
    smearLength(M.idSmearInRemoved) = []; 
    smearSegLenMax(M.idSmearInRemoved) = [];
end

% Initial variables for easy access
Psm = M.Psmear(M.Psmear<1);
P = zeros(1,numel(M.Psmear)); 
P(M.Psmear==1) = 1;
cunitsAll = find(M.isclay);
cunits = cunitsAll(M.Psmear<1);              % Clay units with Psmear < 1
Lsmear = smearLength(M.Psmear<1);
cDiagBound = [M.DiagBot(cunits);M.DiagTop(cunits)];
cDiagMain  = min(abs(cDiagBound));           % Number of diagonals to subtract from G.cartDims(1)
if any(all([cDiagBound(1,:)<0; cDiagBound(2,:)>0]))
    idLay = all([cDiagBound(1,:)<0; cDiagBound(2,:)>0]);
    cDiagMain(idLay) = 0;
end
cnDiag = M.nDiag(cunits);                    % N of diags of each cunit (Psmear < 1)
maxLSmearSeg = smearSegLenMax(M.Psmear<1);         

% Tolerances 
tolP = tolerance;              % tolerance in deviation from P(smear)
if tolP < 1e-3
    tolP = 1e-3;
    disp('Tolerance for probability matching is too low. Re-set to 1e-3')
elseif tolP > 0.2
    tolP = 0.2;
    disp('Tolerance for probability matching is too high. Re-set to 1e-2')
end

% We need to do the process for each smear with Psmear < 1.
Pdisc = zeros(1,numel(cunits));
for j = 1:numel(cunits)
    if verbose == 1, disp(['Smear from source unit ' num2str(cunits(j)) ' is discontinuous.']), end
    
    % 2.1 Initial parameters to place smear segments in each subdomain
    smearL = min([maxLSmearSeg(j), Lsmear(j)]);
    DiagCellsNum = G.cartDims(1)-cDiagMain(j);
    cellDiagL  = sqrt(sum(G.CellDim.^2));
    if smearL > G.cartDims(1)*cellDiagL
        smearL = G.cartDims(1)*cellDiagL;
    end
    claySegCellNum = round(smearL*DiagCellsNum/smearDomainLen);
    claySegCellNum(claySegCellNum==0) = 1;
    if all([cDiagBound(1,j)<0 cDiagBound(2,j)>0])
        Dvals_Mtest = [G.cartDims(1)+cDiagBound(1,j):G.cartDims(1) ...
            fliplr(G.cartDims(1)-cDiagBound(2,j):G.cartDims(1)-1)];
    elseif abs(cDiagBound(1,j)) > abs(cDiagBound(2,j))
        Dvals_Mtest = G.cartDims(1)-abs(cDiagBound(1,j)):G.cartDims(1)-abs(cDiagBound(2,j));
    else
        Dvals_Mtest = G.cartDims(1)-abs(cDiagBound(2,j)):G.cartDims(1)-abs(cDiagBound(1,j));
    end
    assert(~isempty(Dvals_Mtest))
    
    % 2.2 Iterate to end up with number of cells with smear = P(smear).
    % This outer iterative loop matches probability to a value that is
    % either within the indicated tolerance, or to a value that is higher
    % than the "true" probability (M.Psmear).
    mergeNum   = 0.4*round(smearL*DiagCellsNum/smearDomainLen);
    smearCells = [];
    Pfin       = 0;
    maxIts     = 10;
    itNum      = 0;
    sandCells  = 1:DiagCellsNum;
    m          = 0;
    while itNum < maxIts
        itNum        = itNum + 1;
        smearCellsT  = ((Psm(j)-Pfin)*DiagCellsNum);
        smearNum     = round(smearCellsT/claySegCellNum);
        if smearCellsT < 1
            smearCellsT = round(smearCellsT);
            if smearCellsT == 0
                if verbose == 1, disp(['Tolerance of ' num2str(tolP) ' cannot be met: not enough cell resolution.']), end
                break
            else
                smearCellsT = 1;
            end
        end
        if smearNum == 0
            smearNum = 1;
            claySegCellNum = ceil(smearCellsT);  % avoid round to get same small changes above or below .5, ceil to minimize its.
        end
        if numel(smearCells) == DiagCellsNum
            m = m + 1;
            assert(isempty(sandCells));
            smearCells = 1:numel(smearCells)+1;
            DiagCellsNum = DiagCellsNum +1;
            warning(['diagonal cells in smear ' num2str(j) ' expanded ' ...
                     num2str(m) ' times to tolerance-match P.'])
        else
            assert(~ismember(smearNum, inf))
            % This loops over each smear segment, and places each segment
            % at a random position within the longest sand segment that
            % exists at the end of the previous iteration.
            for k = 1:smearNum              
                % Cells of sand segment where smear segment is placed
                sandLim = find(diff(sandCells) > 1);
                sandBounds = [1 sandLim+1; sandLim numel(sandCells)];
                sandCellNum = diff(sandBounds)+1;
                %sandSegNum = numel(sandBounds) - 1;
                [~, sandSegId] = max(sandCellNum);
                sandSegCells = sandCells(sandBounds(1,sandSegId)):sandCells(sandBounds(2,sandSegId));
                sandSegCellNum = numel(sandSegCells);
                
                % Smear segment cells to add
                if claySegCellNum == 1 && ~isempty(smearCells) || ...
                   claySegCellNum < mergeNum
                   % Very short segments, so we attach them to existing
                   % smear segments
                    if numel(sandSegCells) <= claySegCellNum
                        possibleUpperCell = sandSegCells(1);
                    else
                        possibleUpperCell = sandSegCells(1):(sandSegCells(end) - (claySegCellNum-1));
                    end
                    if min(possibleUpperCell) == 1
                        cellsToAdd = max(possibleUpperCell):(max(possibleUpperCell)+(claySegCellNum-1));
                    elseif (max(possibleUpperCell) + (claySegCellNum-1)) == max(sandCells)
                        cellsToAdd = min(possibleUpperCell):(min(possibleUpperCell) + (claySegCellNum-1));
                    else
                        if rand(1) > 0.5
                            cellsToAdd = min(possibleUpperCell):(min(possibleUpperCell) + (claySegCellNum-1));
                        else
                            cellsToAdd = max(possibleUpperCell):(max(possibleUpperCell)+(claySegCellNum-1));
                        end
                    end
                else % Randomly placed, each cell has equal P(smear)
                    initCell = sandSegCells(randi(sandSegCellNum, 1));
                    if (max(sandSegCells) - initCell) + 1 < claySegCellNum
                        cellsToAddBot = initCell:sandSegCells(end);
                        cellsToAddTop = sandSegCells(1):(sandSegCells(1) + claySegCellNum-(numel(cellsToAddBot)+1));
                        cellsToAdd = [cellsToAddTop cellsToAddBot];
                    else
                        cellsToAdd = initCell:(initCell + (claySegCellNum-1));
                    end
                end
                
                % Updates
                smearCells = unique([cellsToAdd'; smearCells]);
                sandCells = 1:DiagCellsNum;
                if max(sandCells) < max(smearCells)
                    smearCells(smearCells > max(sandCells)) = [];
                end
                sandCells(smearCells) = [];
            end
        end
        
        % Vector values in each diagonal of each smear's domain. This is
        % for a single diagonal, will be repmat to the correct number of
        % diagonals for each smear later.
        dvals          = false(DiagCellsNum,1);
        dvals(smearCells) = true;
        
        % Divide segments of length > inL, which result from the
        % random placing of smears.
        if numel(smearCells) < DiagCellsNum
            pathLength = sqrt(faultDisp^2+faultThick^2)*(DiagCellsNum/G.cartDims(1));
            start1 = strfind([0, dvals'==1],[0 1]);
            end1 = strfind([dvals'==1, 0],[1 0]);
            segCells = end1 - start1 + 1;
            segLen = segCells*(pathLength/DiagCellsNum);
            itMaxDiv = 5;
            itDiv = 1;
            while any(segLen > 1.05*smearL)
                ids = segLen > smearL;
                segCenters = start1(ids) + fix(segCells(ids)/2);
                dvals(segCenters) = false;
                start1 = strfind([0, dvals'==1],[0 1]);
                end1 = strfind([dvals'==1, 0],[1 0]);
                segCells = end1 - start1 + 1;
                segLen = segCells*(pathLength/DiagCellsNum);
                if all(segLen <= 1.05*smearL)
                    break
                elseif itDiv == itMaxDiv
                    warning('max num It reached, check smear segment length.')
                    break
                else
                    itDiv = itDiv + 1;
                end
                    
            end
            smearCells = find(dvals);  % update smearCells
        end
        
        % Prepare for matrix and compute P(smear) = P = prob. = fraction
        addZeros = zeros(G.cartDims(1)-numel(dvals), 1);
        numLow = -cDiagBound(2,j); numUp = -cDiagBound(1,j);
        if -cDiagBound(2,j) > 0
            dvals = repmat([addZeros; dvals],1,cnDiag(j));
        else
            dvals = repmat([dvals; addZeros],1,cnDiag(j));
        end
        diagIds = numLow:numUp;
        Mtest = transpose(full(spdiags(dvals, diagIds, false(G.cartDims(1)))));
        %figure(9); spy(Mtest); pause           % visualize current smear.
        
        %Pfin  = numel(smearCells)/DiagCells;   % old 1d matching of P.
        Pfin = sum(sum(Mtest)) / sum(Dvals_Mtest);
        
        % If P already within tolerance of true value, then break.
        if abs(Psm(j) - Pfin) < tolP 
            break
         
        % If it is larger than true value, then enter second iterative
        % loop. Here we just remove cells one row at a time (iteration),
        % until P is within tolerance. Alternatively, if P becomes lower 
        % than true value again, it means that true P cannot be matched 
        % with this cell resolution. So, we break anyway.
        % Values are removed from the top of longest smear segment for 
        % smears below the main diagonal, from the bottom of longest smear
        % segment for smears above the main diagonal, and random for smears
        % with part below and part above the main diagonal.
        elseif Pfin > Psm(j)
            if verbose == 1, disp('2nd iterative loop (removal) required to match Psmear in the 2D subdomain:'), end
            itNum2  = 0;
            maxIts2 = 50;
            while itNum2 < maxIts2
                % Prepare
                itNum2 = itNum2 + 1;
                
                % Find longest smear 
                idBreaks = find(diff(smearCells) > 1)';
                idStartEnd = [unique([1, idBreaks+1]); unique([idBreaks, numel(smearCells)])];
                segCellNum = diff(idStartEnd)+1;
                [~, idLongest] = max(segCellNum);
                longSegCells = smearCells(idStartEnd(1,idLongest):idStartEnd(2,idLongest));
                
                % Remove
                if itNum2 < 2 % only once for each smear to be consistent.
                    val = round(rand(1)); % either 0 (end) or 1 (first cell)
                end
                cellRemove(val==1) = longSegCells(1);
                cellRemove(val==0) = longSegCells(end);
%                 if all(cDiagBound(:,j) < 0)  % remove from top
%                     cellRemove = longSegCells(1);
%                 elseif all(cDiagBound(:,j) > 0)  % remove from bot
%                     cellRemove = longSegCells(end);
%                 else % involves above and below main diag.
%                     if itNum2 < 2 % only once to be consistent.
%                         val = round(rand(1)); % either 0 (end) or 1 (first cell)
%                     end
%                     cellRemove(val==1) = longSegCells(1);
%                     cellRemove(val==0) = longSegCells(end);
%                 end
                smearCells(smearCells == cellRemove)= [];
                dvals2              = false(DiagCellsNum,1);
                dvals2(smearCells)  = true;
                if -cDiagBound(2,j) > 0
                    dvals2 = repmat([addZeros; dvals2],1,cnDiag(j));
                else
                    dvals2 = repmat([dvals2; addZeros],1,cnDiag(j));
                end
                Mtest = transpose(full(spdiags(dvals2, diagIds, ...
                    false(G.cartDims(1)))));
                Pfin2 = sum(Mtest(Mtest==true)) / sum(Dvals_Mtest);
                if abs(Pfin2 - Psm(j)) < tolP
                    if verbose == 1, disp(['Tolerance of ' num2str(tolP) ' met in ' num2str(itNum2) ' iterations.']), end
                    break
                elseif itNum2 == maxIts2
                    if verbose == 1, disp('Tolerance cannot be met: max number of iterations reached'), end
                elseif Pfin2 < Psm(j)
                    if verbose == 1, disp(['Tolerance of ' num2str(tolP) ' cannot be met: not enough cell resolution.']), end
                    if verbose == 1, disp(['Number of iterations = ' num2str(itNum2)]), end
                    break
                end
            end
            
            % In case of segment breaks, check and correct very short smear
            % segments that may be the result of (1) the initial divisions 
            % to avoid segment length above inL + (2) the succesive removal
            % starting from the segment top or bottom. P(smear) does not
            % need to be recomputed because the number of smear cells is
            % the same.
%             if any(diff(smearCells) > 1)
%                 idBreaks = find(diff(smearCells) > 1)';
%                 idStartEnd = [unique([1, idBreaks+1]); unique([idBreaks, numel(smearCells)])];
%                 segmentCells = diff(idStartEnd)+1;
%                 idSeg = segmentCells < mergeNum;
%                 if any(idSeg) && sum(any(idSeg))==1
%                     smearCellsOk = [];
%                     for n=1:numel(idSeg)-1
%                         idOk = find(~idSeg);
%                         smearCellsOk = [smearCellsOk; smearCells(idStartEnd(1,idOk(n)):idStartEnd(2,idOk(n)))];
%                     end
%                     if min(smearCellsOk) == 1
%                         smearCells2 = [smearCellsOk; ((smearCellsOk(end)+1):smearCellsOk(end)+segmentCells(idSeg))'];
%                     elseif max(smearCellsOk) == DiagCellsNum
%                         smearCells2 = [((smearCellsOk(1)-segmentCells(idSeg)):smearCellsOk(1)-1)'; smearCellsOk];
%                     else
%                         smearCellsBroken = smearCells(idStartEnd(1,idSeg):idStartEnd(2,idSeg));
%                         if max(smearCellsBroken) < min(smearCellsOk)
%                             smearCells2 = [((smearCellsOk(1)-segmentCells(idSeg)):smearCellsOk(1)-1)'; smearCellsOk];
%                         else
%                             smearCells2 = [smearCellsOk; ((smearCellsOk(end)+1):smearCellsOk(end)+segmentCells(idSeg))'];
%                         end
%                     end
%                     assert(numel(smearCells) == numel(smearCells2));
%                     dvals2              = false(DiagCellsNum,1);
%                     dvals2(smearCells2)  = true;
%                     if -cDiagBound(2,j) > 0
%                         dvals2 = repmat([addZeros; dvals2],1,cnDiag(j));
%                     else
%                         dvals2 = repmat([dvals2; addZeros],1,cnDiag(j));
%                     end
%                 else
%                     assert(sum(any(idSeg))==0)
%                 end
%             end
            dvals = dvals2;
            Pfin  = Pfin2;
            break
        
        elseif itNum == maxIts && abs(Psm(j) - Pfin) > tolP 
            % less smear than required. Add one cell at a time, picking a
            % random cell with 0 that is ATTACHED to a smear (avoid placing
            % extra short smear segments). Shortest smear is chosen.
            if verbose == 1, disp('2nd iterative loop (addition) required to match Psmear in the 2D subdomain:'), end
            itNum3  = 0;
            maxIts3 = 25;
            while itNum3 < maxIts3
                itNum3  = itNum3 + 1;
                dvals3   = false(DiagCellsNum,1);
                dvals3(smearCells) = true;
                start1 = strfind([0, dvals3'==1],[0 1]);
                end1 = strfind([dvals3'==1, 0],[1 0]);
                segCells = end1 - start1 + 1;
                [~, idMinL] = min(segCells);
                addToStartOrEnd = rand(1);
                if addToStartOrEnd >= 0.5 && end1(idMinL) < DiagCellsNum || ...
                   addToStartOrEnd < 0.5 && start1(idMinL) == 1
                    idToAdd = end1(idMinL) + 1;
                else
                    idToAdd = start1(idMinL) - 1;
                end
                assert(idToAdd > 0)    
                dvals3(idToAdd) = 1;
                smearCells = find(dvals3);

                if -cDiagBound(2,j) > 0
                    dvals3 = repmat([addZeros; dvals3],1,cnDiag(j));
                else
                    dvals3 = repmat([dvals3; addZeros],1,cnDiag(j));
                end
                Mtest = transpose(full(spdiags(dvals3, diagIds, ...
                                               false(G.cartDims(1)))));
                Pfin3 = sum(Mtest(Mtest==true)) / sum(Dvals_Mtest);
                if abs(Pfin3 - Psm(j)) < tolP
                    if verbose == 1, disp(['Tolerance of ' num2str(tolP) ' met in ' num2str(itNum3) ' iterations.']), end
                    break
                elseif itNum3 == maxIts3
                    if verbose == 1, disp('Tolerance cannot be met: max number of iterations reached'), end
                elseif Pfin3 > Psm(j)
                    if verbose == 1, disp(['Tolerance of ' num2str(tolP) ' cannot be met: not enough cell resolution.']), end
                    if verbose == 1, disp(['Number of iterations = ' num2str(itNum3)]), end
                    break
                end
            end
            dvals = dvals3;
            Pfin  = Pfin3;
            break
            
        end
        
    end

    if verbose == 1
        disp(['PSSF-derived P = ' num2str(Psm(j)) '. '...
        'Final smear ' num2str(j) ' P = ' ...
        num2str(Pfin) ' (' num2str(itNum) ' main iterations)'])
    end
    if Pfin > 0
        if verbose == 1 
            numSmearSegFin = ceil(numel(find(diff(dvals(:,1))))*0.5);
            disp(['Input smear length [m] = ' num2str(maxLSmearSeg(j)) ...
            '. Number of final smear segments = ' num2str(numSmearSegFin)]);
        end
        
        % 2.3 Assign 0s and 1s to matrix that maps to grid.
        M.vals = full(spdiags(dvals, diagIds, M.vals));
    else
        if verbose == 1
            disp(['Input smear length [m] = ' num2str(maxLSmearSeg(j)) ...
                '. Number of final smear segments = ' num2str(0)]);
        end
        dvals = false(DiagCellsNum,1);
        numLow = -cDiagBound(2,j); numUp = -cDiagBound(1,j);
        addZeros = zeros(G.cartDims(1)-numel(dvals), 1);
        dvals = repmat([addZeros; dvals],1,cnDiag(j));
        diagIds = numLow:numUp;
        M.vals = full(spdiags(dvals, diagIds, M.vals));
    end
    if verbose == 1, disp('___________________________________________________'), end
    
    % Save final smear probability
    Pdisc(j) = Pfin;
end

% To compromise spdiags truncation rules with smears chopped parallel
% to f.T (what we want) and not to f.L, the last step is to transpose. 
% Note that this needs to be done since the diagonal indices (diagIds) used 
% in spdiags during M.vals and Mtest assignment are flipped. Mtest is 
% directly transposed above, so that P is computed and matched with the 
% same exact configuration as in the definitive mapping matrix (M.vals).
M.vals = transpose(M.vals);

% Pass desired and final probabilities to output
P(M.Psmear<1) = Pdisc;
M.P = [M.Psmear; ...        % desired (calculated)
       P];                  % obtained

end