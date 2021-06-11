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

% Add Psmear and check if  removed/repeated smears
sLAll = zeros(1, max(M.unitIn));    segLMaxAll = zeros(1, max(M.unitIn)); 
sLAll(M.isclayIn) = smear.Length;   segLMaxAll(M.isclayIn) = smear.SegLenMax;
sLAll = sLAll(M.unit);              segLMaxAll = segLMaxAll(M.unit);
sLAll(sLAll == 0) = [];             segLMaxAll(segLMaxAll == 0) = [];


% Initial variables for easy access
pCalc = M.Psmear(M.Psmear<1);
pFin = zeros(1, numel(M.Psmear)); 
pFin(M.Psmear==1) = 1;
idcAll = find(M.isclay);
idc = idcAll(M.Psmear<1);           % Clay units with Psmear < 1
sL = sLAll(M.Psmear<1);
cDiagBound = [M.DiagBot(idc);M.DiagTop(idc)];
cDiagMain  = min(abs(cDiagBound)); % Ndiag to subtract from G.cartDims(1)
if any(all([cDiagBound(1,:)<0; cDiagBound(2,:)>0]))
    idLay = all([cDiagBound(1,:)<0; cDiagBound(2,:)>0]);
    cDiagMain(idLay) = 0;
end
cnDiag = M.nDiag(idc);              % N of diags of each cunit (Psmear < 1)
segLMax = segLMaxAll(M.Psmear<1);       
winBot = M.windowBot(idc);     
winTop = M.windowTop(idc);

% Tolerances 
tolP = tolerance;              % tolerance in deviation from pCalc
if tolP < 0.01
    tolP = 0.01;
    disp('Tolerance for probability matching is too low. Re-set to 0.01')
elseif tolP > 0.2
    tolP = 0.2;
    disp('Tolerance for probability matching is too high. Re-set to 0.2')
end

% We need to do the process for each smear with Psmear < 1.
pObj = zeros(1,numel(idc));
for j = 1:numel(idc)
    if verbose == 1, disp(['Material subdomain ' num2str(idc(j)) ...
                           ' is discontinuous smear.']), end
    
    % 1. Get max number of diag entries (cells) within smear window, 
    % initial estimate of number of cells corresponding to segLMax and
    % number of available cells for smear for the domain j in M.
    segLMax(j) = min([segLMax(j), sL(j)]);
    %diagCellsNum = G.cartDims(1)-cDiagMain(j);
    if ismember(M.unit(idc(j)), FS.HW.Id) && cDiagBound(1,j) < 0 || ...
       ismember(M.unit(idc(j)), FS.FW.Id) && cDiagBound(2,j) > 0
        diagCellsNum = numel(winBot(j):winTop(j)) - cDiagMain(j);
        if all([cDiagBound(1,j)<0 cDiagBound(2,j)>0]) && ...
           ismember(M.unit(idc(j)), FS.HW.Id) 
            dvals_Mtest = [diagCellsNum+cDiagBound(1,j):diagCellsNum ...
                           fliplr(G.cartDims(1)-cDiagBound(2,j):G.cartDims(1)-1)];
        elseif all([cDiagBound(1,j)<0 cDiagBound(2,j)>0]) && ...
               ismember(M.unit(idc(j)), FS.FW.Id) 
           dvals_Mtest = [G.cartDims(1)+cDiagBound(1,j):G.cartDims(1) ...
                          fliplr(diagCellsNum-cDiagBound(2,j):diagCellsNum-1)];
        elseif abs(cDiagBound(1,j)) >= abs(cDiagBound(2,j))
            dvals_Mtest = diagCellsNum-(cnDiag(j) - 1):diagCellsNum;
        else
            dvals_Mtest = fliplr(diagCellsNum-(cnDiag(j) - 1):diagCellsNum);
        end 
    else
        diagCellsNum = numel(winBot(j):winTop(j));
%       if all([cDiagBound(1,j)<0 cDiagBound(2,j)>0])
%          dvals_Mtest = [G.cartDims(1)+cDiagBound(1,j):G.cartDims(1) ...
%                         fliplr(G.cartDims(1)-cDiagBound(2,j):G.cartDims(1)-1)];
        if abs(cDiagBound(1,j)) > abs(cDiagBound(2,j))
            dvals_Mtest = G.cartDims(1)-abs(cDiagBound(1,j)):...
                          G.cartDims(1)-abs(cDiagBound(2,j));
        else
            dvals_Mtest = fliplr(G.cartDims(1)-abs(cDiagBound(2,j)):...
                                 G.cartDims(1)-abs(cDiagBound(1,j)));
        end
    end
    dvals_Mtest(dvals_Mtest > diagCellsNum) = diagCellsNum;
    assert(~isempty(dvals_Mtest))
    assert(all(dvals_Mtest >= 0))
    
    cellDiagL  = sqrt(sum(G.cellDim.^2));
    if segLMax(j) > G.cartDims(2)*cellDiagL
        segLMax(j) = G.cartDims(2)*cellDiagL;
    end
    %claySegCellNum = round(segLMax(j)*DiagCellsNum/smear.DomainLength);
    claySegCellNum = min(round(segLMax(j)*G.cartDims(2)/...
                               smear.DomainLength), diagCellsNum);
    claySegCellNum(claySegCellNum==0) = 1;
    
    % 2 Iterate to end up with number of cells with smear = P(smear).
    % This outer iterative loop matches the smear fraction to a value that 
    % is either within the indicated tolerance, or to a value that is 
    % higher than the "true" probability (M.Psmear).
    mergeNum   = 0.4*round(segLMax(j)*diagCellsNum/smear.DomainLength);
    smearCells = [];
    pObjIt     = 0;
    maxIts     = 10;
    itNum      = 0;
    sandCells  = 1:diagCellsNum;
    %expCount   = 0;
    while itNum < maxIts
        % 2.1 Updates and initial checks
        itNum        = itNum + 1;
        smearCellsT  = (pCalc(j)-pObjIt)*diagCellsNum;
        smearNum     = round(smearCellsT/claySegCellNum);
        
        if smearCellsT < 1
            smearCellsT = round(smearCellsT);
            if smearCellsT == 0
                if verbose == 1
                    disp(['Tolerance of ' num2str(tolP) ...
                          ' cannot be met: not enough cell resolution.'])
                end
                break
            else
                smearCellsT = 1;
            end
        end
        
        if smearNum == 0
            smearNum = 1;
            claySegCellNum = ceil(smearCellsT);  % ceil to minimize its.
        end
        
        % 2.2 Check if all diag entries in smear window are already full of
        % smear. If not, loop over each smear segment, and place each 
        % segment at a random position within the longest sand segment that
        % exists at the end of the previous iteration.
%        if numel(smearCells) == diagCellsNum
%             expCount = expCount + 1;
%             assert(isempty(sandCells));
%             smearCells = 1:numel(smearCells)+1;
%             diagCellsNum = diagCellsNum +1;
%             warning(['diagonal cells in smear ' num2str(j) ' expanded ' ...
%                      num2str(expCount) ' times to tolerance-match P.'])   
%        
%        else
        if numel(smearCells) < diagCellsNum
            for k = 1:smearNum              
                
                % Cells of sand segment where smear segment is placed (the
                % longest).
                sandLim = find(diff(sandCells) > 1);
                sandBounds = [1 sandLim+1; sandLim numel(sandCells)];
                sandCellNum = diff(sandBounds)+1;
                [~, sandSegId] = max(sandCellNum);
                sandSegCells = sandCells(sandBounds(1,sandSegId)):...
                                         sandCells(sandBounds(2,sandSegId));
                sandSegCellNum = numel(sandSegCells);
                
                % Smear segment cells to add
                if claySegCellNum == 1 && ~isempty(smearCells) || ...
                   claySegCellNum < mergeNum
                   % Very short segments, so we attach them to existing
                   % smear segments
                    if numel(sandSegCells) <= claySegCellNum
                        possibleUpperCell = sandSegCells(1);
                    else
                        possibleUpperCell = sandSegCells(1):...
                                            (sandSegCells(end) - ...
                                             (claySegCellNum-1));
                    end
                    if min(possibleUpperCell) == 1
                        cellsToAdd = max(possibleUpperCell):...
                            (max(possibleUpperCell) + (claySegCellNum-1));
                    elseif (max(possibleUpperCell) + ...
                            (claySegCellNum-1)) == max(sandCells)
                        cellsToAdd = min(possibleUpperCell):...
                            (min(possibleUpperCell) + (claySegCellNum-1));
                    else
                        if rand(1) > 0.5
                            cellsToAdd = min(possibleUpperCell):...
                                (min(possibleUpperCell) + (claySegCellNum-1));
                        else
                            cellsToAdd = max(possibleUpperCell):...
                                (max(possibleUpperCell)+(claySegCellNum-1));
                        end
                    end
                    
                else    % Randomly placed, each cell has equal P(smear)
                    initCell = sandSegCells(randi(sandSegCellNum, 1));
                    if (max(sandSegCells) - initCell) + 1 < claySegCellNum
                        cellsToAddBot = initCell:sandSegCells(end);
                        cellsToAddTop = sandSegCells(1):...
                            (sandSegCells(1) + claySegCellNum - ...
                             (numel(cellsToAddBot)+1));
                        cellsToAdd = [cellsToAddTop cellsToAddBot];
                    else
                        cellsToAdd = initCell:(initCell + (claySegCellNum-1));
                    end
                end
                
                % Updates
                smearCells = unique([cellsToAdd'; smearCells]);
                sandCells = 1:diagCellsNum;
                if max(sandCells) < max(smearCells)
                    smearCells(smearCells > max(sandCells)) = [];
                end
                sandCells(smearCells) = [];
            end
        end
        
        % 2.3 Vector values in each diagonal of each smear's domain. This 
        % is for a single diagonal, will be repmat to the correct number of
        % diagonals for each smear later.
        dvals             = false(diagCellsNum, 1);
        dvals(smearCells) = true;
        
        % 2.4 Divide segments of length > segLMax, which result from the
        % random placing of smears.
        if numel(smearCells) < diagCellsNum
            pathLength = smear.DomainLength*(diagCellsNum/G.cartDims(1));
            start1 = strfind([0, dvals'==1],[0 1]);
            end1 = strfind([dvals'==1, 0],[1 0]);
            segCells = end1 - start1 + 1;
            segLen = segCells*(pathLength/diagCellsNum);
            itMaxDiv = 5;
            itDiv = 1;
            while any(segLen > 1.05*segLMax(j))
                ids = segLen > segLMax(j);
                segCenters = start1(ids) + fix(segCells(ids)/2);
                dvals(segCenters) = false;
                start1 = strfind([0, dvals'==1],[0 1]);
                end1 = strfind([dvals'==1, 0],[1 0]);
                segCells = end1 - start1 + 1;
                segLen = segCells*(pathLength/diagCellsNum);
                if all(segLen <= 1.05*segLMax(j))
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
        
        % 2.5 Prepare for matrix, compute iteration probability and compare
        % with calculated probability
        %if -cDiagBound(2,j) > 0
        if cDiagBound(1,j) < 0 && ismember(M.unit(idc(j)), FS.HW.Id)
            addZerosBot = zeros(winBot(j)-1, 1);
            addZerosTop = zeros(G.cartDims(2)-numel([dvals; addZerosBot]), 1);
            dvals = repmat([addZerosTop; dvals; addZerosBot],1,cnDiag(j));
        elseif cDiagBound(2,j) > 0 && ismember(M.unit(idc(j)), FS.FW.Id)
            addZerosTop = zeros(G.cartDims(2)-winTop(j), 1);
            addZerosBot = zeros(G.cartDims(2)-numel([addZerosTop; dvals]), 1);
            dvals = repmat([addZerosTop; dvals; addZerosBot],1,cnDiag(j));
        elseif -cDiagBound(2,j) >= 0
            addZeros = zeros(G.cartDims(1)-numel(dvals), 1);
            dvals = repmat([addZeros; dvals],1,cnDiag(j));
        else
            addZeros = zeros(G.cartDims(1)-numel(dvals), 1);
            dvals = repmat([dvals; addZeros],1,cnDiag(j));
        end
        numLow = -cDiagBound(2,j); numUp = -cDiagBound(1,j);
        diagIds = numLow:numUp;
        Mtest = transpose(full(spdiags(dvals, diagIds, ...
                                       false(G.cartDims(1)))));
        %figure(9); spy(Mtest); pause           % visualize current smear.
        pObjIt = sum(sum(Mtest)) / sum(dvals_Mtest);
        
        % 2.6 If P already within tolerance of true value, we are done 
        % here.
        if abs(pCalc(j) - pObjIt) < tolP 
            break
   
        elseif pObjIt > pCalc(j)
            % If it is larger than true value, then enter second iterative
            % loop. Here we just remove cells one row at a time (iteration),
            % until P is within tolerance. Alternatively, if P becomes lower 
            % than true value again, it means that true P cannot be matched 
            % with this cell resolution. So, we break anyway.
            % Values are removed randomly (but maintaining the same choice
            % for a given segment) from either the top or the bottom of the
            % longest segment.
            if verbose == 1
                disp(['2nd iterative loop (removal) required to match ' ...
                      'Psmear in the 2D subdomain:'])
            end
            itNum2  = 0;
            maxIts2 = 50;
            while itNum2 < maxIts2
                % Prepare
                itNum2 = itNum2 + 1;
                
                % Find longest smear 
                idBreaks = find(diff(smearCells) > 1)';
                idStartEnd = [unique([1, idBreaks+1]); ...
                              unique([idBreaks, numel(smearCells)])];
                segCellNum = diff(idStartEnd)+1;
                [~, idLongest] = max(segCellNum);
                longSegCells = smearCells(idStartEnd(1,idLongest):...
                               idStartEnd(2,idLongest));
                
                % Remove and compute 2nd iterative loop probability
                if itNum2 < 2             % only once/smear for consistency
                    val = round(rand(1)); % either 0 (end) or 1 (first cell)
                end
                if val == 1
                    cellRemove = longSegCells(1);
                else
                    cellRemove = longSegCells(end);
                end
                smearCells(smearCells == cellRemove)= [];                  %#ok
                dvals2              = false(diagCellsNum,1);
                dvals2(smearCells)  = true;
                %if -cDiagBound(2,j) > 0
                if cDiagBound(1,j) < 0 && ismember(M.unit(idc(j)), FS.HW.Id)
                    addZerosTop = zeros(G.cartDims(2)-numel([dvals2; addZerosBot]), 1);
                    dvals2 = repmat([addZerosTop; dvals2; addZerosBot],1,cnDiag(j));
                elseif cDiagBound(2,j) > 0 && ismember(M.unit(idc(j)), FS.FW.Id)
                    addZerosBot = zeros(G.cartDims(2)-numel([addZerosTop; dvals2]), 1);
                    dvals2 = repmat([addZerosTop; dvals2; addZerosBot],1,cnDiag(j));
                elseif -cDiagBound(2,j) >= 0
                    addZeros = zeros(G.cartDims(1)-numel(dvals2), 1);
                    dvals2 = repmat([addZeros; dvals2],1,cnDiag(j));
                else
                    addZeros = zeros(G.cartDims(1)-numel(dvals2), 1);
                    dvals2 = repmat([dvals2; addZeros],1,cnDiag(j));
                end
                Mtest = transpose(full(spdiags(dvals2, diagIds, ...
                                               false(G.cartDims(1)))));
                pObjIt2 = sum(Mtest(Mtest==true)) / sum(dvals_Mtest);
                
                % Break if needed and verbose messages
                if abs(pObjIt2 - pCalc(j)) < tolP
                    if verbose == 1
                        disp(['Tolerance of ' num2str(tolP) ' met in ' ...
                              num2str(itNum2) ' iterations.'])
                    end
                    break
                elseif itNum2 == maxIts2
                    if verbose == 1
                        disp(['Tolerance cannot be met: max number of ' ...
                              'iterations reached'])
                    end
                elseif pObjIt2 < pCalc(j)
                    if verbose == 1
                        disp(['Tolerance of ' num2str(tolP) ...
                              ' cannot be met: not enough cell resolution.'])
                        disp(['Number of iterations = ' num2str(itNum2)])
                    end
                    break
                end
            end
            
            % pass dvals and object probabilities from 2nd iterative loop
            dvals = dvals2;
            pObjIt  = pObjIt2;
            break
        
        elseif itNum == maxIts && abs(pCalc(j) - pObjIt) > tolP 
            % less smear than required. Add one cell at a time, picking a
            % random cell with 0 that is ATTACHED to a smear (avoid placing
            % extra short smear segments). Shortest smear is chosen.
            if numel(smearCells) == diagCellsNum
                if verbose == 1
                    disp(['Psmear cannot be matched, not enough cell '...
                        'resolution'])
                end
            else
                if verbose == 1
                    disp(['2nd iterative loop (addition) required to match' ...
                        ' Psmear in the 2D subdomain:'])
                end
                itNum3  = 0;
                maxIts3 = 25;
                while itNum3 < maxIts3
                    itNum3  = itNum3 + 1;
                    dvals3   = false(diagCellsNum,1);
                    dvals3(smearCells) = true;
                    start1 = strfind([0, dvals3'==1],[0 1]);
                    end1 = strfind([dvals3'==1, 0],[1 0]);
                    segCells = end1 - start1 + 1;
                    [~, idMinL] = min(segCells);
                    addToStartOrEnd = rand(1);
                    if addToStartOrEnd >= 0.5 && end1(idMinL) < diagCellsNum || ...
                            addToStartOrEnd < 0.5 && start1(idMinL) == 1
                        idToAdd = end1(idMinL) + 1;
                    else
                        idToAdd = start1(idMinL) - 1;
                    end
                    assert(idToAdd > 0)
                    dvals3(idToAdd) = 1;
                    smearCells = find(dvals3);
                    
                    % Compute 3rd iterative loop probability
                    %if -cDiagBound(2,j) > 0
                    if cDiagBound(1,j) < 0 && ismember(M.unit(idc(j)), FS.HW.Id)
                        addZerosTop = zeros(G.cartDims(2)-numel([dvals3; addZerosBot]), 1);
                        dvals3 = repmat([addZerosTop; dvals3; addZerosBot],1,cnDiag(j));
                    elseif cDiagBound(2,j) > 0 && ismember(M.unit(idc(j)), FS.FW.Id)
                        addZerosBot = zeros(G.cartDims(2)-numel([addZerosTop; dvals3]), 1);
                        dvals3 = repmat([addZerosTop; dvals3; addZerosBot],1,cnDiag(j));
                    elseif -cDiagBound(2,j) >= 0
                        addZeros = zeros(G.cartDims(1)-numel(dvals3), 1);
                        dvals3 = repmat([addZeros; dvals3],1,cnDiag(j));
                    else
                        addZeros = zeros(G.cartDims(1)-numel(dvals3), 1);
                        dvals3 = repmat([dvals3; addZeros],1,cnDiag(j));
                    end
                    Mtest = transpose(full(spdiags(dvals3, diagIds, ...
                                                   false(G.cartDims(1)))));
                    pObjIt3 = sum(Mtest(Mtest==true)) / sum(dvals_Mtest);
                    
                    % Breaks and verbose messages
                    if abs(pObjIt3 - pCalc(j)) < tolP
                        if verbose == 1
                            disp(['Tolerance of ' num2str(tolP) ' met in ' ...
                                num2str(itNum3) ' iterations.'])
                        end
                        break
                    elseif itNum3 == maxIts3
                        if verbose == 1
                            disp(['Tolerance cannot be met: max number of ' ...
                                'iterations reached'])
                        end
                    elseif pObjIt3 > pCalc(j)
                        if verbose == 1
                            disp(['Tolerance of ' num2str(tolP) ...
                                ' cannot be met: not enough cell resolution.'])
                            disp(['Number of iterations = ' num2str(itNum3)])
                        end
                        break
                    end
                end
                
                % pass dvals and object probabilities from 3rd iterative loop
                dvals = dvals3;
                pObjIt  = pObjIt3;
                break
            end
        end
    end

    % 3. Final messages and assign values to mapping matrix
    if verbose == 1
        disp(['SSFc-derived P = ' num2str(pCalc(j)) '. '...
              'Final smear ' num2str(j) ' P = ' ...
              num2str(pObjIt) ' (' num2str(itNum) ' main iterations)'])
    end
    if pObjIt > 0
        if verbose == 1 
            numSmearSegFin = ceil(numel(find(diff(dvals(:,1))))*0.5);
            disp(['Input smear length [m] = ' num2str(segLMax(j)) ...
                  '. Number of final smear segments = ' ...
                  num2str(numSmearSegFin)]);
        end
        
        % Assign 0s and 1s to matrix that maps to grid.
        M.vals = full(spdiags(dvals, diagIds, M.vals));
    else
        if verbose == 1
            disp(['Input smear length [m] = ' num2str(segLMax(j)) ...
                '. Number of final smear segments = ' num2str(0)]);
        end
        dvals = false(diagCellsNum,1);
        numLow = -cDiagBound(2,j); numUp = -cDiagBound(1,j);
        addZeros = zeros(G.cartDims(1)-numel(dvals), 1);
        dvals = repmat([addZeros; dvals],1,cnDiag(j));
        diagIds = numLow:numUp;
        M.vals = full(spdiags(dvals, diagIds, M.vals));
    end
    if verbose == 1
        disp('___________________________________________________')
    end
    
    % Save final smear probability
    pObj(j) = pObjIt;
end

% To compromise spdiags truncation rules with smears chopped parallel
% to f.T (what we want) and not to f.L, the last step is to transpose. 
% Note that this needs to be done since the diagonal indices (diagIds) used 
% in spdiags during M.vals and Mtest assignment are flipped. Mtest is 
% directly transposed above, so that P is computed and matched with the 
% same exact configuration as in the definitive mapping matrix (M.vals).
M.vals = transpose(M.vals);

% Pass desired and final probabilities to output
pFin(M.Psmear<1) = pObj;
M.P = [M.Psmear; ...        % desired (calculated)
       pFin];               % obtained (including continuous smears)

end