function M = faultMaterialMap(G, FS, faultDisp, faultThick, ...
                              smearThickInFault, Psmear)
%
% -----------------------------SUMMARY------------------------------------
% This function takes as inputs the Grid structure (G), footwall,
% hangingwall and fault information and creates a structure, M, which
% contains the information regarding the smears and sand distribution in
% the modeled domain (the fault). The mapping matrix itself, stored in
% M.vals and M.units, provides a direct map to the materials in the
% simulation grid.
%
% The matrix M.vals is initialized with all potential 1s (all cells in the 
% domain (the fault) which can potentially have smear, and 0s (all cells in
% the domain that will surely contain sand.
%
% Each smear (which occupies a given number of diagonals in M.vals as well
% as the simulation grid) is centered around the middle of the
% corresponding unit in the FW or HW. The number of diagonals are decided
% based on Ts (Ls).
%
% Note that MRST grid indexing (G) starts at bottom left, columns (x) move 
% faster. Standar MATLAB matrices start counting at top left, and rows (y) 
% move faster.
%
% Also note that drectly superposed sand layers (consecutive in FW or HW)
% should be collapsed in FW and HW variables. That is, FW, for example,
% should be passed as [1 0 1] and not [1 0 0 0 1] for eg thicknesses [20 60
% 20] and not [20 20 20 20 20]. Above the smearing threshold, consecutive
% clay layers can be passed (eg. FW.IsClay = [1 0 1 1 1]).
%
% ------------------------------INPUTS------------------------------------
% G   = MRST Grid structure
% fw  = footwall information structure.
% hw  = hangingwall information structure
% f   = fault information structure
% Ls  = length of smears in the direction parallel to a diagonal from the 
%       top left or top right of fault to bottom right or bottom left of 
%       fault.
% Tap = Apparent layer thickness along the fault [fw hw]
%
% ------------------------------OUTPUT------------------------------------
% M = matrix structure. For arrays, the FW goes first (left) and then the 
%     HW. Contains the following fields:
%   nDiagTot        = number of total diagonals in grid of size n*n, i.e. 
%                     2*n -1.
%   vals            = matrix of n*n that contains the mapping of 1s (smear)
%                     and 0s (sand) to the grid. Note that the smears are
%                     just initialized with all 1s here, and the final
%                     configuration is defined in another function, based 
%                     on the chosen geostatistical method.
%   units           = matrix of n*n that indicates which final unit is
%                     present in each cell
%   unit            = array with final unit, i.e. 1 (bottom of fw) to n
%                     (top of hw). Note that, in some cases, this number
%                     can be different from the input number of units given
%                     that if 2 or more sand units are one on top of
%                     another in the stratigraphy, they may be collapsed
%                     into 1 since it is assumed that they are equal and so
%                     generate the same type of sand-based FZ material.
%   isclay          = logical array indicating whether each final unit is
%                     sand (0) or clay (1).
%   nDiag           = number of diagonals that each subdomain has in the
%                     grid.
%   layerCenter     = center height of each layer in the FW and HW, 
%                     relative to the bottom of the modeled fault portion. 
%   layerTop        = top height of each layer in the FW and HW.
%   layerBot        = base height of each layer in the FW and HW.
%   layerDiagCenter = Grid diagonal number corresponding to the center of 
%                     each layer in the FW and HW. Follows the same
%                     convention as MATLAB's spdiags, i.e. 0 is the main
%                     diagonal, lower diagonals are negative and upper are
%                     positive. 
%   layerDiagTop    = Grid diagonal number corresponding to the top of 
%                     each layer in the FW and HW.
%   layerDiagBot    = Grid diagonal number corresponding to the base of 
%                     each layer in the FW and HW.
%   nDiagLayer      = Number of diagonals between layerDiagTop and
%                     layerDiagBot, i.e. the number of diagonals that each 
%                     layer in the FW and HW would take up.
%   DiagTop         = Grid diagonal number corresponding to the top of 
%                     each subdomain in the fault.
%   DiagBot         = Grid diagonal number corresponding to the bottom of 
%                     each subdomain in the FW and HW.
%   clayDiagBot     = Grid diagonal number corresponding to the bottom of 
%                     each clay source subdomain in the FW and HW.
%   unitIn          = array with input units, i.e. 1 (bottom of fw) to n
%                     (top of hw).
%   isclayIn        = logical array indicating whether each input unit is
%                     sand (0) or clay (1).
%   divLayerDiag    = Number of diagonals in the lower and upper triangle
%                     of the grid, when there is a layer whose 
%                     corresponding subdomain in the fault crosses the 
%                     main diagonal. Otherwise [0, 0];
%__________________________________________________________________________

% Initial values to Matrix structure
M.nDiagTot    = 2*G.cartDims(1) - 1;            % total number of diagonals
M.vals        = false(G.cartDims(1));           % actual matrix of 0s and 1s
M.units       = zeros(size(M.vals));            % Unit domain of each cell  (parent Id)
M.unit        = FS.ParentId;                    % Unit domain of each group (parent Id)
M.isclay      = [FS.FW.IsClay, FS.HW.IsClay];   % total units and clay or not  
M.unitIn      = M.unit;                         % For reference (unchanged)
M.isclayIn    = M.isclay;                       % "
        

% 1. Diagonals with potential smear
fDiagL = sqrt(faultDisp^2 + faultThick^2);
layerLs = zeros(1, max(FS.ParentId));
layerLs(M.isclay) = smearThickInFault;
M.nDiag = round((layerLs./fDiagL)*M.nDiagTot);

% fix number of diagonals if already out of bounds (smears "too thick")
c = find(M.isclay);
smearThickAsFault = 0;
if sum(M.nDiag) > sum(M.nDiagTot)
    smearThickAsFault = 1;
    M.nDiag = fix(M.nDiag ./ sum(M.nDiag).*M.nDiagTot);
    c1 = find(M.nDiag > 1);
    if any(M.nDiag == 1)
       M.nDiag(M.nDiag == 1) = 0;
    end
    toadd = M.nDiagTot - sum(M.nDiag);
    if toadd > 0
        toaddEach = fix(M.nDiag(c1).*toadd/sum(M.nDiag));
        M.nDiag(c1) = M.nDiag(c1) + toaddEach;
        idx = randi(numel(c1),1);
        if M.nDiagTot > sum(M.nDiag)
            M.nDiag(c1(idx)) = M.nDiag(c1(idx)) + ...
                               (M.nDiagTot - sum(M.nDiag));          
        end
        %M.nDiag(c1(idx)) = M.nDiag(c1(idx)) + toadd;
    end
end

% Position of stratigraphic layers with respect to diagonals in M
M.layerCenter = [cumsum(FS.Tap(FS.FW.Id))-FS.Tap(FS.FW.Id)/2 ...
                 cumsum(FS.Tap(FS.HW.Id))-FS.Tap(FS.HW.Id)/2];
M.layerTop    = [cumsum(FS.Tap(FS.FW.Id)) cumsum(FS.Tap(FS.HW.Id))];
M.layerBot    = [0 M.layerTop(FS.FW.Id(1:end-1)) ...
                 0 M.layerTop(FS.HW.Id(1:end-1))];
M.layerDiagCenter = round(M.layerCenter.*(G.cartDims(1)/faultDisp)) - ...
                    G.cartDims(1);
M.layerDiagCenter(FS.HW.Id) = M.layerDiagCenter(FS.HW.Id) + G.cartDims(1);
M.layerDiagTop = round(M.layerTop.*(G.cartDims(1)/faultDisp)) - ...
                 G.cartDims(1);
M.layerDiagTop(FS.HW.Id) = M.layerDiagTop(FS.HW.Id) + G.cartDims(1);
M.layerDiagTop(M.layerDiagTop>0) = M.layerDiagTop(M.layerDiagTop>0)-1;
M.layerDiagBot = round(M.layerBot.*(G.cartDims(1)/faultDisp)) - ...
                 G.cartDims(1);
M.layerDiagBot(FS.HW.Id) = M.layerDiagBot(FS.HW.Id) + G.cartDims(1)+1;
M.layerDiagBot(M.layerDiagBot==0) = 1;
M.layerDiagBot(1) = -(G.cartDims(1)-1);
M.nDiagLayer = (M.layerDiagTop - M.layerDiagBot)+1;

M.DiagTop = M.layerDiagCenter + fix((M.nDiag-1)/2); % fix to be conservative, and both fix for symmetry.
M.DiagBot = M.layerDiagCenter - fix((M.nDiag-1)/2);

% If any M.nDiag is 2 we consider 3 diagonals; if any M.nDiag is 1 we neglect it.
if any(all([M.nDiag(M.isclay) < 3; M.nDiag(M.isclay) > 1])) 
    nid = find(all([M.nDiag < 3; M.nDiag > 1]));
    M.DiagTop(nid) = M.layerDiagCenter(nid) + 1;
    M.DiagBot(nid) = M.layerDiagCenter(nid) - 1;
    M.nDiag(nid) = 3;
end

% Indicate if the fault material is all smear
% -> This block needed because in case of thick smears not near center, they
%    would likely have part outside of this fault throw interval, hence
%    this throw interval should not necessarily be all smear.
if smearThickAsFault == 1
    if ~any(M.DiagBot(c1) <= -(G.cartDims(1) - 1)) || ...
       ~any(M.DiagTop(c1) >= (G.cartDims(1) - 1))
        smearThickAsFault = 0;
    end
end

% Here we may need to adjust nDiag for top and bottom smears, since they were
% centered around the source layer and may no longer cover all diagonals and/or
% be out of bounds (if smearThickAsFault == 1).
if smearThickAsFault == 1
   M.DiagBot(c(1)) = -(G.cartDims(1) - 1);
   M.nDiag(c(1)) = M.DiagTop(c(1)) - M.DiagBot(c(1)) + 1;
   M.DiagTop(c(end)) = G.cartDims(1) - 1;
   M.nDiag(c(end)) = M.DiagTop(c(end)) - M.DiagBot(c(end)) + 1;
end

% Check and remove potential smear overlaps (we don't consider this here).
if sum(M.isclay)>1        
    for n=1:numel(c)-1
        if M.DiagTop(c(n)) >= M.DiagBot(c(n+1)) && ...
                M.DiagTop(c(n)) <= M.DiagTop(c(n+1))
            M.DiagTop(c(n)) = M.DiagBot(c(n+1))-1;
            
        elseif M.DiagTop(c(n)) >= M.DiagBot(c(n+1))
            M.DiagBot(c(n+1)) = M.DiagTop(c(n))+1;
            M.DiagTop(c(n+1)) = M.DiagTop(c(n))+1;
        end
    end
end
if any(M.DiagTop < M.DiagBot)
    itnum = 0;
    while (any(M.DiagTop < M.DiagBot))
        itnum = itnum + 1;
        idc = find(M.DiagTop(c) < M.DiagBot(c));
        if idc > 1
            M.DiagTop(c(idc-1)) = M.DiagTop(c(idc));
            M.DiagBot(c(idc)) = M.DiagTop(c(idc));
        else
            M.DiagTop(c(idc)) = M.DiagBot(c(idc));
        end
        if itnum > 50
           error('Maximum number of iterations surpassed.')
        end
    end
end

% Make sure no values out of bounds
M.DiagBot(M.DiagBot<-(G.cartDims(1)-1)) = -(G.cartDims(1)-1); 
M.DiagBot(M.DiagBot>(G.cartDims(1)-1)) = G.cartDims(1)-1;
M.DiagTop(M.DiagTop<-(G.cartDims(1)-1)) = -(G.cartDims(1)-1);
M.DiagTop(M.DiagTop>(G.cartDims(1)-1)) = G.cartDims(1)-1;
M.nDiag = abs(M.DiagTop - M.DiagBot) + 1;
if any((M.DiagTop - M.DiagBot) >= M.nDiagTot-1)
    idmx = M.nDiag > M.nDiagTot;
    if sum(idmx) == 1
        M.DiagBot(~idmx) = M.DiagTop(~idmx);
    elseif sum(idmx) > 1
        idmxc = idmx;
        idmxc(idmxc==0) = [];
        idmxc = min(c(idmxc));
        idmx2 = false(1,numel(idmx));
        idmx2(idmxc) = true;
        M.DiagBot(~idmx2) = M.DiagTop(~idmx2);
    end
    if M.DiagBot(idmx) < -(G.cartDims(1) - 1)
        M.DiagBot(idmx) = -(G.cartDims(1) - 1);
    end
    if M.DiagTop(idmx) > G.cartDims(1)
        M.DiagTop(idmx) = G.cartDims(1);
    end
end
M.nDiag(c) = abs(M.DiagBot(c)-M.DiagTop(c)) + 1;
M.nDiag(M.DiagBot == M.DiagTop) = 0;

% Make sure no layers with 1 diag at top and bottom
if M.DiagBot(c(1))+(G.cartDims(1)-1) == 1 
   M.DiagBot(c(1)) = -(G.cartDims(1)-1);
   M.nDiag(c(1)) = M.nDiag(c(1))+1;
end
if (G.cartDims(1)-1)-M.DiagTop(c(end)) == 1 
   M.DiagTop(c(end)) = (G.cartDims(1)-1);
   M.nDiag(c(end)) = M.nDiag(c(end))+1;
end
if any(M.DiagBot(c) == M.DiagTop(c))
   ideq =  M.DiagBot(c) == M.DiagTop(c);
   M.DiagBot(c(ideq)) = 0;
   M.DiagTop(c(ideq)) = 0;
end

% Make sure no single diagonals without smear are left
if smearThickAsFault == 1 && sum(M.nDiag) ~= M.nDiagTot
    isClayWithDiag = all([M.isclay; M.nDiag > 0]); 
    cWDiag = M.unit(isClayWithDiag);     % only those with nDiag > 0
    if M.DiagBot(cWDiag(1)) ~= -(G.cartDims(1) - 1)
        M.DiagBot(cWDiag(1)) = -(G.cartDims(1) - 1);
        M.nDiag(cWDiag(1)) = (M.DiagTop(cWDiag(1)) - M.DiagBot(cWDiag(1))) + 1;
    end
    if M.DiagTop(cWDiag(end)) ~= G.cartDims(1) - 1
        M.DiagTop(cWDiag(end)) = G.cartDims(1) - 1;
        M.nDiag(cWDiag(end)) = (M.DiagTop(cWDiag(end)) - M.DiagBot(cWDiag(end))) + 1;
    end
    for n=1:numel(cWDiag)-1
        gap = abs(M.DiagTop(cWDiag(n))-M.DiagBot(cWDiag(n+1)));
       if gap > 1 && M.nDiag(cWDiag(n)) ~= 0 && M.nDiag(cWDiag(n+1)) ~= 0
           M.DiagTop(cWDiag(n)) = M.DiagTop(cWDiag(n)) + (gap-1);
           M.nDiag(cWDiag(n)) = M.nDiag(cWDiag(n)) + (gap-1);
       end
    end
end

% Track which smear domains will be removed (if 0 diag)
M.Psmear = Psmear;
if any(M.nDiag(c) == 0)
    M.idSmearInRemoved = find(M.nDiag(c) == 0);
    M.Psmear(M.idSmearInRemoved) = [];
    idr = c(M.idSmearInRemoved);
    M.unit(idr) = [];
    M.isclay(idr) = [];
    M.nDiag(idr) = [];
    M.DiagBot(idr) = []; M.DiagTop(idr) = [];
end

% 2. Diagonals with sand
diagIds = -G.cartDims(1)+1:G.cartDims(1)-1;
idc = find(M.isclay);
vals = zeros(1,numel(diagIds));
%nOverlap = [0 0];
for n=1:sum(M.isclay)
   clayDiag = M.DiagBot(idc(n)):M.DiagTop(idc(n));
   if clayDiag(1) < 0
       vals = [zeros(1,(G.cartDims(1)-1)-abs(clayDiag(1))) clayDiag ...
               zeros(1,(G.cartDims(1)-1)-clayDiag(end))] + vals;
   elseif clayDiag(1) == 0 && numel(clayDiag) == 1
       vals = zeros(1, M.nDiagTot) + vals;
   elseif clayDiag(1) == 0
       vals = [zeros(1,G.cartDims(1)-1) clayDiag ...
               zeros(1,(G.cartDims(1)-1)-clayDiag(end))] + vals;
   else
       vals = [zeros(1,G.cartDims(1)) zeros(1,clayDiag(1)-1) clayDiag ...
               zeros(1,(G.cartDims(1)-1)-clayDiag(end))] + vals;
   end
end
idRem = abs(vals-diagIds) == 0;
diagIds(idRem) = [];
flag = 0;
if numel(diagIds)>1
        if ~any(all([M.DiagBot(M.isclay)<=0; M.DiagTop(M.isclay)>=0]))
            diagIds = [diagIds(diagIds<0) 0 diagIds(diagIds>0)];    
        end
        diffsf = diff(diagIds)>1;
        diffsi = [false diffsf(1:end-1)];
        if diffsf(end) == 1
            flag   = 1;
            valRem = diagIds(end);
            diffsf(end)  = [];
            diagIds(end) = [];
        end
        sandIds = [diagIds(1) diagIds(diffsi); diagIds(diffsf) diagIds(end)];
        M.DiagBot = M.DiagBot(M.isclay);
        M.DiagTop = M.DiagTop(M.isclay);
        M.clayDiagBot = M.DiagBot;
        M.DiagBot = sort([M.DiagBot sandIds(1,:)]);
        M.DiagTop = sort([M.DiagTop sandIds(2,:)]);
        [~, ids] = intersect(M.DiagBot, sandIds(1,:));
        idSandInFault = M.unit(~M.isclay);
        sandParent = zeros(1, numel(ids));
        if numel(ids) ~= numel(idSandInFault)
            for n=1:numel(ids)
                [~, idClosestSand] = min(abs(ids(n) - idSandInFault));
                sandParent(n) = idSandInFault(idClosestSand);
            end
        else
            sandParent = idSandInFault;
        end 
        nunits = sum([sum(M.isclay), numel(ids)]);
        unitAll = zeros(1, nunits);
        unitAll(ids) = sandParent;
        idclay = unitAll == 0;
        unitAll(idclay) = M.unit(M.isclay);
        M.unit = unitAll;
        M.isclay = false(1, numel(M.unit));
        M.isclay(idclay) = true;
        if any(M.DiagTop(M.DiagBot == 0) == 0)      % remove clay smear domain with 0 diags.
            %error('To code.')
            idTopZ = find(M.DiagTop == 0);
            idBotZ = find(M.DiagBot == 0);
            idTopBotZ = idTopZ(idTopZ == idBotZ);
            M.DiagTop(idTopBotZ) = [];
            M.DiagBot(idTopBotZ) = [];
            M.unit(idTopBotZ) = [];
            M.isclay(idTopBotZ) = [];
        end
        if any(M.DiagBot == 0) && any(M.DiagTop == 0) %&& ~all(M.DiagTop(M.DiagBot==0)==0)
            error('To code.')
            %M.DiagTop(M.DiagTop == 0) = [];
            %M.DiagBot(M.DiagBot == 0) = [];
        end
        M.nDiag = abs(M.DiagBot - M.DiagTop)+1;
        %M.unitIn = M.unit;
        %M.unit   = 1:numel(M.DiagBot);
        %M.isclayIn = M.isclay;
        %M.isclay = false(1,M.unit(end));
        %for n=1:numel(M.clayDiagBot)
        %    M.isclay(M.DiagBot==M.clayDiagBot(n)) = true;
        %end
end

if flag == 1
   if any(M.DiagBot == valRem+1)
       ij = M.DiagBot == valRem+1;
       M.DiagBot(ij) = valRem;
       M.nDiag(ij) = M.nDiag(ij) + 1;
   else
       ij = M.DiagTop == valRem+1;
       M.DiagTop(ij) = valRem;
       M.nDiag(ij) = M.nDiag(ij) + 1;
   end
end
M.nDiag(M.DiagBot == M.DiagTop) = 0;
M.DiagTop(M.nDiag==0) = 0;
M.DiagBot(M.nDiag==0) = 0;

% 3. Add field divLayerDiag for layer with lower and upper diags
if any(all([M.DiagBot<0; M.DiagTop>0]))
    idLay = all([M.DiagBot<0; M.DiagTop>0]);
    %M.DiagBot(idLay) = M.DiagBot(idLay)+1;
    M.divLayerDiag = [abs(M.DiagBot(idLay))+1 M.DiagTop(idLay)];
else
    %if (any(M.DiagTop == -1) && ~any(M.DiagBot==0)) || ...
    %   (any(M.DiagTop == -1) && ~isempty(M.DiagBot(M.DiagTop==0)==0) ...
    %   && sum(M.DiagBot(M.DiagTop==0)==0)==1 && M.DiagBot(M.DiagTop==0)==0)
    idbt = 1:numel(M.DiagBot);
    if sum(M.DiagBot(M.DiagTop == 0) == 0 > 0)
        idbt(all([M.DiagBot == 0; M.DiagTop == 0])) = [];
    end
    id = M.DiagTop == -1;
    if ~any(M.DiagTop(idbt) == 0) && ~any(M.DiagBot(idbt) == 0)
        if ~isempty(id)
            M.DiagTop(id) = 0;
        else
            id = M.DiagBot == 1;
            M.DiagBot(id) = 0;
        end
        M.nDiag(id) = M.nDiag(id)+1;
    end
    [v, id] = min(abs(M.DiagBot(idbt)));
    [v(2), id(2)] = min(abs(M.DiagTop(idbt)));
    [~, idLay] = min(v);
    idLay = idbt(id(idLay));
    if M.DiagBot(idLay) < 0
        M.divLayerDiag = [abs(M.DiagBot(idLay))+1 M.DiagTop(idLay)];
    elseif M.DiagBot(idLay) > 0
        M.divLayerDiag = [0 abs(M.DiagTop(idLay))];
    else
        M.divLayerDiag = [0 abs(M.DiagTop(idLay))+1];
    end
end

c = find(M.isclay);
if any(M.nDiag(c)==0)
   ij = M.nDiag(c) == 0;
   M.nDiag(c(ij)) = [];
   M.DiagTop(c(ij)) = [];
   M.DiagBot(c(ij)) = [];
   M.isclay(c(ij)) = [];
   M.unit(c(ij)) = [];
end

% In case of multiple, very thin smears we may need to adjust and add 1
% diagonal to some entries.
if sum(M.nDiag) ~= M.nDiagTot
   zerBotTop = find(all([M.DiagBot == 0; M.DiagTop == 0]));
   zerBotTop(zerBotTop == 1) = [];
   zerBotTop(zerBotTop == numel(M.isclay)) = [];
   assert(all(M.isclay(zerBotTop) == 0))
   diagPerLyr = (M.nDiagTot - sum(M.nDiag)) / numel(zerBotTop);
   if mod(diagPerLyr, 1) ~= 0
      nDiagMis = M.nDiagTot - sum(M.nDiag);
      diagPerLyr = 1;
      if nDiagMis < numel(zerBotTop)
          zerBotTop = zerBotTop(1:nDiagMis);
%      else
%           diagPerLyr = repelem(fix(nDiagMis / numel(zerBotTop)), numel(zerBotTop));
%           maxIts = 20;
%           nits = 0;
%           idLyrs = repmat(1:numel(diagPerLyr), 1, 50);
%           k = 1;
%           if sum(diagPerLyr) < nDiagMis
%               modifier = 1;
%           elseif sum(diagPerLyr) > nDiagMis
%               modifier = -1;
%           end
%           while sum(diagPerLyr) ~= nDiagMis && nits < maxIts
%               diagPerLyr(idLyrs(k)) = diagPerLyr(k) + modifier;
%               nits = nits + 1;
%               if nits == maxIts
%                   error('Loop did not converge within maxIts')
%               end
%           end
      end
   end
   M.nDiag(zerBotTop) = diagPerLyr;
   M.DiagBot(zerBotTop) = M.DiagTop(zerBotTop - 1) + 1;
   M.DiagTop(zerBotTop) = M.DiagBot(zerBotTop + 1) - 1;
end

% Remove any sand entries with 0 diagonals
if any(M.nDiag==0)
   id = M.nDiag == 0;
   M.nDiag(id) = [];
   M.DiagTop(id) = [];
   M.DiagBot(id) = [];
   M.isclay(id) = [];
   M.unit(id) = [];
end

% Adjust DiagBot and DiagTop if needed (very coarse grids, few nDiag)
if sum((abs(M.DiagBot - M.DiagTop) + 1) - M.nDiag) ~= 0
    id1 = -(G.cartDims(1)-1);
    M.DiagBot = [id1 id1+cumsum(M.nDiag(1:end-1))];
    M.DiagTop = [M.DiagBot(2:end)-1 G.cartDims(1)-1];
end

% Final check if we are missing a 0 in either DiagBot or DiagTop;
% [~, idb] = min(abs(M.DiagBot));
% [~, idt] = min(abs(M.DiagTop));
% if ~any([M.DiagBot(idb), M.DiagBot(idt), M.DiagTop(idt), ...
%          M.DiagTop(idb)] == 0) && 

% Check that we have the correct number of diagonals for each layer.
try
    assert(sum(M.nDiag) == M.nDiagTot)
catch
    error('Material mapping was unsuccessful.') 
end
%assert(sum(M.nDiag(1:fw.n-1))+M.divLayerDiag(1) == G.cartDims(1))
%assert(sum(M.nDiag(hw.id))+M.divLayerDiag(2) == G.cartDims(1)-1)

% 4. Populate mapping matrix with all potential 1s and sure 0s
for n = 1:numel(M.nDiag)
    if M.DiagTop(n) < M.DiagBot(n), M.DiagTop(n) = M.DiagBot(n); end    
    if M.isclay(n) == 1
        M.vals = full(spdiags(true(G.cartDims(1), M.nDiag(n)), ...
                              -M.DiagTop(n):-M.DiagBot(n), M.vals));   
    else
        M.vals = full(spdiags(false(G.cartDims(1), M.nDiag(n)), ...
                              -M.DiagTop(n):-M.DiagBot(n), M.vals));
    end   
    M.units = full(spdiags(M.unit(n)*ones(G.cartDims(1), M.nDiag(n)), ...
                           -M.DiagTop(n):-M.DiagBot(n), M.units));
end
M.units = transpose(M.units); 
% M.vals is also transposed later.

end