function faultPoro = idealPackingPoro(G, M, FS, endMemberPoro)
% Compute based on Vcl and number of cells of each unit domain in fault,
% and get average Vcl based on total n of fault cells. Then use the ideal
% packing model of Revil et al. (2002) to get the porosity as a fcn of Vcl
% in fault, and end-member porosities.

% Get fault vcl
vcl = [FS.FW.Vcl FS.HW.Vcl];
if isfield(M, 'idSmearInRemoved')
    id = find(M.isclayIn);
    id_removed = id(M.idSmearInRemoved);
    vcl(id_removed) = [];
    endMemberPoro(id_removed, :) = [];
end
vcl_f = sum(sum(vcl(M.units)))/G.cells.num;

% Ideal packing model porosity
poroSand = mean(endMemberPoro(~M.isclay, 1));
idc = find(M.isclay);
for n=1:sum(M.isclay)
   poro(n) = endMemberPoro(idc(n), 1) + rand(1)*(endMemberPoro(idc(n), 2) - ...
                                                 endMemberPoro(idc(n), 1));
   x(n) = sum(ismember(M.vals, idc(n))) / sum(ismember(M.vals, idc));
end
poroClay = sum(poro.*x);

if sum(M.nDiag(~M.isclay)) > 0
    shalySand = vcl_f < poroSand;
    faultPoro(shalySand) = poroSand - vcl_f*(1-poroClay);
    faultPoro(~shalySand) = vcl_f*poroClay;
else
    faultPoro = poroClay;
end