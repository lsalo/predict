function SGR = getSGR(thick, vsh, res)
%
% Get SGR for a given throw interval in a faulted section.
%
% INPUTS
%   thick: true thickness of each layer in the FW and HW. 1x2 cell.
%   vsh: shale fraction of each layer in the FW and HW. 1x2 cell.
%   res: resolution, value relative to thickness unit in thick. For
%        example, if thick is in m, res = 1 means 1 value per m, 0.01 means
%        1 value per cm. Or, if thick is in cm, res = 1 means 1 value per 
%        cm.
% 
% OUTPUT:
% SGR: Yielding's (1997) shale gouge ratio, at the indicated resolution
% (e.g. one value per m). Values are given as a fraction, not percent.
%

if nargin < 3 
    res = 1;
end

% Get number of outputs and heights
throw = sum(thick{1});                  % throw 
n     = throw/res;                      % n SGR values
h     = linspace(0,throw,n)';           % height at which SGR is computed

% Get footwall thicknesses for SGR computation
hfw   = repmat(cumsum(thick{1}), n, 1) - h;
Tfw   = repmat(thick{1}, n, 1);
Tfw(hfw<0) = 0;
Tfw = min(abs(hfw), Tfw);    

% Get HW thicknesses for SGR computation
nlay = numel(thick{2});
hw_thick = repmat(throw - sum(Tfw, 2), 1, nlay);
thw_layer = repmat(thick{2}, n , 1);
hw_thick(:, 2:end) = hw_thick(:, 2:end) - cumsum(thw_layer(:, 1:nlay-1), 2);
hw_thick(hw_thick < 0) = 0;
Thw = min(thw_layer, hw_thick); 

% calculate
vsh = repmat(cell2mat(vsh), n, 1);
T = [Tfw, Thw];
assert(all(sum(T, 2) == throw));
SGR = sum(vsh.*T, 2)/throw;


end