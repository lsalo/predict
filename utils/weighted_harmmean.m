function H = weighted_harmmean(vals,weights)
%
% to compute weighted harmmean
%
% Example:
%   K = [10^-11, 4*10^-15];
%   w = [0.9643, 0.0357];
%   num = sum(w.*(1./K));
%   den = sum(w);
%   H = 1/(num/den);
num = sum(weights.*(1./vals));
den = sum(weights);
H = 1/(num/den);
end