function H = weighted_harmmean(vals, weights)
%
% compute weighted harmmean
%
% Example:
%   K = [10^-11, 4*10^-15];
%   w = [0.9643, 0.0357];
%   out = weighted_harmmean(K, w);
%
num = sum(weights.*(1./vals));
den = sum(weights);
H = 1/(num/den);
end