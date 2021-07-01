function out = getAveragingPerm(vals, opt)
%
% vals = m x n array of values, where m is the permeability at each
% location of a given section (e.g. perm in each cell) and n is the number
% of sections (e.g. number of realizations).
%
% reference: Lie (2019) chapter 15.3 (MRST book)
%
out = zeros(1, numel(opt));
for n=1:numel(opt)
   if strcmp(opt{n}, 'ha')                  % harmonic-arithmetic avging
       out(n) = mean(harmmean(vals, 2));
   elseif strcmp(opt{n}, 'ah')              % arithmetic-harmonic avging
       out(n) = harmmean(mean(vals, 1));
   else
       error('method not supported')
   end
end


end
