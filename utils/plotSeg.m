function plotSeg(faults, nSeg)
% Plot strike-parallel segmentation (number of segments)
%
nSegv = cellfun(@(x) numel(x.SegLen), faults);
nbins = 10;
edg = linspace(nSeg.range(1), nSeg.range(2), nbins+1);
latx = {'Interpreter', 'latex'};
sz = [14, 12];

fh1 = figure(randi(10000, 1, 1));
histogram(nSegv, edg, 'Normalization', 'probability','FaceColor', [0.3 0.3 0.3])
xlabel('$N_\mathrm{seg}$ [-]', latx{:}, 'fontSize', sz(2))
ylabel('P [-]', latx{:}, 'fontSize', sz(2))
title(['$N_\mathrm{sim}$ = ' num2str(numel(faults))], ...
      latx{:}, 'fontSize', sz(1))
xlim(nSeg.range)
ylim([0 1])
grid on
xticks(nSeg.range(1):2:nSeg.range(2))
set(fh1, 'position', [500, 200, 275, 250]);


end