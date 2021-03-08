function plotUpscaledPerm(faults)
%
%
%

%% Preliminaries

% Utilities
latx = {'Interpreter', 'latex'};
sz = [16, 14];

% Fault MatProps
perms = cell2mat(cellfun(@(x) x.Perm, faults, ...
                         'UniformOutput', false)) ./ (milli*darcy);
                     
% Hist params
nbins = 25;
logMinP = log10(min(min(perms)));
logMaxP = log10(max(max(perms)));
edges = logspace(fix(logMinP)-1, fix(logMaxP)+1, nbins);


% Plot
fh = figure(randi(1000, 1, 1));
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile
histogram(perms(:, 1), edges, 'Normalization', 'probability', ...
          'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 1)
xlabel('$\hat{k}_{xx}$ [mD]', latx{:}, 'fontSize', sz(2))
ylabel('P [-]', latx{:}, 'fontSize', sz(2))
xlim([0.7*10^logMinP 2*10^logMaxP])
ylim([0 1])
grid on
set(gca,'XScale','log')
xticks(10.^(fix(logMinP)-1:2:fix(logMaxP)+1))
%xticklabels({'10^{-5}' '10^{-4}' '10^{-3}' '0.01' '0.1' '1' '10'})
text(2*10^logMinP, 0.9, ['$N_\mathrm{sim} = $ ' num2str(numel(faults))], ...
     latx{:}, 'fontSize', sz(2))

nexttile
histogram(perms(:, 2), edges, 'Normalization', 'probability', ...
          'FaceColor', [1 0.5 0.5], 'FaceAlpha', 1)
xlabel('$\hat{k}_{yy}$ [mD]', latx{:}, 'fontSize', sz(2))
%ylabel('P [-]', latx{:}, 'fontSize', sz(2))
xlim([0.7*10^logMinP 2*10^logMaxP])
ylim([0 1])
grid on
set(gca,'XScale','log')
xticks(10.^(fix(logMinP)-1:2:fix(logMaxP)+1))

nexttile
histogram(perms(:, 3), edges, 'Normalization', 'probability', ...
          'FaceColor', [0.5 0.5 1], 'FaceAlpha', 1)
xlabel('$\hat{k}_{zz}$ [mD]', latx{:}, 'fontSize', sz(2))
%ylabel('P [-]', latx{:}, 'fontSize', sz(2))
xlim([0.7*10^logMinP 2*10^logMaxP])
ylim([0 1])
grid on
set(gca,'XScale','log')
xticks(10.^(fix(logMinP)-1:2:fix(logMaxP)+1))
set(fh, 'position', [500, 200, 600, 250]);



end