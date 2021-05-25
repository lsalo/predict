function plotUpscaledPerm(faults)
%
%
%

%% Preliminaries

% Utilities
latx = {'Interpreter', 'latex'};
sz = [14, 12];

% Fault MatProps
perms = cell2mat(cellfun(@(x) x.Perm, faults, ...
                         'UniformOutput', false)) ./ (milli*darcy);
if any(any(perms < 0))
    id = unique([find(perms(:, 1)<0), find(perms(:, 2)<0), find(perms(:, 3)<0)]);
    warning(['Negative upscaled perms found in ' num2str(numel(id))...
             ' simulations (ignored).'])
    perms(id, :) = [];
end
    
% Hist params
K = log10(perms);
nbins = 25;
logMinP = min(min(K));
logMaxP = max(max(K));
edges = linspace(fix(logMinP)-1, fix(logMaxP)+1, nbins);


% Plot

% Histograms
fh = figure(randi(1000, 1, 1));
tiledlayout(3, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
labls = ["$\log_{10}(k_{xx}$ [mD])", ...
        "$\log_{10}(k_{yy}$ [mD])", ...
        "$\log_{10}(k_{zz}$ [mD])"];
nexttile(1)
histogram(K(:, 1), edges, 'Normalization', 'probability', ...
          'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 1)
ylabel(labls(1), latx{:}, 'fontSize', sz(2))
%ylabel('P [-]', latx{:}, 'fontSize', sz(2))
xlim([fix(logMinP)-1 fix(logMaxP)+1])
ylim([0 1]); yticks(0:.2:1)
grid on
%xticks(10.^(fix(logMinP)-1:2:fix(logMaxP)+1))

nexttile(5)
rr = [255, 125, 125]/255;
histogram(K(:, 2), edges, 'Normalization', 'probability', ...
          'FaceColor', rr, 'FaceAlpha', 1)
%xlabel('$\hat{k}_{yy}$ [mD]', latx{:}, 'fontSize', sz(2))
%ylabel('P [-]', latx{:}, 'fontSize', sz(2))
xlim([fix(logMinP)-1 fix(logMaxP)+1])
ylim([0 1]); yticks(0:.2:1)
grid on
%xticks(10.^(fix(logMinP)-1:2:fix(logMaxP)+1))

nexttile(9)
bb = [125, 125, 255]/255;
histogram(K(:, 3), edges, 'Normalization', 'probability', ...
          'FaceColor', bb, 'FaceAlpha', 1)
xlabel(labls(3), latx{:}, 'fontSize', sz(2))
%ylabel('P [-]', latx{:}, 'fontSize', sz(2))
xlim([fix(logMinP)-1 fix(logMaxP)+1])
ylim([0 1]); yticks(0:.2:1)
grid on
%xticks(10.^(fix(logMinP)-1:2:fix(logMaxP)+1))

% Scatters
tidss = [2 3 4 6 7 8];
tidsh = [1 5 9];
x = [2 3 1 3 1 2];
y = [1 1 2 2 3 3];
colrs = [128, 0, 0; 0, 119, 128; 128, 0, 0; 128, 119, 128; 0, 119, 128; ...
         128, 119, 128]./255;
[R, P] = corrcoef(K);           % corrcoeff and pval matrices
a = 0.05;                       % significance level
pvals = P'; pvals(tidsh) = [];
r = R'; r(tidsh) = [];
idlaby = [4 7];
lablsy = labls(2:end);
idlabx = [7 8];
lablsx = labls(1:end-1);
for n=1:numel(tidss)
    nexttile(tidss(n))
    colormap(hot);
    histogram2(K(:,x(n)), K(:,y(n)), edges, edges, ...
               'Normalization', 'Probability', 'DisplayStyle','tile', ...
               'ShowEmptyBins','off');
    c = colorbar;
    if n == 1
        c.Label.Interpreter = 'latex'; 
        c.Label.String = 'P [-]';
        c.Label.FontSize = 9;
    end
    %scatter(K(:,x(n)), K(:,y(n)), 4, 'o', 'MarkerEdgeColor', colrs(n, :))
     if pvals(n) < a
         colr = 'm';
         fw = 'bold';
     else
         colr = 'k';
         fw = 'normal';
     end
     text(edges(1) + 0.05*(edges(end) - edges(1)), ...
          edges(end) - 0.05*(edges(end) - edges(1)), ...
          num2str(round(r(n), 3)),  'color', colr, 'fontSize', 10, 'fontWeight', fw); 
    xlim([fix(logMinP)-1 fix(logMaxP)+1])
    ylim([fix(logMinP)-1 fix(logMaxP)+1])
    grid on
    if ismember(tidss(n), idlaby)
        idly = ismember(idlaby, tidss(n));
        ylabel(lablsy(idly), latx{:}, 'fontSize', sz(2))
    end
    if ismember(tidss(n), idlabx)
        idlx = ismember(idlabx, tidss(n));
        xlabel(lablsx(idlx), latx{:}, 'fontSize', sz(2))
    end
end
set(fh, 'position', [200, 200, 575, 400]);

end