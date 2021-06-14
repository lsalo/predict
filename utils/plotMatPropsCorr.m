function [R, P] = plotMatPropsCorr(faults, FS, idPlot)
%
%
%

%% Preparations
% Fault MatProps
layerId = FS.ParentId;
thick = cell2mat(cellfun(@(x) x.MatProps.thick, faults, ...
                 'UniformOutput', false));
%DTratios = faults{1}.Disp ./ thick;
N = numel(thick);
resFric = zeros(N, numel(layerId));
SSFc = zeros(N, numel(layerId));
kprime = zeros(N, numel(layerId));
perm = zeros(N, numel(layerId));
poro = zeros(N, numel(layerId));
for n = 1:numel(layerId)
    resFric(:, n) = cell2mat(cellfun(@(x) x.MatProps.resFric(n), faults, ...
                         'UniformOutput', false));
    SSFc(:, n) = cell2mat(cellfun(@(x) x.MatProps.ssfc(n), faults, ...
                    'UniformOutput', false));
    kprime(:, n) = cell2mat(cellfun(@(x) x.MatProps.permAnisoRatio(n), faults, ...
                            'UniformOutput', false));
    perm(:, n) = cell2mat(cellfun(@(x) x.MatProps.perm(n), faults, ...
                            'UniformOutput', false))/(milli*darcy);
    poro(:, n) = cell2mat(cellfun(@(x) x.MatProps.poro(n), faults, ...
                            'UniformOutput', false));
end
%thickRep = repmat(thick, numel(layerId), 1);
%resFricResh = reshape(resFric, N*numel(layerId), 1);
%SSFcResh = reshape(SSFc, N*numel(layerId), 1);
%kprime = reshape(kprime, N*numel(layerId), 1);
%permResh = reshape(perm, N*numel(layerId), 1);
%poro = reshape(poro, N*numel(layerId), 1);
%layerIds = repelem(layerId', N);
%clayId = find(FS.Vcl >= FS.IsClayVcl);
%id = ismember(layerIds, clayId);

%% Scatterhist Plots
% Utilities
latx = {'Interpreter', 'latex'};
sz = [14, 12];

% resFric and SSFc
% mrk = '.';
% lst = {'-', '-.', ':', '--', '-', '-.'};
% colr = [0 0 0; 1 0 0; 0 0 1; 0.5 0.5 0.5; 1 0 1; 0 1 1];

% figure(36)
% sh = scatterhist(resFricResh(id), SSFcResh(id), 'Group', layerIds(id), 'Marker', '.', ...
%             'Kernel', 'on', 'location', 'southeast', 'Direction', 'out', ...
%             'Color', colr, 'LineStyle', lst, 'LineWidth', 1.5, ...
%             'Marker', mrk, 'MarkerSize', 4);
% grid on
% xlabel('$\phi_\mathrm{r}$ [deg.]', latx{:}, 'fontSize', sz(2))
% ylabel('SSFc [-]', latx{:}, 'fontSize', sz(2))
% sh(1).Legend.String = [['Vcl = ' num2str(FS.Vcl(clayId(1)))] ...
%                        string(FS.Vcl(clayId(2:end)))];
% sh(1).Legend.Interpreter = 'latex';
% sh(1).Legend.FontSize = 10;
% sh(1).Legend.Box = 'on';
% 
% % thickness and perm
% figure(37)
% sh = scatterhist(log10(faults{1}.Disp./thickRep), log10(permResh), 'Group', layerIds, ...
%             'Marker', '.', ...
%             'Kernel', 'on', 'location', 'southeast', 'Direction', 'out', ...
%             'Color', colr, 'LineStyle', lst, 'LineWidth', 1.5, ...
%             'Marker', mrk, 'MarkerSize', 4);
% grid on
% xlabel('$\log_{10}(\frac{f_\mathrm{D}}{f_\mathrm{T}})$ [-]', latx{:}, 'fontSize', sz(2))
% ylabel('$\log_{10}(k_{xx}$ [mD]$)$', latx{:}, 'fontSize', sz(2))
% sh(1).Legend.String = [['Vcl = ' num2str(FS.Vcl(1))] ...
%                        string(FS.Vcl(2:end))];
% sh(1).Legend.Interpreter = 'latex';
% sh(1).Legend.FontSize = 10;
% sh(1).Legend.Box = 'on';
% xlim([1 3])


%% Correlation plot
% Delete data if constant
varn = ["t","phi","ssfc","n","k","k'"];
idn = 4;
nbins = [20, 20, 20, 20, 20, 20];
if FS.Vcl(idPlot) < FS.IsClayVcl
    varn(3) = [];   % remove ssfc
    idn = 3;
    nbins(3) = [];
    T = [log10(thick), resFric(:, idPlot), poro(:,idPlot), ...
         log10(perm(:,idPlot)), kprime(:,idPlot)];
    idk = 4;
    if max(T(:,4)) == min(T(:,4))
        T(:, idk) = [];
        varn(4) = [];
        nbins(4) = [];
    end
else
    T = [log10(thick), resFric(:, idPlot), SSFc(:,idPlot), poro(:,idPlot), ...
         log10(perm(:,idPlot)), kprime(:,idPlot)];
end

nv = size(T, 2);
edges = cell(1, nv);
for n=1:nv
    % Edges for histograms
    m = min(T(:,n));
    M = max(T(:,n));
    edges{n} = linspace(m, M, nbins(n));
end

% Correlation analysis
[R, P] = corrcoef(T);           % corrcoeff and pval matrices
a = 0.01;                       % significance level

% Plot
fh = figure(38);
t = tiledlayout(nv, nv, 'Padding', 'compact', 'TileSpacing', 'compact');
title(t,['Layer id: ' num2str(idPlot) ...
         ' $\mid V_\mathrm{cl}=$ ' num2str(FS.Vcl(idPlot)) ....
         ' $\mid z_\mathrm{f} =$ ' num2str(unique(FS.DepthFaulting(1))) ...
         ' m $\mid z_\mathrm{max} =$ ' num2str(unique(FS.DepthBurial(1))) ...
         ' m'], latx{:}, 'fontSize', sz(1));
% histograms
tidsh = 1:nv+1:nv^2;
labls = ["$\log_{10}(\mathrm{f}_\mathrm{T}$[m])", ...
         "$\phi_\mathrm{r}$[deg.]", "SSFc[-]", "$n$[-]", ...
         "$\log_{10}(k$[mD])", "$k^\prime$"];
if ~ismember('ssfc', varn), labls(3) = []; end
if ~ismember('k', varn), labls(4) = []; end
for n=1:numel(tidsh)
    nexttile(tidsh(n))
    histogram(T(:, n), edges{n}, 'Normalization', 'probability', ...
              'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 1)
    if n == 1
        ylabel(labls{n}, latx{:}, 'fontSize', sz(2))
    elseif n==numel(tidsh)
        xlabel(labls{n}, latx{:}, 'fontSize', sz(2))
    end
    xlim([min(edges{n}) max(edges{n})])
    p = histcounts(T(:,n), edges{n}, 'normalization', 'probability');
    ylim([0 round(max(p)+0.05, 1)])
    grid on
end

% scatter plots
tidss = 1:nv^2;
tidss(ismember(tidss,tidsh)) = [];
x = repmat(1:numel(varn), 1, nv); x(tidsh) = [];
y = repelem(1:numel(varn), 1, nv); y(tidsh) = [];
pvals = P'; pvals(tidsh) = [];
r = R'; r(tidsh) = [];
idlaby = numel(varn)+1:nv:nv^2;
lablsy = labls(2:end);
idlabx = nv^2-numel(varn)+1:nv^2-1;
lablsx = labls(1:end-1);
for n=1:numel(tidss)
    nexttile(tidss(n))
    scatter(T(:,x(n)), T(:,y(n)), '.', 'MarkerEdgeColor', [0.7 0.7 0.7])
    if pvals(n) < a
        colr = 'm';
        fw = 'bold';
    else
        colr = 'k';
        fw = 'normal';
    end
    text(edges{x(n)}(1) + 0.05*(edges{x(n)}(end) - edges{x(n)}(1)), ...
         edges{y(n)}(end) - 0.05*(edges{y(n)}(end) - edges{y(n)}(1)), ...
         num2str(round(r(n), 3)),  'color', colr, 'fontSize', 10, 'fontWeight', fw); 
    xlim([min(edges{x(n)}) max(edges{x(n)})])
    ylim([min(edges{y(n)}) max(edges{y(n)})])
    grid on
    if ismember(tidss(n), idlaby)
        idly = ismember(idlaby, tidss(n));
        ylabel(lablsy{idly}, latx{:}, 'fontSize', sz(2))
    end
    if ismember(tidss(n), idlabx)
        idlx = ismember(idlabx, tidss(n));
        xlabel(lablsx{idlx}, latx{:}, 'fontSize', sz(2))
    end
end
set(fh, 'position', [200, 0, 900, 600]);

end