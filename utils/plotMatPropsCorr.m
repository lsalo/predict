function plotMatPropsCorr(faults, FS)
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
                            'UniformOutput', false))/(milli*darcy);
end
thickRep = repmat(thick, numel(layerId), 1);
resFric = reshape(resFric, N*numel(layerId), 1);
SSFc = reshape(SSFc, N*numel(layerId), 1);
kprime = reshape(kprime, N*numel(layerId), 1);
perm = reshape(perm, N*numel(layerId), 1);
poro = reshape(poro, N*numel(layerId), 1);
layerIds = repelem(layerId', N);
clayId = find(FS.Vcl >= FS.IsClayVcl);
id = ismember(layerIds, clayId);

%% Plots
% Utilities
latx = {'Interpreter', 'latex'};
sz = [14, 12];

% resFric and SSFc
mrk = '.';
lst = {'-', '-.', ':', '--', '-', '-.'};
colr = [0 0 0; 1 0 0; 0 0 1; 0.5 0.5 0.5; 1 0 1; 0 1 1];

figure(36)
sh = scatterhist(resFric(id), SSFc(id), 'Group', layerIds(id), 'Marker', '.', ...
            'Kernel', 'on', 'location', 'southeast', 'Direction', 'out', ...
            'Color', colr, 'LineStyle', lst, 'LineWidth', 1.5, ...
            'Marker', mrk, 'MarkerSize', 4);
grid on
xlabel('$\phi_\mathrm{r}$ [deg.]', latx{:}, 'fontSize', sz(2))
ylabel('SSFc [-]', latx{:}, 'fontSize', sz(2))
sh(1).Legend.String = [['Vcl = ' num2str(FS.Vcl(clayId(1)))] ...
                       string(FS.Vcl(clayId(2:end)))];
sh(1).Legend.Interpreter = 'latex';
sh(1).Legend.FontSize = 10;
sh(1).Legend.Box = 'on';

% thickness and perm
figure(37)
sh = scatterhist(log10(faults{1}.Disp./thickRep), log10(perm), 'Group', layerIds, ...
            'Marker', '.', ...
            'Kernel', 'on', 'location', 'southeast', 'Direction', 'out', ...
            'Color', colr, 'LineStyle', lst, 'LineWidth', 1.5, ...
            'Marker', mrk, 'MarkerSize', 4);
grid on
xlabel('$\log_{10}(\frac{f_\mathrm{D}}{f_\mathrm{T}})$ [-]', latx{:}, 'fontSize', sz(2))
ylabel('$\log_{10}(k_{xx}$ [mD]$)$', latx{:}, 'fontSize', sz(2))
sh(1).Legend.String = [['Vcl = ' num2str(FS.Vcl(1))] ...
                       string(FS.Vcl(2:end))];
sh(1).Legend.Interpreter = 'latex';
sh(1).Legend.FontSize = 10;
sh(1).Legend.Box = 'on';
xlim([1 3])

% poro and perm

end