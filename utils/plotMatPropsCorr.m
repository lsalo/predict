function plotMatPropsCorr(faults, FS)
%
%
%

%% Preparations
% Fault MatProps
layerId = FS.ParentId;
thick = cell2mat(cellfun(@(x) x.MatProps.Thick, faults, ...
                 'UniformOutput', false));
%DTratios = faults{1}.Disp ./ thick;
N = numel(thick);
phi = zeros(N, numel(layerId));
SSFc = zeros(N, numel(layerId));
kprime = zeros(N, numel(layerId));
perm = zeros(N, numel(layerId));
poro = zeros(N, numel(layerId));
for n = 1:numel(layerId)
    phi(:, n) = cell2mat(cellfun(@(x) x.MatProps.ResFric(n), faults, ...
                         'UniformOutput', false));
    SSFc(:, n) = cell2mat(cellfun(@(x) x.MatProps.SSFc(n), faults, ...
                    'UniformOutput', false));
    kprime(:, n) = cell2mat(cellfun(@(x) x.MatProps.PermAnisoRatio(n), faults, ...
                            'UniformOutput', false));
    perm(:, n) = faults{1}.MatProps.Perm{n}(N)/(milli*darcy);
    poror = faults{1}.MatProps.Poro(n, :);
    if diff(poror) == 0
        poro(:, n) = poror(1)*ones(N, 1);
    else
        poro(:, n) = poror(1) + rand(N, 1) .* abs(diff(poror));
    end
end
phi = reshape(phi, N*numel(layerId), 1);
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

% Scatterhist
mrk = '.';
lst = {'-', '-.', ':', '--', '-', '-.'};
colr = [0 0 0; 1 0 0; 0 0 1; 0.5 0.5 0.5; 1 0 1; 0 1 1];
figure(36)
sh = scatterhist(phi(id), SSFc(id), 'Group', layerIds(id), 'Marker', '.', ...
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
% Histogram2

end