function plotPoroVclThick(faults, dim, fignum)
%
%
%

%% Preliminaries

% Utilities
latx = {'Interpreter', 'latex'};
sz = [14, 12];

% Fault MatProps
poro = cell2mat(cellfun(@(x) x.Poro, faults, 'UniformOutput', false));
if nargin < 2 || dim == 2
    thick = cell2mat(cellfun(@(x) x.MatProps.thick, faults, ...
                             'UniformOutput', false));
elseif dim == 3
    thick = cell2mat(cellfun(@(x) x.Thick, faults, 'UniformOutput', false));
end
disp =  cell2mat(cellfun(@(x) x.Disp, faults, 'UniformOutput', false));
vcl = cell2mat(cellfun(@(x) x.Vcl, faults, 'UniformOutput', false));

    
% Hist params
nbins = 25;
minV = 0;
maxV = 1;
edges = linspace(minV, maxV, nbins);
edg_T = logspace(1, 3, nbins);

% Plot
if nargin < 3
    fh1 = figure(randi(1000, 1, 1));
    tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
else
    fh1 = figure(fignum);
end

nexttile(1)
hold on
histogram(vcl , edges, 'Normalization', 'probability', ...
          'FaceColor', [0.3 0.3 0.3])
xlabel('$\mathrm{f}_{V_\mathrm{cl}}$ [-]', latx{:}, 'fontSize', sz(2))      
ylabel('P [-]', latx{:}, 'fontSize', sz(2))
xlim([0 1])
ylim([0 1])
grid on
xticks(0:.2:1)
yticks(0:.2:1)
hold off

nexttile(2)
hold on
histogram(disp./thick , edg_T, 'Normalization', 'probability','FaceColor', [0.3 0.3 0.3])
xlabel('$\mathrm{f}_\mathrm{D} / \mathrm{f}_\mathrm{T}$ [-]', latx{:}, 'fontSize', sz(2))
%ylabel('P [-]', latx{:}, 'fontSize', sz(2))
%title(['$N_\mathrm{sim} =$ ' num2str(numel(faults))], ...
%         latx{:}, 'fontSize', sz(1))
xlim([8 1100])
% ylim([0 0.2])
% yticks(0:.04:.2)
ylim([0 1])
yticks(0:.2:1)
grid on
set(gca,'XScale','log')
xticks([10 100 1000])
xticklabels({'10' '10^{2}' '10^{3}'})
%text(10, 0.9, ['$N_\mathrm{sim} =$ ' num2str(numel(faults))], ...
%     latx{:}, 'fontSize', sz(2))
hold off

nexttile(3)
hold on
histogram(poro , edges, 'Normalization', 'probability', ...
          'FaceColor', [0.3 0.3 0.3])
xlabel('$\mathrm{f}_\mathrm{n}$ [-]', latx{:}, 'fontSize', sz(2))
%ylabel('P [-]', latx{:}, 'fontSize', sz(2))
xlim([0 1])
ylim([0 1])
grid on
xticks(0:.2:1)
yticks(0:.2:1)
hold off
set(fh1, 'position', [100, 100, 600, 200]);
end