function plotMatPropsHist(faults, smears, FS, id)
%
%
%


%% Preparations
% Utilities
latx = {'Interpreter', 'latex'};
sz = [14, 12];

% Fault MatProps
thick = cell2mat(cellfun(@(x) x.MatProps.thick, faults, ...
                         'UniformOutput', false));
DTratios = faults{1}.Disp ./ thick;
phi = cell2mat(cellfun(@(x) x.MatProps.resFric(id), faults, ...
                         'UniformOutput', false));
SSFc = cell2mat(cellfun(@(x) x.MatProps.ssfc(id), faults, ...
                         'UniformOutput', false));
kprime = cell2mat(cellfun(@(x) x.MatProps.permAnisoRatio(id), faults, ...
                         'UniformOutput', false));
perm = cell2mat(cellfun(@(x) x.MatProps.perm(id), faults, ...
                         'UniformOutput', false))./(milli*darcy);
poro = cell2mat(cellfun(@(x) x.MatProps.poro(id), faults, ...
                         'UniformOutput', false));
%N = 1000;
%perm = faults{1}.MatProps.perm{id}(N)/(milli*darcy);
%poro = faults{1}.MatProps.poro(id, :);
% if diff(poro) == 0
%     poro = poro(1)*ones(N, 1);
% else
%     poro = poro(1) + rand(N, 1) .* abs(diff(poro));
% end

% Fault Histogram params
nbins = 25;
edg_T = logspace(1, 3, nbins);
edg_phi = linspace(0, 35, nbins);
edg_poro = linspace(0, 0.6, nbins);
lim_perm = [fix(log10(min(perm)))-1, fix(log10(max(perm)))+1];
edg_k = logspace(lim_perm(1), lim_perm(2), nbins);
edg_SSFc = linspace(0, 12, nbins);
edg_kprime = linspace(0, 10, nbins);

% Smear props
if faults{1}.MatMap.isclayIn(id) == 1
    sid = find(find(faults{1}.MatMap.isclayIn) == id);
    if id <= max(FS.FW.Id)
       layerThick = FS.FW.Thickness(id); 
    else
       layerThick = FS.HW.Thickness(FS.HW.Id == id);
    end
    sthick = cell2mat(cellfun(@(x) x.Thick(sid), smears, 'UniformOutput', false));
    Trat = sthick ./ layerThick;
    sfrac = cell2mat(cellfun(@(x) x.Psmear(sid), smears, 'UniformOutput', false));
    segLenMax = cell2mat(cellfun(@(x) x.SegLenMax(sid), smears, 'UniformOutput', false));
    slength = cell2mat(cellfun(@(x) x.Length(sid), smears, 'UniformOutput', false));
    segLenFrac = segLenMax ./ slength;

    % Smear historgram params
    edg_Trat = linspace(0, 0.3, nbins);
    edg_sfrac = linspace(0, 1, nbins);
    edg_segLenFrac = linspace(0, 1, nbins);
end

% Section info
if id <= max(FS.FW.Id)
    layer.vcl = FS.FW.Vcl(id);
    layer.thick = FS.FW.Thickness(id);
    layer.zf = FS.FW.DepthFaulting;
    layer.zmax = FS.FW.DepthBurial(id);
    layer.clayMine = FS.FW.ClayMine;
else
    idhw = find(FS.HW.Id == id);
    layer.vcl = FS.HW.Vcl(idhw);
    layer.thick = FS.HW.Thickness(idhw);
    layer.zf = FS.HW.DepthFaulting;
    layer.zmax = FS.HW.DepthBurial(idhw);
    layer.clayMine = FS.HW.ClayMine;
end
    

%% Plot Matprops
fh1 = figure(randi(1000, 1, 1));
histogram(DTratios, edg_T, 'Normalization', 'probability','FaceColor', [0.3 0.3 0.3])
xlabel('$f_\mathrm{D} / f_\mathrm{T}$ [-]', latx{:}, 'fontSize', sz(2))
ylabel('P [-]', latx{:}, 'fontSize', sz(2))
title(['$N_\mathrm{sim} =$ ' num2str(numel(faults))], ...
         latx{:}, 'fontSize', sz(1))
xlim([8 1100])
ylim([0 1])
grid on
set(gca,'XScale','log')
xticks([10 100 1000])
xticklabels({'10' '10^{2}' '10^{3}'})
%text(10, 0.9, ['$N_\mathrm{sim} =$ ' num2str(numel(faults))], ...
%     latx{:}, 'fontSize', sz(2))
set(fh1, 'position', [500, 200, 275, 250]);

fh2 = figure(randi(1000, 1, 1));
if faults{1}.MatMap.isclayIn(id)
    tiledlayout(1, 5, 'Padding', 'none', 'TileSpacing', 'compact');
else
    tiledlayout(1, 3, 'Padding', 'none', 'TileSpacing', 'compact');
end
nexttile
sgtitle(['$N_\mathrm{sim} =$ ' num2str(numel(faults)), ...
         ' $\vert$ Id = ' num2str(id) ' $\vert$ $V_\mathrm{cl}$ = ' num2str(layer.vcl) ...
         ' $\vert$ $T$ = ' num2str(layer.thick) ' m $\vert$ $z_\mathrm{f}$ = ' ...
         num2str(layer.zf) ' m $\vert$ $z_\mathrm{max}$ = ' num2str(layer.zmax) ...
         ' m $\vert$ $m$ = ' num2str(layer.clayMine)], latx{:}, 'fontSize', sz(2))
histogram(phi, edg_phi, 'Normalization', 'probability','FaceColor', [0.3 0.3 0.3])
xlabel('$\phi_\mathrm{r}$ [deg.]', latx{:}, 'fontSize', sz(2))
ylabel('P [-]', latx{:}, 'fontSize', sz(2))
%title('Residual friction', 'fontsize', sz(2))
xlim([0, 35])
ylim([0 1])
grid on
xticks([0 10 20 30])
%text(3, 0.9, ['Id = ' num2str(id)], latx{:}, 'fontSize', sz(2))

nexttile
histogram(poro, edg_poro, 'Normalization', 'probability','FaceColor', [0.3 0.3 0.3])
xlabel('$n$ [-]', latx{:}, 'fontSize', sz(2))
%ylabel('P [-]', latx{:}, 'fontSize', sz(2))
%title('Porosity', 'fontsize', sz(2))
xlim([0, 0.6])
ylim([0 1])
grid on
xticks([0 0.1 0.2 0.3 0.4 0.5 0.6])
%text(0.05, 0.9, ['N = ' num2str(N)], latx{:}, 'fontSize', sz(2))

nexttile
histogram(perm, edg_k, 'Normalization', 'probability','FaceColor', [0.3 0.3 0.3])
xlabel('$\hat{k}_{xx}$ [mD]', latx{:}, 'fontSize', sz(2))
%ylabel('P [-]', latx{:}, 'fontSize', sz(2))
%title('Permeability', 'fontsize', sz(2))
xlim([10^(lim_perm(1)) 10^(lim_perm(2))])
ylim([0 1])
grid on
set(gca,'XScale','log')
xticks(10.^(lim_perm(1):lim_perm(2)))
%xticklabels({'10^{-7}' '10^{-5}' '10^{-3}' '0.1' '10' '10^3'})
%text(1.3*10^(lim_perm(1)), 0.9,['N = ' num2str(N)],  latx{:}, 'fontSize', sz(2))

if faults{1}.MatMap.isclayIn(id)
    nexttile
    histogram(kprime, edg_kprime, 'Normalization', 'probability','FaceColor', [0.3 0.3 0.3])
    xlabel('$k^\prime$ [-]', latx{:}, 'fontSize', sz(2))
    %ylabel('P [-]', latx{:}, 'fontSize', sz(2))
    %title('Perm. anisotropy', 'fontsize', sz(2))
    xlim([0, 10])
    ylim([0 1])
    grid on
    xticks([0 2 4 6 8 10])
    
    nexttile
    histogram(SSFc, edg_SSFc, 'Normalization', 'probability','FaceColor', [0.3 0.3 0.3])
    xlabel('SSF$_\mathrm{c}$ [-]', latx{:}, 'fontSize', sz(2))
    %ylabel('P [-]', latx{:}, 'fontSize', sz(2))
    %title('Critical SSF', 'fontsize', sz(2))
    xlim([0, 12])
    ylim([0 1])
    grid on
    xticks(0:2:12)
end
%if faults{1}.MatMap.isclayIn(id)
%    set(fh2, 'position', [500, 200, 1100, 250]);
%else
set(fh2, 'position', [500, 200, 660, 250]);
%end


%% Plot Smears
if faults{1}.MatMap.isclayIn(id) == 1
    fh3 = figure(randi(1000, 1, 1));
    tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile
    sgtitle(['$N_\mathrm{sim} =$ ' num2str(numel(faults)), ...
         ' $\vert$ Id = ' num2str(id) ' $\vert$ $V_\mathrm{cl}$ = ' num2str(layer.vcl) ...
         ' $\vert$ $T$ = ' num2str(layer.thick) ' m $\vert$ $z_\mathrm{f}$ = ' ...
         num2str(layer.zf) ' m $\vert$ $z_\mathrm{max}$ = ' num2str(layer.zmax) ...
         ' m $\vert$ $m$ = ' num2str(layer.clayMine)], latx{:}, 'fontSize', sz(2))
    histogram(Trat, edg_Trat, 'Normalization', 'probability','FaceColor', [0.3 0.3 0.3])
    xlabel('$s_\mathrm{T} / T$ [-]', latx{:}, 'fontSize', sz(2))
    ylabel('P [-]', latx{:}, 'fontSize', sz(2))
    %title('Residual friction', 'fontsize', sz(2))
    xlim([0, 0.3])
    ylim([0 1])
    grid on
    xticks([0 0.05 0.1 0.15 0.2 0.25 0.3])
    
    nexttile
    histogram(sfrac, edg_sfrac, 'Normalization', 'probability','FaceColor', [0.3 0.3 0.3])
    xlabel('$s_\mathrm{L} / f_\mathrm{L}$ [-]', latx{:}, 'fontSize', sz(2))
    %ylabel('P [-]', latx{:}, 'fontSize', sz(2))
    %title('Residual friction', 'fontsize', sz(2))
    xlim([0, 1])
    ylim([0 1])
    grid on
    xticks(0:0.2:1)
    
    nexttile
    histogram(segLenFrac, edg_segLenFrac, 'Normalization', ...
              'probability','FaceColor', [0.3 0.3 0.3])
    xlabel('$s_\mathrm{l} / s_\mathrm{L}$ [-]', latx{:}, 'fontSize', sz(2))
    %ylabel('P [-]', latx{:}, 'fontSize', sz(2))
    %title('Residual friction', 'fontsize', sz(2))
    xlim([0, 1])
    ylim([0 1])
    grid on
    xticks(0:0.2:1)
    set(fh3, 'position', [500, 200, 660, 250]);
end

end