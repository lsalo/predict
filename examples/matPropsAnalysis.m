%% Analysis of material properties
%
% 
% 

clear
close all

% Stratigraphic inputs
vcl = [0 0.15 0.25 0.35 0.4 0.55 0.7 0.9];
isClayVcl = 0.4;
thick = [2, 15, 50, 100];
zf = [50, 500, 1000, 2000];
zmax = [500, 1000, 2000, 3000];
clayMine = {'kao', 'mic', 'sme'};

f.throw = 100;
f.dip = 60;
f.disp = f.throw/sind(f.dip);

Nsim = 1000;

% Plotting utils
latx = {'Interpreter', 'latex'};
lw = {'linewidth', 1};
sz = [14, 12];

%% Fault thickness

f.thick = getFaultThickness(f.disp, Nsim);
DTratios = f.disp./f.thick;

% Histogram parameters
nbins = 25;
edg_T = logspace(1, 3, nbins);

% Figure
fh1 = figure(1);
histogram(DTratios, edg_T, 'Normalization', 'probability','FaceColor', ...
          [0.3 0.3 0.3])
xlabel('$f_\mathrm{D} / f_\mathrm{T}$ [-]', latx{:}, 'fontSize', sz(2))
ylabel('P [-]', latx{:}, 'fontSize', sz(2))
title(['$N_\mathrm{sim} =$ ' num2str(Nsim)], latx{:}, 'fontSize', sz(1))
xlim([8 1100]); ylim([0 1]); grid on; set(gca,'XScale','log')
xticks([10 100 1000]); xticklabels({'10' '10^{2}' '10^{3}'})
set(fh1, 'position', [500, 200, 275, 250]);


%% Residual friction angle of clay and sand sources

phi = zeros(Nsim, numel(vcl));
for n=1:Nsim
    phi(n, :) = getResidualFrictionAngle(vcl);
end

nbins = 25;
edg_phi = linspace(0, 35, nbins);

fh2 = figure(2);
tiledlayout(1, numel(vcl), 'Padding', 'compact', 'TileSpacing', 'compact');
for n = 1:numel(vcl)
    nexttile
    histogram(phi(:, n), edg_phi, 'Normalization', 'probability','FaceColor', ...
              [0.3 0.3 0.3])
    if n == 1
        xlabel('$\phi_\mathrm{r}$ [deg.]', latx{:}, 'fontSize', sz(2))
        ylabel('P [-]', latx{:}, 'fontSize', sz(2))
        title(['$N_\mathrm{sim} =$ ' num2str(Nsim), ' $\vert$ $V_\mathrm{cl}$ = ' ...
               num2str(round(vcl(n), 2))], latx{:}, 'fontSize', sz(2))
    else
        title(['$V_\mathrm{cl}$ = ' num2str(round(vcl(n), 2))], latx{:}, ...
               'fontSize', sz(2))
    end
    xlim([0, 35]); ylim([0 1]); grid on; xticks([0 10 20 30]);
end
set(fh2, 'position', [500, 200, 150*numel(vcl), 175]);


%% SSFc 

N = numel(zf)*numel(thick);
SSFc = zeros(Nsim, numel(vcl), N);
id = find(vcl >= isClayVcl);
SSFcBounds = zeros(2, numel(id), N);
zf_all = repelem(zf', numel(thick), 1);
thick_all = repmat(thick', numel(zf), numel(vcl));
throw = f.throw;
parfor j=1:N
    [~, SSFcBounds(:, :, j)] = getSSFc(vcl, isClayVcl, zf_all(j), ...
                                       thick_all(j, :), throw);
    for n=1:Nsim
        SSFc(n, :, j) = getSSFc(vcl, isClayVcl, zf_all(j), ...
                                thick_all(j, :), throw);
    end
end

% histogram params
nbins = 25;
edg_SSFc = linspace(0, 20, nbins);
colrs = repmat(hsv(numel(thick)), numel(zf), 1);
nextzf = 1:numel(thick):N;
rat = round(f.throw ./ thick, 1);

fh3 = figure(3);
tiledlayout(numel(id), numel(zf), 'Padding', 'compact', 'TileSpacing', 'compact');
for n = 1:numel(id)
    for j=1:N
        if any(ismember(nextzf, j))
            nexttile
            hold on
            plot(repelem(SSFcBounds(1, n, j), 1, 2), [0 1], '-', 'color', ...
                [0.8 0.8 0.8], 'linewidth', 2, 'HandleVisibility','off')
            plot(repelem(SSFcBounds(2, n, j), 1, 2), [0 1], '-', 'color', ...
                [0.8 0.8 0.8], 'linewidth', 2, 'HandleVisibility','off') 
        end
        if j==1 && n == 1
            text(SSFcBounds(2, n, j)+0.5, 0.9, 'SSFc lim.', latx{:}, ...
                 'fontsize', sz(2))
        end
        if j <= numel(thick) && n == 1
                histogram(SSFc(:, id(n), j), edg_SSFc, 'Normalization', ...
                  'probability', 'DisplayStyle', 'stairs', 'DisplayName', ...
                  ['$t / T =$ ' num2str(rat(j), '%1.1f')], 'EdgeColor', ...
                  colrs(j, :));
            xlabel('SSFc [-]', latx{:}, 'fontSize', sz(2))
            ylabel('P [-]', latx{:}, 'fontSize', sz(2))
            title(['$V_\mathrm{cl}$ = ' num2str(vcl(id(n))) ...
                   ' $\vert$ $z_\mathrm{f}$ = ' num2str(zf_all(j)) ' m'], ...
                   latx{:}, 'fontSize', sz(2))
            legend(latx{:}, 'fontSize', sz(2), 'location', 'southeast')
        else
            histogram(SSFc(:, id(n), j), edg_SSFc, 'Normalization', ...
                  'probability', 'DisplayStyle', 'stairs', ...
                  'EdgeColor', colrs(j, :));
            title(['$V_\mathrm{cl}$ = ' num2str(vcl(id(n))) ...
                   ' $\vert$ $z_\mathrm{f}$ = ' num2str(zf_all(j)) ' m'], ...
                   latx{:}, 'fontSize', sz(2))
        end
        plot([0 24], [0 0], '-k', 'HandleVisibility','off')
        xlim([0, 20]); xticks([0 4 8 12 16 20]); yticks([0 0.2 0.4 0.6 0.8 1]);
        ylim([0 1]); grid on;
    end
end
hold off
set(fh3, 'position', [500, 200, 200*numel(zf), 150*numel(id)]);


%% Porosity

ids = find(vcl < isClayVcl);
vcls = vcl(ids);
idc = find(vcl >= isClayVcl);
vclc = vcl(idc);
zmax_all = repmat(repmat(zmax', 1, numel(vclc)), numel(zf), 1);
zf_all = repelem(zf', numel(zmax), 1);
idRem = zf_all > zmax_all(:, 1);  % faulting depth must be <= max burial
zmax_all(idRem, :) = []; zf_all(idRem, :) = [];
zmaxc = repmat(zmax', 1, 1);
poros = zeros(Nsim, numel(vcls), numel(zf_all));     % depends on zf, zmax, vcl
poroc = zeros(Nsim, numel(zmax));                    % only depends on zmax
for j=1:numel(zf_all)
    for n=1:Nsim
        poros_range = getPorosity(vcls, isClayVcl, zf_all(j), ...
                                  zmax_all(j, :), 'zmax', 0);
        poros(n, :, j) =  poros_range(:, 1);
    end
end
for j=1:numel(zmax)
    for n=1:Nsim
        poroc_range = getPorosity(vclc(1), isClayVcl, zf(1), zmaxc(j), ...
                                  'zmax', 0);
        poroc(n, j) =  poroc_range(1) + rand(1, 1) ...
                                      .* (poroc_range(2) - poroc_range(1));
    end
end

% histogram params
nbins = 25;
edg_poro = linspace(0, 0.6, nbins);
colrs = repmat(hsv(numel(zmax)), numel(zf), 1);
colrs(idRem, :) = [];
N = numel(zf_all);
nextzf = [1; 1+find(diff(zf_all) > 0)];        

% Plots
fh4 = figure(4);
tiledlayout(numel(vcls), numel(zf), 'Padding', 'compact', 'TileSpacing', 'compact');
for n = 1:numel(vcls)
    for j=1:N
        if any(ismember(nextzf, j))
            nexttile
            hold on
            %plot(repelem(SSFcBounds(1, n, j), 1, 2), [0 1], '-', 'color', ...
            %    [0.8 0.8 0.8], 'linewidth', 2, 'HandleVisibility','off')
            %plot(repelem(SSFcBounds(2, n, j), 1, 2), [0 1], '-', 'color', ...
            %    [0.8 0.8 0.8], 'linewidth', 2, 'HandleVisibility','off') 
        end
        if j <= numel(zmax) && n == 1
                histogram(poros(:, ids(n), j), edg_poro, 'Normalization', ...
                  'probability', 'DisplayStyle', 'stairs', 'DisplayName', ...
                  ['$z_\mathrm{max} =$ ' num2str(zmax_all(j))], 'EdgeColor', ...
                  colrs(j, :));
            xlabel('$n$ [-]', latx{:}, 'fontSize', sz(2))
            ylabel('P [-]', latx{:}, 'fontSize', sz(2))
            title(['$V_\mathrm{cl}$ = ' num2str(vcls(n)) ...
                   ' $\vert$ $z_\mathrm{f}$ = ' num2str(zf_all(j)) ' m'], ...
                   latx{:}, 'fontSize', sz(2))
            h = legend(latx{:}, 'fontSize', sz(2), 'location', 'northwest');
            set(h.BoxFace, 'ColorType','truecoloralpha', ...
                'ColorData', uint8(255*[1;1;1;.5]));  
        else
            histogram(poros(:, ids(n), j), edg_poro, 'Normalization', ...
                  'probability', 'DisplayStyle', 'stairs', ...
                  'EdgeColor', colrs(j, :));
            title(['$V_\mathrm{cl}$ = ' num2str(vcls(n)) ...
                   ' $\vert$ $z_\mathrm{f}$ = ' num2str(zf_all(j)) ' m'], ...
                   latx{:}, 'fontSize', sz(2))
        end
        plot([0 24], [0 0], '-k', 'HandleVisibility','off')
        xlim([0, 0.6]); xticks([0 0.1 0.2 0.3 0.4 0.5 0.6]); 
        yticks([0 0.2 0.4 0.6 0.8 1]); ylim([0 1]); grid on;
    end
end
hold off
set(fh4, 'position', [500, 200, 200*numel(zf), 150*numel(ids)]);

fh5 = figure(5);
hold on
for n=1:numel(zmax)
    histogram(poroc(:, n), edg_poro, 'Normalization', ...
              'probability', 'DisplayStyle', 'stairs', 'DisplayName', ...
              ['$z_\mathrm{max} =$ ' num2str(zmax(n))], 'EdgeColor', ...
              colrs(n, :));
end
plot([0 0.6], [0 0], '-k', 'HandleVisibility','off')
hold off
xlabel('$n$ [-]', latx{:}, 'fontSize', sz(2))
ylabel('P [-]', latx{:}, 'fontSize', sz(2))
%title('Clay Porosity', latx{:}, 'fontSize', sz(2))
legend(latx{:}, 'fontSize', sz(2), 'location', 'northeast')
xlim([0, 0.5]); xticks([0 0.1 0.2 0.3 0.4 0.5]); 
yticks([0 0.2 0.4 0.6 0.8 1]); ylim([0 1]); grid on;
set(fh5, 'position', [500, 200, 275, 225]);


%% Permeability Anisotropy ratio
idc    = find(vcl >= isClayVcl);
vclc   = vcl(idc);
zmax   = repelem(1500, 1, numel(idc));
porozf = zeros(numel(zf), 1);
for n=1:numel(zf)
    r = getPorosity(vclc(1), isClayVcl, zf(n), zmax, 'zf', 0);
    porozf(n) = r(1) + rand(1, 1) .* (r(2) - r(1));
end

shear_strain = [10 50 100 500 1000];           % fault throw / fault thick.
silt  = [0 1];
N = numel(zf)*numel(clayMine)*numel(shear_strain)*numel(silt);
cm_all     = repmat(repelem(clayMine', 1, 1), N/numel(clayMine), 1);
ss_all     = repmat(repelem(shear_strain', numel(clayMine), 1), ...
                    numel(zf)*numel(silt), 1);
zf_all     = repmat(repelem(zf', numel(clayMine)*numel(shear_strain), 1), ...
                    numel(silt), 1);
porozf_all = repmat(repelem(porozf, numel(clayMine)*numel(shear_strain), 1), ...
                    numel(silt), 1);
silt_all   = repelem(silt', N/numel(silt), 1);
krat = zeros(Nsim, N);
for j=1:N
    for n=1:Nsim
        krat(n, j) = getAnisotropyRatio(0.5, isClayVcl, zf_all(j), ...
                                        cm_all(j), ss_all(j), ...
                                        silt_all(j), porozf);
    end
end

% Histogram params
nbins = 50;
edg_krat = logspace(0, 4, nbins);
colrs = repmat(hsv(numel(clayMine)), N/numel(clayMine), 1);
nextcm = 1:numel(clayMine):N/2;  

% Plots
fh6 = figure(6);
tiledlayout(numel(zf), numel(shear_strain), 'Padding', 'compact', 'TileSpacing', 'compact');
idSilt = find(diff(silt_all) > 0);
for j=1:N/numel(silt)
    if any(ismember(nextcm, j))
        nexttile
        hold on
    end
    if j <= numel(clayMine)
        histogram(krat(:, j), edg_krat, 'Normalization', ...
            'probability', 'DisplayStyle', 'stairs', 'DisplayName', ...
            ['$m =$ ' cm_all{j}], 'EdgeColor', colrs(j, :));
        if j == numel(clayMine)
            histogram(krat(:, j+idSilt), edg_krat, 'Normalization', ...
                'probability', 'DisplayStyle', 'stairs', 'DisplayName', ...
            'silt = 0.1', 'EdgeColor', colrs(j, :), ...
            'LineStyle', '-.');
        else
            histogram(krat(:, j+idSilt), edg_krat, 'Normalization', ...
                'probability', 'DisplayStyle', 'stairs', 'EdgeColor', ...
                colrs(j, :), 'LineStyle', '-.', 'HandleVisibility','off');
        end
        xlabel('$k^\prime$ [-]', latx{:}, 'fontSize', sz(2))
        ylabel('P [-]', latx{:}, 'fontSize', sz(2))
        title(['  $z_\mathrm{f}$ = ' num2str(zf_all(j)) ...
               ' $\vert$ $\gamma = $' num2str(ss_all(j))], latx{:}, 'fontSize', sz(2))
        h = legend(latx{:}, 'fontSize', sz(2), 'location', 'northeast');
         set(h.BoxFace, 'ColorType','truecoloralpha', ...
                'ColorData', uint8(255*[1;1;1;.7])); 
    else
        histogram(krat(:, j), edg_krat, 'Normalization', ...
            'probability', 'DisplayStyle', 'stairs', ...
            'EdgeColor', colrs(j, :));
        histogram(krat(:, j+idSilt), edg_krat, 'Normalization', ...
            'probability', 'DisplayStyle', 'stairs', 'EdgeColor', ...
            colrs(j, :), 'LineStyle', '-.');
        title(['  $z_\mathrm{f}$ = ' num2str(zf_all(j)) ...
               ' $\vert$ $\gamma = $' num2str(ss_all(j))], latx{:}, 'fontSize', sz(2))
    end
    h = plot([1 10000], [0 0], '-k', 'HandleVisibility','off');
    set(gca,'XScale','log')
    xticks([1 10 100 1000 10000]); xticklabels({'1' '10' '10^{2}' '10^{3}' '10^{4}'})
    xlim([0, 10000]); yticks([0 0.2 0.4 0.6 0.8 1]); ylim([0 1]); grid on;
end
hold off
set(fh6, 'position', [500, 200, 200*numel(shear_strain), ...
    150*numel(zf)]);



%% Permeability

% Sands
vcls   = [0.1 0.25 0.35];
ids    = 1:numel(vcls);
zmax_all = repmat(repmat(zmax', 1, numel(vcls)), numel(zf), 1);
zf_all   = repelem(zf', numel(zmax), 1);
idRem    = zf_all > zmax_all(:, 1);  % faulting depth must be <= max burial
zmax_all(idRem, :) = []; zf_all(idRem, :) = [];
N = numel(zf_all);
perms =  zeros(Nsim, numel(vcls), N);
for n=1:N
   permf = getPermeability(vcls, isClayVcl, zf_all(n), zmax_all(n, :), ...
                           [], []);
   for k = 1:numel(vcls)
   perms(:, k, n) = permf{k}(Nsim) / (milli*darcy);     % mD
   end
end

% hist params
nbins = 25;
edg_perms = logspace(-4, 6, nbins);
colrs = repmat(hsv(numel(zmax)), numel(zf), 1);
colrs(idRem, :) = [];
nextzf = [1; 1 + find(diff(zf_all))];  

% Plot
fh7 = figure(7);
tiledlayout(numel(vcls), numel(zf), 'Padding', 'compact', 'TileSpacing', 'compact');
for n = 1:numel(vcls)
    for j=1:N
        if any(ismember(nextzf, j))
            nexttile
            hold on
        end
        if j <= numel(zmax) && n == 1
                histogram(perms(:, ids(n), j), edg_perms, 'Normalization', ...
                  'probability', 'DisplayStyle', 'stairs', 'DisplayName', ...
                  ['$z_\mathrm{max} =$ ' num2str(zmax_all(j))], 'EdgeColor', ...
                  colrs(j, :));
            xlabel('$k_\mathrm{xx}$ [mD]', latx{:}, 'fontSize', sz(2))
            ylabel('P [-]', latx{:}, 'fontSize', sz(2))
            title(['$V_\mathrm{cl}$ = ' num2str(vcls(n)) ...
                   ' $\vert$ $z_\mathrm{f}$ = ' num2str(zf_all(j)) ' m'], ...
                   latx{:}, 'fontSize', sz(2))
            h = legend(latx{:}, 'fontSize', sz(2), 'location', 'northwest');
            set(h.BoxFace, 'ColorType','truecoloralpha', ...
                'ColorData', uint8(255*[1;1;1;.7])); 
        else
            histogram(perms(:, ids(n), j), edg_perms, 'Normalization', ...
                  'probability', 'DisplayStyle', 'stairs', ...
                  'EdgeColor', colrs(j, :));
            title(['$V_\mathrm{cl}$ = ' num2str(vcls(n)) ...
                   ' $\vert$ $z_\mathrm{f}$ = ' num2str(zf_all(j)) ' m'], ...
                   latx{:}, 'fontSize', sz(2))
        end
        plot([0.0001 10^6], [0 0], '-k', 'HandleVisibility','off')
        set(gca,'XScale','log')
        xlim([0.0001, 10^6]); xticks([1e-4 0.01 1 10^2 10^4 10^6]); ...
        xticklabels({'10^{-4}' '10^{-2}' '1' '10^2' '10^{4}' '10^{6}'})
        yticks([0 0.2 0.4 0.6 0.8 1]); ylim([0 1]); grid on;
    end
end
hold off
set(fh7, 'position', [500, 200, 200*numel(zf), 150*numel(vcls)]);


% Clays
vclc   = vcl(vcl >= isClayVcl);
idc = 1:numel(vclc);
zmaxc = repmat(zmax', 1, numel(vclc));
permc =  zeros(Nsim, numel(vclc), numel(zmax));
for n=1:numel(zmax)
    poro = getPorosity(vclc, isClayVcl, zf(1), zmaxc(n, :), 'zmax', 0);
    permf = getPermeability(vclc, isClayVcl, zf(1), zmaxc(n, :), ...
                           [], poro);
    for k = 1:numel(vclc)
   permc(:, k, n) = permf{k}(Nsim) / (micro*darcy);    
   end
end


% hist params
nbins = 25;
edg_permc = logspace(-5, 1, nbins);
colrs = hsv(numel(zmax));
%nextzf = [1; 1 + find(diff(zf_all))];  

% Plot
fh8 = figure(8);
tiledlayout(1, numel(vclc), 'Padding', 'compact', 'TileSpacing', 'compact');
for n = 1:numel(vclc)
    nexttile
    hold on
    for j=1:numel(zmax)
        if n == 1
            histogram(permc(:, n, j), edg_permc, 'Normalization', ...
                  'probability', 'DisplayStyle', 'stairs', 'DisplayName', ...
                  ['$z_\mathrm{m} =$ ' num2str(zmax(j))], 'EdgeColor', ...
                  colrs(j, :));
            xlabel('$k_\mathrm{xx}$ [$\mu$D]', latx{:}, 'fontSize', sz(2))
            ylabel('P [-]', latx{:}, 'fontSize', sz(2))
            title(['$V_\mathrm{cl}$ = ' num2str(vclc(n))], ...
                   latx{:}, 'fontSize', sz(2))
            legend(latx{:}, 'fontSize', sz(2), 'location', 'northwest')
        else
            histogram(permc(:, n, j), edg_permc, 'Normalization', ...
                  'probability', 'DisplayStyle', 'stairs', ...
                  'EdgeColor', colrs(j, :));
            title(['$V_\mathrm{cl}$ = ' num2str(vclc(n))], ...
                   latx{:}, 'fontSize', sz(2))
        end
        plot([10^-5 10], [0 0], '-k', 'HandleVisibility','off')
        set(gca,'XScale','log')
        xlim([10^-5, 10]); xticks([1e-5 1e-4 0.001 0.01 0.1 1 10]); ...
        xticklabels({'10^{-5}' '10^{-4}' '10^{-3}' '10^{-2}' '0.1' '1' '10'})
        yticks([0 0.2 0.4 0.6 0.8 1]); ylim([0 1]); grid on;
    end
end
hold off
set(fh8, 'position', [500, 200, 200*numel(vclc), 200]);