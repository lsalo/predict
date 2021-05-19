%% Analysis of material properties (marginal distributions)
%
% 
% 

clear
close all

% Stratigraphic inputs
vcl = [0 0.15 0.25 0.35 0.4 0.55 0.7 0.9];
isClayVcl = 0.4;
thick = [2, 15, 50, 100];
zf = [50, 500, 1000, 2000, 3000];
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
f.thick = getFaultThickness();
thickVal = f.thick.fcn(repelem(f.disp, Nsim, 1));
DTratios = f.disp./thickVal;

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
phi = getResidualFrictionAngle(vcl);
phiv = zeros(Nsim, numel(vcl));
for n=1:numel(vcl)
    phiv(:, n) = phi.fcn{n}(Nsim);
end

nbins = 25;
edg_phi = linspace(0, 35, nbins);

fh2 = figure(2);
tiledlayout(1, numel(vcl), 'Padding', 'compact', 'TileSpacing', 'compact');
for n = 1:numel(vcl)
    nexttile
    histogram(phiv(:, n), edg_phi, 'Normalization', 'probability','FaceColor', ...
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
SSFcVals = zeros(Nsim, numel(vcl), N);
id = find(vcl >= isClayVcl);
zf_all = repelem(zf', numel(thick), 1);
thick_all = repmat(thick', numel(zf), numel(vcl));
fdisp = f.disp;
for j=1:N
    SSFc = getSSFc(vcl, isClayVcl, zf_all(j), thick_all(j, :), fdisp);
    SSFcVals(:, id, j) = cell2mat(cellfun(@(x) x(Nsim), SSFc.fcn(id), ...
                                          'uniformOutput', false));
end

% histogram params
nbins = 25;
edg_SSFc = linspace(0, 20, nbins);
colrs = repmat([0 0 0; 1 0 0; 0 0 1; 0.6 0.6 0.6], numel(zf), 1);
nextzf = 1:numel(thick):N;
tileids = 1:numel(id):N;
tileids = (tileids + (0:numel(id)-1)')';
rat = round(thick ./ f.disp, 2);

fh3 = figure(3);
tiledlayout(numel(zf), numel(id), 'Padding', 'compact', 'TileSpacing', 'compact');
for n = 1:numel(id)
    for j=1:N
        idt = find(ismember(nextzf, j));
         if any(idt)
             nexttile(tileids(idt, n))
             hold on
         end
        if j <= numel(thick) && n == 1
                histogram(SSFcVals(:, id(n), j), edg_SSFc, 'Normalization', ...
                  'probability', 'DisplayStyle', 'stairs', 'DisplayName', ...
                  ['$T / f.D =$ ' num2str(rat(j), '%1.2f')], 'EdgeColor', ...
                  colrs(j, :));
            xlabel('SSFc [-]', latx{:}, 'fontSize', sz(2))
            ylabel('P [-]', latx{:}, 'fontSize', sz(2))
            title(['$V_\mathrm{cl}$ = ' num2str(vcl(id(n))) ...
                   ' $\vert$ $z_\mathrm{f}$ = ' num2str(zf_all(j)) ' m'], ...
                   latx{:}, 'fontSize', sz(2))
            legend(latx{:}, 'fontSize', sz(2), 'location', 'southeast')
        else
            histogram(SSFcVals(:, id(n), j), edg_SSFc, 'Normalization', ...
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
set(fh3, 'position', [500, 200, 225*numel(id), 150*numel(zf)]);


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
    poro = getPorosity(vcls, isClayVcl, zmax_all(j, :), 'zmax', ...
                       zf_all(j), 0);
    poros(:, :, j) =  cell2mat(cellfun(@(x) x(Nsim), poro.fcn, ...
                                          'uniformOutput', false));
end
for j=1:numel(zmax)
    poro = getPorosity(vclc(1), isClayVcl, zmaxc(j), 'zmax', zf(1), 0);
    poroc(:, j) = cell2mat(cellfun(@(x) x(Nsim), poro.fcn, ...
                           'uniformOutput', false));
end

% histogram params
nbins = 25;
edg_poro = linspace(0, 0.6, nbins);
colrs = repmat([0 0 0; 1 0 0; 0 0 1; 0.6 0.6 0.6], numel(zf), 1);
colrs(idRem, :) = [];
N = numel(zf_all);
nextzf = [1; 1+find(diff(zf_all) > 0)];    
tileids = 1:numel(vcls):numel(zf)*numel(vcls);
tileids = (tileids + (0:numel(vcls)-1)')';

% Plots
fh4 = figure(4);
tiledlayout(numel(zf), numel(vcls), 'Padding', 'compact', 'TileSpacing', 'compact');
for n = 1:numel(vcls)
    for j=1:N
        idt = find(ismember(nextzf, j));
        if any(idt)
            nexttile(tileids(idt, n))
            hold on 
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
set(fh4, 'position', [500, 200, 225*numel(ids), 150*numel(zf)]);

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
g = [10 50 100 500 1000];           % shear strain: fault throw / fault thick.
N = numel(zf)*numel(clayMine)*numel(g);
g_all = repmat(g',N/numel(g), 1);
zf_all = repmat(repelem(zf', numel(g), 1), numel(clayMine), 1);
cm_all = repelem(clayMine', N/numel(clayMine), 1);
vclkr = [0.1 0.3 0.5 0.7 0.9];
krVals = zeros(Nsim, numel(vclkr), N);
for n=1:N
    poro = getPorosity(vclkr, isClayVcl, zf_all(n), 'zf', zf_all(n), 0);
    porov = cell2mat(cellfun(@(x) x(Nsim), poro.fcn, 'uniformOutput', false));
    krat = getAnisotropyRatio(vclkr, zf_all(n), cm_all{n});
    gv = repelem(g_all(n), Nsim, 1);
    for j=1:numel(vclkr)
        krVals(:, j, n) = krat.fcn{j}(gv,porov(:,j));
    end
end

% Histogram params
nbins = 25;
edg_krat = linspace(0, 10, nbins);
colrs = repmat([0 0 0; 1 0 0; 0 0 1; 0.5 0.5 0.5; 0 1 1], N/numel(g), 1);
styls = repmat({'-','-','--','-.',':'}', N/numel(g), 1);
lw = repmat([2 1.5 1 0.5 0.5]', N/numel(g), 1);
nextcm = 1:N/numel(clayMine):N; 
nextzf = [1; 1+find(diff(zf_all) ~= 0)];    
tileids = 1:numel(vclkr):numel(zf)*numel(vclkr);
tileids = (tileids + (0:numel(vclkr)-1)')';

% Plots
fh6 = figure(6);
tiledlayout(numel(zf), numel(vclkr), 'Padding', 'compact', 'TileSpacing', 'compact');
for n=1:numel(vclkr)
    %for j=1:N/numel(clayMine)                           % kaolinite
    %for j=N/numel(clayMine)+1:2*(N/numel(clayMine))     % mica
    for j = 2*(N/numel(clayMine))+1:3*(N/numel(clayMine)) % smectite
        idt = find(ismember(nextzf, j));
        if any(idt)
            if idt > size(tileids, 1) && idt < 2*size(tileids, 1)+1 
                idt = idt - size(tileids, 1);
            elseif idt >= 2*size(tileids, 1)+1
                idt = idt - 2*size(tileids, 1);
            end
            nexttile(tileids(idt, n))
            hold on 
        end
        if j <= numel(g) && n == 1
                histogram(krVals(:, n, j), edg_krat, 'Normalization', ...
                  'probability', 'DisplayStyle', 'stairs', 'DisplayName', ...
                  ['$\gamma =$ ' num2str(g_all(j))], 'EdgeColor', ...
                  colrs(j, :), 'lineStyle', styls{j}, 'lineWidth', lw(j));
            xlabel('$k^\prime$ [-]', latx{:}, 'fontSize', sz(2))
            ylabel('P [-]', latx{:}, 'fontSize', sz(2))
            title(['$V_\mathrm{cl}$ = ' num2str(vclkr(n)) ...
                   ' $\vert$ $z_\mathrm{f}$ = ' num2str(zf_all(j)) ' m'], ...
                   latx{:}, 'fontSize', sz(2))
            h = legend(latx{:}, 'fontSize', sz(2), 'location', 'northwest');
            set(h.BoxFace, 'ColorType','truecoloralpha', ...
                'ColorData', uint8(255*[1;1;1;.5]));  
        else
            histogram(krVals(:, n, j), edg_krat, 'Normalization', ...
                  'probability', 'DisplayStyle', 'stairs', ...
                  'EdgeColor', colrs(j, :), 'lineStyle', styls{j}, 'lineWidth', lw(j));
            title(['$V_\mathrm{cl}$ = ' num2str(vclkr(n)) ...
                   ' $\vert$ $z_\mathrm{f}$ = ' num2str(zf_all(j)) ' m'], ...
                   latx{:}, 'fontSize', sz(2))
        end
        plot([0 10], [0 0], '-k', 'HandleVisibility','off')
        xlim([0, 10]); xticks(0:2:10); 
        yticks([0 0.2 0.4 0.6 0.8 1]); ylim([0 1]); grid on;            
    end
end
hold off
set(fh6, 'position', [500, 200, 200*numel(vclkr), ...
    175*numel(zf)]);



%% Permeability

% Sands
vcls   = [0.1 0.25 0.35];
zmax_all = repmat(repmat(zmax', 1, numel(vcls)), numel(zf), 1);
zf_all   = repelem(zf', numel(zmax), 1);
idRem    = zf_all > zmax_all(:, 1);  % faulting depth must be <= max burial
zmax_all(idRem, :) = []; zf_all(idRem, :) = [];
N = numel(zf_all);
perms =  zeros(Nsim, numel(vcls), N);
for n=1:N
   perm = getPermeability(vcls, isClayVcl, zf_all(n), zmax_all(n, :), ...
                           []);
   perms(:, :, n) =cell2mat(cellfun(@(x) x(Nsim), perm.fcn', ...
                            'uniformOutput', false))/9.8692e-16;     % mD
end

% hist params
nbins = 50;
edg_perms = logspace(-4, 6, nbins);
colrs = repmat([0 0 0; 1 0 0; 0 0 1; 0.5 0.5 0.5], numel(zf), 1);
colrs(idRem, :) = [];
nextzf = [1; 1 + find(diff(zf_all))]; 
tileids = 1:numel(vcls):numel(zf)*numel(vcls);
tileids = (tileids + (0:numel(vcls)-1)')';

% Plot
fh7 = figure(7);
tiledlayout(numel(zf), numel(vcls), 'Padding', 'compact', 'TileSpacing', 'compact');
for n = 1:numel(vcls)
    for j=1:N
        idt = find(ismember(nextzf, j));
        if any(idt)
            nexttile(tileids(idt, n))
            hold on 
        end
        if j <= numel(zmax) && n == 1
                histogram(perms(:, n, j), edg_perms, 'Normalization', ...
                  'probability', 'DisplayStyle', 'stairs', 'DisplayName', ...
                  ['$z_\mathrm{max} =$ ' num2str(zmax_all(j, 1))], 'EdgeColor', ...
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
            histogram(perms(:, n, j), edg_perms, 'Normalization', ...
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
set(fh7, 'position', [500, 200, 200*numel(vcls), 175*numel(zf)]);


% Clays
vclc   = vcl(vcl >= isClayVcl);
zmaxc = repmat(zmax', 1, numel(vclc));
permc =  zeros(Nsim, numel(vclc), numel(zmax));
for n=1:numel(zmax)
    poro = getPorosity(vclc, isClayVcl, zmaxc(n, :), 'zmax', 1000, 0);
    porov = cell2mat(cellfun(@(x) x(Nsim), poro.fcn, 'uniformOutput', false));
    permf = getPermeability(vclc, isClayVcl, zf(1), zmaxc(n, :), []);
    for k = 1:numel(vclc)
        permc(:, k, n) = permf.fcn{k}(porov(:,k), poro.range{k}(1), ...
                                      poro.range{k}(2))/9.8692e-19;  %microD  
    end
end


% hist params
nbins = 50;
edg_permc = logspace(-5, 1, nbins);
colrs = [0 0 0; 1 0 0; 0 0 1; 0.5 0.5 0.5];  

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