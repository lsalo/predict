%%

clear
close all

% Stratigraphic inputs and fault
Nsim = 1000;
isClayVcl = 0.4;
vcl = [0.4 0.55 0.7 0.9];
T = [2 10 20 50 100];
zf = [100 1000 2000];                         

f.Dip = 60;
f.Throw = 100;    
f.Grid.TargetCellDim = [0.1, 1];       % needed for min smear dimensions

% Other dependent params
Tap = T ./ sind(f.Dip);
f.Disp = f.Throw/sind(f.Dip);

% Generate smear samples
N = numel(T)*numel(zf);
thick = repmat(repelem(T', 1, numel(vcl)), numel(zf), 1);
thickap = repmat(repelem(Tap', 1, numel(vcl)), numel(zf), 1);
zf_all = repelem(zf', numel(T), 1);
sThick        = zeros(Nsim, numel(vcl), N);
sThickInFault = zeros(Nsim, numel(vcl), N);
sLen          = zeros(Nsim, numel(vcl), N);
sSegLenMax    = zeros(Nsim, numel(vcl), N);
sDomLen       = zeros(Nsim, N);
sFrac         = zeros(Nsim, numel(vcl), N);
Trat          = zeros(Nsim, numel(vcl), N);
segLenFrac    = zeros(Nsim, numel(vcl), N);
for j=1:N
    [~, f.MatProps.SSFcBounds] = getSSFc(vcl, isClayVcl, zf_all(j), ...
                                         thick(j, :), f.Disp);
    for n=1:Nsim
        % Variable props
        f.MatProps.Thick = getFaultThickness(f.Disp, 1);
        [f.Alpha, f.Delta, f.Zeta] = getFaultAngles(f.Throw, f.Dip, ...
                                                    f.MatProps.Thick);
        f.MatProps.ResFric = getResidualFrictionAngle(vcl);
        f.MatProps.SSFc = getSSFc(vcl, isClayVcl, zf_all(j), thick(j, :), f.Disp);
        
        % smear struct
        smear = Smear(vcl, isClayVcl, thick(j, :), thickap(j, :), ...
                      zf_all(j), f, 1);
        
        % unpack
        sThick(n, :, j) = smear.Thick;
        Trat(n, :, j) = smear.Thick ./ thick(j, :);
        sThickInFault(n, :, j) = smear.ThickInFault;
        sLen(n, :, j) = smear.Length;
        sSegLenMax(n, :, j) = smear.SegLenMax;
        segLenFrac(n, :, j) = smear.SegLenMax ./ smear.Length;
        sDomLen(n, j) = smear.DomainLength;
        sFrac(n, :, j) = smear.Psmear;
    end
end

% Smear historgram params
nbins = 25;
edg_Trat = linspace(0, 0.3, nbins);
edg_sfrac = linspace(0, 1, nbins);
edg_segLenFrac = linspace(0, 1, nbins);
rat = f.Throw ./ T;
colrs = repmat(hsv(numel(T)), numel(vcl), 1);
id_vcl = repelem((1:numel(vcl))', numel(T), 1);
nextvcl = [1; 1 + find(diff(id_vcl))];  
zf_all = repmat(reshape(zf_all, numel(T), numel(zf)), numel(vcl), 1);
idd = repmat((1:numel(T))', numel(vcl), 1);

% Plots
latx = {'Interpreter', 'latex'};
lw = {'linewidth', 1};
sz = [14, 12];

fh1 = figure(1);
tiledlayout(1, numel(vcl), 'Padding', 'compact', 'TileSpacing', 'compact');
for j=1:numel(id_vcl)       %   DOES NOT CHANGE WITH ZF (ONLY T, Vcl)  
    if any(ismember(nextvcl, j))
        nexttile
        hold on
    end
    if j <= numel(rat) 
        histogram(Trat(:, id_vcl(j), idd(j)), edg_Trat, 'Normalization', ...
            'probability', 'DisplayStyle', 'stairs', 'DisplayName', ...
            ['$t / T =$ ' num2str(rat(j), '%1.1f')], 'EdgeColor', ...
            colrs(j, :));
        xlabel('$s_\mathrm{T} / T$ [-]', latx{:}, 'fontSize', sz(2))
        ylabel('P [-]', latx{:}, 'fontSize', sz(2))
        title(['$V_\mathrm{cl}$ = ' num2str(vcl(id_vcl(j)))], ...
            latx{:}, 'fontSize', sz(2))
        leg = legend(latx{:}, 'fontSize', sz(2), 'location', 'northeast');
        set(leg.BoxFace, 'ColorType','truecoloralpha', ...
            'ColorData', uint8(255*[1;1;1;.6])); 
    else
        histogram(Trat(:, id_vcl(j), idd(j)), edg_Trat, 'Normalization', ...
            'probability', 'DisplayStyle', 'stairs', ...
            'EdgeColor', colrs(j, :));
        title(['$V_\mathrm{cl}$ = ' num2str(vcl(id_vcl(j)))], ...
            latx{:}, 'fontSize', sz(2))
    end
    plot([0 0.3], [0 0], '-k', 'HandleVisibility','off')
    xlim([0, 0.3]); ylim([0 1]); grid on
    xticks([0 0.05 0.1 0.15 0.2 0.25 0.3])
    yticks([0 0.2 0.4 0.6 0.8 1]);
end
hold off
set(fh1, 'position', [500, 200, 250*numel(vcl), 225]);

lst = repmat({'-', '-', '-.', '--', ':'}', numel(vcl), 1);
lwi = repmat([0.5, .5, .5, .5, .5]', numel(vcl), 1);
fh2 = figure(2);
tiledlayout(numel(zf), numel(vcl), 'Padding', 'compact', 'TileSpacing', 'compact');
k = 0;
for n = 1:numel(zf)
    for j=1:numel(id_vcl)
        if any(ismember(nextvcl, j))
            nexttile
            hold on
        end
        if j <= numel(rat) && n == 1
            histogram(sFrac(:, id_vcl(j), idd(j)+numel(T)*k), edg_sfrac, 'Normalization', ...
                  'probability', 'DisplayStyle', 'stairs', 'DisplayName', ...
                  ['$t / T =$ ' num2str(rat(j), '%1.1f')], 'EdgeColor', ...
                  colrs(j, :), 'LineStyle', lst{j}, 'lineWidth', lwi(j));
            xlabel('$s_\mathrm{L} / f_\mathrm{L}$ [-]', latx{:}, 'fontSize', sz(2))
            ylabel('P [-]', latx{:}, 'fontSize', sz(2))
            title(['$V_\mathrm{cl}$ = ' num2str(vcl(id_vcl(j))) ...
                   ' $\vert$ $z_\mathrm{f}$ = ' num2str(zf_all(j, n)) ' m'], ...
                   latx{:}, 'fontSize', sz(2))
            leg = legend(latx{:}, 'fontSize', sz(2), 'location', 'northeast');
            set(leg.BoxFace, 'ColorType','truecoloralpha', ...
                'ColorData', uint8(255*[1;1;1;.6])); 
        else
            histogram(sFrac(:, id_vcl(j), idd(j)+numel(T)*k), edg_sfrac, 'Normalization', ...
                  'probability', 'DisplayStyle', 'stairs', ...
                  'EdgeColor', colrs(j, :), 'LineStyle', lst{j}, 'lineWidth', lwi(j));
            title(['$V_\mathrm{cl}$ = ' num2str(vcl(id_vcl(j))) ...
                   ' $\vert$ $z_\mathrm{f}$ = ' num2str(zf_all(j, n)) ' m'], ...
                   latx{:}, 'fontSize', sz(2))
        end
        plot([0 1], [0 0], '-k', 'HandleVisibility','off')
        xlim([0, 1]); ylim([0 1]); grid on
        xticks([0 0.2 0.4 0.6 0.8 1])
        yticks([0 0.2 0.4 0.6 0.8 1]);
    end
    k = k+1;
end
hold off
set(fh2, 'position', [500, 200, 250*numel(vcl), 225*numel(zf)]);

fh3 = figure(3);            % THIS RATIO JUST DEPENDS ON T, not zf or vcl
hold on
for n=1:numel(T)
    histogram(segLenFrac(:, 2, n), edg_segLenFrac, 'Normalization', ...
             'probability', 'DisplayStyle', 'stairs', 'DisplayName', ...
             ['$t / T =$ ' num2str(rat(n), '%1.1f')], 'EdgeColor', ...
             colrs(n, :), 'lineWidth', lwi(n));
    xlabel('$s_\mathrm{l} / s_\mathrm{L}$ [-]', latx{:}, 'fontSize', sz(2))
    ylabel('P [-]', latx{:}, 'fontSize', sz(2))
    legend(latx{:}, 'fontSize', sz(2), 'location', 'northeast')
end
plot([0 1], [0 0], '-k', 'HandleVisibility','off')
hold off
xlim([0, 1]); ylim([0 0.5]); grid on
xticks([0 0.2 0.4 0.6 0.8 1])
yticks([0 0.1 0.2 0.3 0.4 0.5]);
set(fh3, 'position', [500, 200, 250, 200]);

