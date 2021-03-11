%%

clear
close all

% Stratigraphic inputs and fault
Nsim = 1000;
isClayVcl = 0.4;
vcl = [0.4 0.55 0.7 0.9];
T = [2 10 20 50 100];
zf = 1000;                         

f.Dip = 60;
f.Throw = 100;    
f.Grid.Resolution = [0.1, 1];       % needed for min smear dimensions

% Other dependent params
Tap = T ./ sind(f.Dip);
f.Disp = f.Throw/sind(f.Dip);

% Generate smear samples
N = numel(T);
thick = repelem(T', 1, numel(vcl));
thickap = repelem(Tap', 1, numel(vcl));
sThick        = zeros(Nsim, numel(vcl), N);
sThickInFault = zeros(Nsim, numel(vcl), N);
sLen          = zeros(Nsim, numel(vcl), N);
sSegLenMax    = zeros(Nsim, numel(vcl), N);
sDomLen       = zeros(Nsim, N);
sFrac        = zeros(Nsim, numel(vcl), N);
Trat          = zeros(Nsim, numel(vcl), N);
segLenFrac    = zeros(Nsim, numel(vcl), N);
for j=1:N
    [~, f.MatProps.SSFcBounds] = getSSFc(vcl, isClayVcl, zf, ...
                                         thick(j, :), f.Throw);
    for n=1:Nsim
        % Variable props
        f.MatProps.Thick = getFaultThickness(f.Disp, 1);
        [f.Alpha, f.Delta, f.Zeta] = getFaultAngles(f.Throw, f.Dip, ...
                                                    f.MatProps.Thick);
        f.MatProps.ResFric = getResidualFrictionAngle(vcl);
        f.MatProps.SSFc = getSSFc(vcl, isClayVcl, zf, thick(j, :), f.Throw);
        
        % smear struct
        smear = Smear(vcl, isClayVcl, thick(j, :), thickap(j, :), ...
                      zf, f, 1);
        
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
colrs = hsv(numel(T));


% Plots
latx = {'Interpreter', 'latex'};
lw = {'linewidth', 1};
sz = [14, 12];

fh1 = figure(1);
tiledlayout(1, numel(vcl), 'Padding', 'compact', 'TileSpacing', 'compact');
for n = 1:numel(vcl)
    nexttile
    hold on
    for j=1:N
        if n == 1
            histogram(Trat(:, n, j), edg_Trat, 'Normalization', ...
                  'probability', 'DisplayStyle', 'stairs', 'DisplayName', ...
                  ['$t / T =$ ' num2str(rat(j), '%1.1f')], 'EdgeColor', ...
                  colrs(j, :));
            xlabel('$s_\mathrm{T} / T$ [-]', latx{:}, 'fontSize', sz(2))
            ylabel('P [-]', latx{:}, 'fontSize', sz(2))
            title(['$V_\mathrm{cl}$ = ' num2str(vcl(n)) ...
                   ' $\vert$ $z_\mathrm{f}$ = ' num2str(zf) ' m'], ...
                   latx{:}, 'fontSize', sz(2))
            legend(latx{:}, 'fontSize', sz(2), 'location', 'northeast')
        else
            histogram(Trat(:, n, j), edg_Trat, 'Normalization', ...
                  'probability', 'DisplayStyle', 'stairs', ...
                  'EdgeColor', colrs(j, :));
            title(['$V_\mathrm{cl}$ = ' num2str(vcl(n)) ...
                   ' $\vert$ $z_\mathrm{f}$ = ' num2str(zf) ' m'], ...
                   latx{:}, 'fontSize', sz(2))
        end
        plot([0 24], [0 0], '-k', 'HandleVisibility','off')
        xlim([0, 0.3]); ylim([0 1]); grid on
        xticks([0 0.05 0.1 0.15 0.2 0.25 0.3])
        yticks([0 0.2 0.4 0.6 0.8 1]);
    end
end
hold off
set(fh1, 'position', [500, 200, 200*N, 200]);

