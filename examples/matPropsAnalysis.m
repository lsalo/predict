%% Example 1: Analysis of material properties
%
% 
% 

clear
close all

% Stratigraphic inputs
faultDisp = 100;
Nsim = 1000;
vcl = linspace(0, 1, 10);

% Plotting utils
latx = {'Interpreter', 'latex'};
sz = [14, 12];

%% Fault thickness

ft = getFaultThickness(faultDisp, Nsim);
DTratios = faultDisp./ft;

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
tiledlayout(2, 5, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile
histogram(phi(:, 1), edg_phi, 'Normalization', 'probability','FaceColor', ...
          [0.3 0.3 0.3])
xlabel('$\phi_\mathrm{r}$ [deg.]', latx{:}, 'fontSize', sz(2))
ylabel('P [-]', latx{:}, 'fontSize', sz(2))
title(['$N_\mathrm{sim} =$ ' num2str(Nsim), ' $\vert$ $V_\mathrm{cl}$ = ' ...
       num2str(vcl(1))], latx{:}, 'fontSize', sz(2))
xlim([0, 35]); ylim([0 1]); grid on; xticks([0 10 20 30]);
set(fh2, 'position', [500, 200, 1000, 400]);