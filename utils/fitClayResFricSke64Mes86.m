%
% (1) Fit data from Skempton (1964) and Mesri & Cepeda Diaz (1986)
% (2) Get 10% confidence intervals (covers visual scatter)
% (3) Add stochastic term modeled by beta distribution bounded by selected
%     confidence intervals.
%
clc, clear, close all

%% Fit data
% Skempton, Geotechnique (1964), Fig. 7
limits = [0 1; 0 40];
x      = [98.73 916.88];
y      = [83 736];
vals   = [235.33 283.54 391.21 441.76 480.73 548.39 621.39 661.84 727.52;
          571.48 488.80 392.06 342.58 376.66 343.40 233.80 295.68 263.68]; % [xcoord; ycoord], left to right.
ske64 = [(vals(1, :) - x(1)).*(limits(1,2) - limits(1,1)) / (x(2) - x(1)); ...
         (vals(2, :) - y(1)).*(limits(2,2) - limits(2,1)) / (y(2) - y(1))];

% Mesri & Cepeda-Diaz, Geotechnique (1986), Table 1 and Fig. 2
limits = [0 40];
y      = [15.3 133.64];
vals   = [92.29 44.68 77.70 57.42 56.80 94.85 44.58 41.68 44.92 62.92 50.81 ...
          34.31 29.49 30.21 41.85 30.08 42.17 29.77 29.96 28.58];
mes86 = [0.19 0.33 0.33 0.59 0.68 0.24 0.64 0.51 0.59 0.42 0.77 0.63 0.67 0.67 0.43 0.71 0.72 0.84 0.75 0.88;
         (vals - y(1)).*(limits(2) - limits(1)) / (y(2) - y(1))];

% Fit
[vcl, id] = sort([ske64(1,:) mes86(1,:)]);
phi = [ske64(2,:) mes86(2,:)];
phi = phi(id);
[xData, yData] = prepareCurveData(vcl, phi);
ft = fittype( 'exp2' );
[fitresult, gof] = fit(xData, yData, ft);

% Plot fit with data
fh1 = figure(1);
subplot(1,2,1)
plot(xData, yData, 'ok', 'markerSize', 3)
hold on
plot(0.2:0.02:1, fitresult(0.2:0.02:1), '-k', 'lineWidth', 1)
xlabel( 'V$_\mathrm{cl}$', 'Interpreter', 'latex', 'fontSize', 14);
ylabel( '$\phi_\mathrm{r}$ [deg.]', 'Interpreter', 'latex', 'fontSize', 14);
xlim([0.2 1])
xticks([0.2 0.4 0.6 0.8 1])
ylim([0 35])
%title('Fit for data with Vcl \in [0.2, 1]')
grid on

%% Confidence intervals
ci = confint(fitresult, 0.1);
lowLim = @(x) ci(1,1).*exp(ci(1,2).*x) + ci(1,3).*exp(ci(1,4).*x);
upLim = @(x) ci(2,1)*exp(ci(2,2)*x) + ci(2,3)*exp(ci(2,4)*x);


%% Final model with stochastic part
% betarnd to sample random beta values
% rescale to move to each interval
% Looking at the plot and fit, we have more values towards and below
% lowLim, so a < b. The actual values of a and b are subject to change.
a = 3;
b = 5;

% Generate random data
x = 0.25:0.05:1;                % examplve vcl input
x(x < 0.2) = [];                % this model only valid for vcl >= 0.2
lower = lowLim(x);              
upper = upLim(x);
normVals = betarnd(a,b,1,numel(x));
phiVar = lower + normVals.*(upper - lower);


% Plot
subplot(1,2,1)
plot(0.2:0.02:1, lowLim(0.2:0.02:1), '--k')
plot(x, phiVar, 'sb', 'markerFaceColor', 'b', 'markerSize', 3.5)
plot(0.2:0.02:1, upLim(0.2:0.02:1), '--k')
legend('Data', 'Fit', 'Limits', 'Samples', 'Location', 'NorthEast');

% Check distribution
phiDist = rescale(betarnd(a, b, 1, 10000), lowLim(0.6), upLim(0.6));
subplot(1,2,2)
histogram(phiDist, 10, 'Normalization', 'probability', 'faceColor', [0.3 0.3 0.3], 'edgeColor', 'none', 'faceAlpha', 1)
xlabel('$\phi_\mathrm{r}$ [deg.]', 'Interpreter', 'latex', 'fontSize', 14);
ylabel('$P$ [-]', 'Interpreter', 'latex', 'fontSize', 14);
text(8.5, 0.235, '$N = 10^4, V_{cl} = 0.6$', 'Interpreter', 'latex', 'fontSize', 12)
set(fh1, 'position', [500, 500, 500, 300]);

