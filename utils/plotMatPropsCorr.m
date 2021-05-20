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
                            'UniformOutput', false));
end
thickRep = repmat(thick, numel(layerId), 1);
resFricResh = reshape(resFric, N*numel(layerId), 1);
SSFcResh = reshape(SSFc, N*numel(layerId), 1);
%kprime = reshape(kprime, N*numel(layerId), 1);
permResh = reshape(perm, N*numel(layerId), 1);
%poro = reshape(poro, N*numel(layerId), 1);
layerIds = repelem(layerId', N);
clayId = find(FS.Vcl >= FS.IsClayVcl);
id = ismember(layerIds, clayId);

%% Scatterhist Plots
% Utilities
latx = {'Interpreter', 'latex'};
sz = [14, 12];

% resFric and SSFc
mrk = '.';
lst = {'-', '-.', ':', '--', '-', '-.'};
colr = [0 0 0; 1 0 0; 0 0 1; 0.5 0.5 0.5; 1 0 1; 0 1 1];

figure(36)
sh = scatterhist(resFricResh(id), SSFcResh(id), 'Group', layerIds(id), 'Marker', '.', ...
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
sh = scatterhist(log10(faults{1}.Disp./thickRep), log10(permResh), 'Group', layerIds, ...
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


%% Correlation plots
% corrplot requires the Econometrics Toolbox, which requires the Statistics 
% and ML toolbox as well as the Optimization Toolbox.

idfig = 38;
for n=1:FS.ParentId(end)
    if FS.Vcl(n) < FS.IsClayVcl
        T = table(thick, resFric(:,1), poro(:,1), log10(perm(:,1)), kprime(:,1));
        T.Properties.VariableNames = ["Thick", "ResFric", "Poro", ...
                                      "Perm", "k'"];
    else
        T = table(thick, resFric(:,1), SSFc(:,n), poro(:,1), ...
                  log10(perm(:,1)), kprime(:,1));
        T.Properties.VariableNames = ["Thick", "ResFric", "SSFc", "n", ...
                                      "log k", "k'"]; 
    end
    figure(idfig+(n-1))                          
    [R, pval, h] = corrplot(T); grid on
    title(['Layer Id = '  num2str(n) ' $ \vert$ $V_\mathrm{cl}$ = ' num2str(FS.Vcl(n)) ...
           ' $ \vert$ $z_\mathrm{f}$ = ' num2str(FS.DepthFaulting(1)) ' m' ...
           ' $ \vert$ $z_\mathrm{max}$ = ' num2str(FS.DepthBurial(1)) ' m' ], ...
           latx{:}, 'fontSize', sz(2))
    
    % ------------------------------Fix axes-------------------------------
    % https://www.mathworks.com/matlabcentral/answers/172723-corrplot-plotting-on-strange-x-and-y-axes-values
    lineHandles = h(strcmp(get(h, 'type'), 'line'));       %get handles for scatter plots only
    % Loop through each scatter plot
    for i = 1:numel(lineHandles)
        % Axes limits
        x = lineHandles(i).XData;                         %x data 
        y = lineHandles(i).YData;                         %y data
        %limx = [fix(min(x)), fix(max(x))+1];
        %limy = [fix(min(y)), fix(max(y))+1];
        xlim(lineHandles(i).Parent, [min(x) max(x)]);    % set x limit to range of x data
        ylim(lineHandles(i).Parent, [min(y) max(y)]);    % set y limit to range of y data
        
        % Axes ticks
        %lineHandles(i).Parent.XTick = unique([limx(1) lineHandles(i).Parent.XTick limx(2)]);
        %lineHandles(i).Parent.YTick = unique([limy(1) lineHandles(i).Parent.YTick limy(2)]);
        
        % Marker type and colors
        lineHandles(i).MarkerEdgeColor = [0.5 0.5 0.5];
        
        % Line and correlation coeff
        lineHandles(i).Parent.Children(1).Color = [0 0 0];
        
        % Grid
        lineHandles(i).Parent.YGrid = 'on';
        lineHandles(i).Parent.XGrid = 'on';
    end
    
    % now take care of the x axis limits of the histogram plots
    histHandles = h(strcmp(get(h, 'type'), 'histogram'));     %handles to all hist plots
    
    % loop through hist plots
    for j = 1:numel(histHandles)
        % Bin limits
        x = histHandles(j).BinEdges;                         %bin edges
        xlim(histHandles(j).Parent, [min(x), max(x)]);       %set x limits
        
        % Colors
        histHandles(j).FaceColor = [0.7 0.7 0.7];
        histHandles(j).EdgeColor = [0 0 0];
    end
    %----------------------------------------------------------------------
end

end