function plotId = selectSimId(arg, data)    
%
%
%

uperm = cell2mat(cellfun(@(x) x.Perm, data, 'UniformOutput', false));
if strcmp(arg, 'maxX')
    [~, plotId] = max(uperm(:, 1));
elseif strcmp(arg, 'maxZ')
    [~, plotId] = max(uperm(:, 3));
elseif strcmp(arg, 'minX')
    [~, plotId] = min(uperm(:, 1));
elseif strcmp(arg, 'minZ')
    [~, plotId] = min(uperm(:, 3));
elseif strcmp(arg, 'randm')
    plotId = randi(Nsim, 1);
end

end