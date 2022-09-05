function plotRunTime(telaps)
%
%
%
fh = figure(321);
for n=2:3
    dim = n;
    subplot(1,2,n-1)
    if dim == 2
        bar(1:6, telaps{n-1}, 'FaceColor', [0.5 0.5 0.5])
        description = '2D | N = 1000 | Intel(R) Xeon(R) Gold 6230 CPU @ 2.1 GHz';
        xticks([1 2 3 4 5 6]); xticklabels([1 2 4 8 16 32])
        xlim([0 7]);
    elseif dim == 3
        bar(1:4, telaps{n-1}, 'FaceColor', [0.5 0.5 0.5])
        description = '3D | N = 1000 | Intel(R) Xeon(R) Gold 6230 CPU @ 2.1 GHz';
        xticks([1 2 3 4]); xticklabels([4 8 16 32])
        xlim([0 5])
    end
    ylim([0 1800]); yticks(0:200:1800)
    xlabel('# of cores')
    ylabel('Elapsed time (s)')
    title(description, 'fontSize', 8, 'fontweight', 'normal')
    grid on
    set(fh, 'position', [100, 100, 850, 300]);
end

end
