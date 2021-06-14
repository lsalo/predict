function plotRunTime(ncores, telaps, description)
%
%
%
figure(321)
bar(ncores, telaps, 'FaceColor', [0.5 0.5 0.5])
xlabel('N cores')
ylabel('Elapsed time (s)')
title(description, 'fontSize', 10)
xlim([0 20]); xticks([1 4 16])
end
