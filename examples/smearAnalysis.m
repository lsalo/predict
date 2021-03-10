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
smear = cell(Nsim, N);
for j=1:N
    [~, f.MatProps.SSFcBounds] = getSSFc(vcl, isClayVcl, zf, ...
                                         thick(j, :), f.Throw);
    for n=1:Nsim
        f.MatProps.Thick = getFaultThickness(f.Disp, 1);
        [f.Alpha, f.Delta, f.Zeta] = getFaultAngles(f.Throw, f.Dip, ...
                                                    f.MatProps.Thick);
        f.MatProps.ResFric = getResidualFrictionAngle(vcl);
        f.MatProps.SSFc = getSSFc(vcl, isClayVcl, zf, thick(j, :), f.Throw);
        smear{n, j} = Smear(vcl, isClayVcl, thick(j, :), thickap(j, :), ...
                            zf, f, 1);
    end
end

