function [alpha, delta, zeta] = getFaultAngles(throw, dip, thick)
%
% DOCUMENTATION TBD
%

disp = throw / sind(dip);
gamma = 90 - dip;
b = thick ./ sind(dip) + throw*cotd(dip);
delta = atand(throw ./ b);
alpha = 90 - gamma - delta;
f1 = atand(disp ./ thick);
zeta = f1 - alpha;

end