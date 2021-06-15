function k = getPermSpe02(vcl, zf, zmax)
%
% Compute permeablity based on Sperrevik et al., NPSSP (2002)
%

a    = [8*10^4, 19.4, 0.00403, 0.0055, 12.5];
k    = a(1)*exp(-(a(2)*vcl + a(3)*zmax + ...
                  (a(4)*zf - a(5)).*(1-vcl).^7));    % [mD]