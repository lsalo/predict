function [TFW, THW] = getThick(FW, HW, faultDip)
%
% get true layer thickness (measured perpendicular to the dip direction).
% It is assumed that FW.Thickness and HW.Thickness are the apparent
% thicknesses on the fault.
%
cFW = FW.Dip;
cHW = HW.Dip;
c = [cFW, cHW];
g = 90 - faultDip;

if g == 0 && all(c == 0) % vertical fault, hzntal layers
    TFW = FW.Thickness;
    THW = HW.Thickness;
elseif g == 90
    error('Fault dip cannot be 0')
elseif all(c == 0) % hzntal layers
    TFW = FW.Thickness .* cosd(g);
    THW = HW.Thickness .* cosd(g);
elseif g == 0      % vertical fault
    TFW = FW.Thickness .* cosd(cFW);
    THW = HW.Thickness .* cosd(cHW);
else               % dipping fault and layers
    % FW
    TFW = zeros(1, numel(FW.Id));
    if any(-cFW == g)
        id = -cFW == g;
        TFW(id) = FW.Thickness(id);
        TFW(~id) = FW.Thickness(~id) .* cosd(g + cFW(~id));
    else
        TFW = FW.Thickness .* cosd(g + cFW);
    end
    
    % HW
    THW = zeros(1, numel(HW.Id));
    if any(-cHW == g)
        id = -cHW == g;
        THW(id) = HW.Thickness(id);
        THW(~id) = HW.Thickness(~id) .* cosd(g + cHW(~id));
    else
        THW = HW.Thickness .* cosd(g + cHW);
    end
end
end