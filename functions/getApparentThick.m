function [TapFW, TapHW] = getApparentThick(FW, HW, faultDip)
%
% get apparent thickness of layers at fault cutoffs
%
%Tap = [FS.FW.Thickness ./ (cosd(FS.FW.Dip)*sind(faultDip)), ...
%       FS.HW.Thickness ./ (cosd(FS.HW.Dip)*sind(faultDip
%
% For dipping layers:
% In the FW, dip angle (FS.FW.Dip) must be negative if it is
% dipping away from the fault (highest point is at the contact
% between the layer and the fault), and positive otherwise.
% In the HW, it has to be the opposite, i.e. the dip angle
% (FS.HW.Dip) must be negative if the layers are dipping
% towards the fault.
%

cFW = FW.Dip;
cHW = HW.Dip;
c = [cFW, cHW];
g = 90 - faultDip;

if g == 0 && all(c == 0) % vertical fault, hzntal layers
    TapFW = FW.Thickness;
    TapHW = HW.Thickness;
elseif g == 90
    error('Fault dip cannot be 0')
elseif all(c == 0) % hzntal layers
    TapFW = FW.Thickness ./ cosd(g);
    TapHW = HW.Thickness ./ cosd(g);
elseif g == 0      % vertical fault
    TapFW = FW.Thickness ./ cosd(cFW);
    TapHW = HW.Thickness ./ cosd(cHW);
else               % dipping fault and layers
    % FW
    TapFW = zeros(1, numel(FW.Id));
    if any(-cFW == g)
        id = -cFW == g;
        TapFW(id) = FW.Thickness(id);
        TapFW(~id) = FW.Thickness(~id) ./ cosd(g + cFW(~id));
    else
        TapFW = FW.Thickness ./ cosd(g + cFW);
    end
    
    % HW
    TapHW = zeros(1, numel(HW.Id));
    if any(-cHW == g)
        id = -cHW == g;
        TapHW(id) = HW.Thickness(id);
        TapHW(~id) = HW.Thickness(~id) ./ cosd(g + cHW(~id));
    else
        TapHW = HW.Thickness ./ cosd(g + cHW);
    end
end
end