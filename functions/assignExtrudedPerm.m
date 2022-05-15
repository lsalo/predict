function extrudedPerm = assignExtrudedPerm(extrudedPerm, faultSection, ...
                                           segLen, cellDim)
%
%
%

assert(mod(segLen/cellDim, 1)==0)
nklayers = segLen/cellDim;
idFirst0 = find(extrudedPerm==0, 1);
kyy = faultSection.Grid.perm(:,1);
kyz = faultSection.Grid.perm(:,2);
kzz = faultSection.Grid.perm(:,3);
kxx = faultSection.Grid.permy;
layerSize = size(kyy, 1);
idPerm3D = [1 4 5 6];                       % rotated about the x axis (along-strike)
extrudedPerm(idFirst0:(idFirst0-1)+nklayers*layerSize, idPerm3D) = ...
            repmat([kxx, kyy, kyz, kzz], nklayers, 1); 

end
     