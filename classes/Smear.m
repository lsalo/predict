classdef Smear
    %
    % SUMMARY:
    % Get a 2D [x, z] smear object with:
    %   - smear thickness (Thick).
    %   - apparent smear thickness in fault (ThickInFault).
    %   - smear down-dip length (Length).
    %   - maximum down-dip smear segment length (SegLenMax).
    %   - Length of smear domains in the fault (Domain Length).
    %   - Smear fraction in the corresponding domain within the fault 
    %     (Psmear). 
    % Further description on the dimensions of clay smears are provided in 
    % the paper supplement
    %
    %
    % SYNOPSIS:
    %   smear = Smear(myFaultedSection, myFault, G, Nsim)
    %   * see examples for usage workflow in 2D and 3D
    %
    % DESCRIPTION:
    %   See constructor for details.
    %
    %
    % REQUIRED PARAMETERS:
    %   myFaultedSection: An instance of FaultedSection.
    %   myFault: Usually, an instance of Fault with valid properties/fields 
    %            (see Fault class documentation).
    %   G: Fault grid (MRST grid structure)
    %   Nsim: Number of smear property realizations (typically 1).
    %
    %
    % RETURNS:
    %   Class instance.
    %   
    % ____________________________________________________________________
    
    properties (SetAccess = protected)
        Thick
        ThickInFault
        Length      
        SegLenMax
        DomainLength
        Psmear
    end
    
    properties (SetAccess = protected, Hidden)
       ParentId
       ThickApparent
    end
    
    methods
        function smear = Smear(FS, fault, G, Nsim)
            %
            % instantiate a Smear object
            % 
            % INPUTS:
            %   FS:     an instance of FaultedSection
            %   fault:  an instance of Fault2D with material property values
            %   G:      fault grid (MRST grid)
            %   Nsim:   number of smear realizations
            %   
            %   * refer to examples for workflow
            %
            % Key references:
            %   Childs et al., GSLSP (2007)
            %   Egholm et al., Geology (2008)
            %   Grant, PG (2017)
            %
            % MODELS:
            %   Vcl < smear threshold: NaN (it does not apply).
            %
            %   Vcl >= smear threshold:
            %           - Thick: Average clay smear thickness of Egholm  
            %                    et al., Geology (2008), Eq. 5.
            %           - ThickInFault: apparent smear thickness in the dir
            %                           parallel to a diagonal from the 
            %                           top left or top right of the fault to 
            %                           bottom right or bottom left of 
            %                           fault.
            %           - Length: Total smear down-dip length within the fault, as 
            %                     obtained from the source layer thickness, 
            %                     smear angle in the fault (alpha), and 
            %                     SSFc.
            %           - SegLenMax: The maximum down-dip length of an individual 
            %                        segment of a smear within the fault.
            %                        The closer the SSFc is to the upper
            %                        bound, the more fragmented (same logic
            %                        as Grant, PG, 2017).
            %           - DomainLength: Length of smear domains within
            %                           fault, from FW cutoff to HW cutoff.
            %           - Psmear: Function of source layer thickness, SSFc, 
            %                     fault thickness, fault dip and fault 
            %                     throw. It can be seen as the probability 
            %                     of randomly finding a clay smear-filled 
            %                     cell within the corresponding smear 
            %                     domain in the fault. Or, the fraction
            %                     of a given domain filled by clay smear.
            %                     Similar logic to Childs et al., GSLSP
            %                     2007).
            %
            % OUTPUT:
            %  smear object with corresponding properties.
            %
            %--------------------------------------------------------------
            
            % Get required vars
            vcl = FS.Vcl;
            isClayVcl = FS.IsClayVcl;
            if ~isempty(FS.TotThick)    % clay layer(s) beyond throw window
                Tap = max(FS.Tap, FS.TotThick{1}, 'omitnan');
                T = max(FS.Thick, FS.TotThick{2}, 'omitnan');
            else
                Tap = FS.Tap;
                T = FS.Thick;
            end
            zf = FS.DepthFaulting;
            
            % Initialize
            idc = find(vcl >= isClayVcl);
            ids = find(vcl < isClayVcl);
            smear.ParentId = FS.ParentId(idc);
            assert(numel(zf) == 2)
            zf = [repelem(zf(1), sum(idc < FS.HW.Id(1))), ...
                  repelem(zf(2), sum(idc >= FS.HW.Id(1)))];
            N = numel(idc);
            [smear.Thick, smear.ThickInFault, smear.Length, ...
             smear.SegLenMax, smear.Psmear, smear.DomainLength] = ...
                        deal(nan(Nsim, N), nan(Nsim, N), nan(Nsim, N), ...
                        nan(Nsim, N), nan(Nsim, N), nan(Nsim, 1));
            SSFc  = fault.MatProps.ssfc(:, idc);
            phi   = fault.MatProps.resFric;
            epoin = cell2mat(FS.MatPropDistr.ssfc.range(idc)')';
            Lf = Tap(idc) + fault.Disp;     % Current Lf
            % Egholm et al. (2008) Lf (after fault formation):
            %Lf_eg  = (fault.throw + thick) ./ sind(theta_s);
            
            minVal = zeros(1, numel(zf)); 
            for n = 1:Nsim
                % (1) Thickness
                theta_c   = 45 + phi(n, idc);
                if numel(ids) > 0
                    theta_s = 45 + sum(phi(n, ids).*Tap(ids)./sum(Tap(ids)));
                else
                    theta_s = 45 + 30;
                end
                termCot = cotd(theta_c) - cotd(theta_s);
                smear.Thick(n, :)  = termCot .* T(idc).^2 ./ Lf;
                
                
                % (2) Thickness in Fault or Apparent thickness in fault
                smear.ThickInFault(n, :) = smear.Thick(n, :) ./ ...
                                           cosd(fault.Zeta(n));
                % if z is very close to 90, ThickInFault is > than true.
                TapMax  = sqrt(fault.Disp^2 + ...
                               fault.MatProps.thick(n)^2);
                smear.ThickInFault(n, smear.ThickInFault(n, :) > TapMax) = TapMax;
                
                
                % (3) Total smear length
                smear.Length(n, :) = (Tap(idc) ./ cosd(fault.Alpha(n))) .* ...
                                      SSFc(n, :);
                
                
                % (4) Maximum length of smear segments (if discontinuous). 
                smear.SegLenMax(n, :) = ((epoin(2, :) - SSFc(n, :)) ./ ...
                                         (epoin(2, :) - epoin(1, :))) .* ...
                                         smear.Length(n, :);
                
                % Max smear length is at least 10% of the total smear length or
                % 3 times the grid resolution in the length dimension,
                % whichever is larger. However, if very shallow faulting (50m
                % or less), smear fragmentation can be increased due to larger
                % shear strength of the clays compared to the sand. In that
                % case, the Max smear length is at least 5% of the total smear
                % length or 2 times the grid resolution in the length
                % dimension, whichever is larger.
                minVal(zf > 50) = max(0.1*smear.Length(n, zf>50), ...
                                      3*G.cellDim(end));
                minVal(zf <= 50) = max(0.05*smear.Length(n, zf<=50), ...
                                       2*G.cellDim(end));
                idBelowLim = smear.SegLenMax(n, :) < minVal;
                smear.SegLenMax(n, idBelowLim) = minVal(idBelowLim);
                
                % (5) Psmear. If < 1 then the smear is discontinuous, if NaN
                % there is no smear, and if Psmear is 1 then the smear is
                % continuous (i.e. it connects the source layer at each side of
                % the fault without breaks).
                smear.DomainLength(n) = fault.Throw / sind(fault.Delta(n));
                smear.Psmear(n, :) = smear.Length(n, :) ./ smear.DomainLength(n);
                if any(smear.Psmear(n, :) > 1)
                    smear.Psmear(n, smear.Psmear(n, :) > 1) = 1;
                end
            end
        
        end
         
    end
    
end


