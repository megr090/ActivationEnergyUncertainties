function [Agbs] = computeGBSFlowRateParameter_TestParams(T,A0gbsplus,A0gbsminus,Qgbsplus,Qgbsminus)

R = 8.314; % J/mol K

% Grain Boundary Sliding
Agbs = zeros(size(T));
for i=1:length(T)
    if T(i) < 262
        A0 = A0gbsminus;
        Qcgbs = Qgbsminus;
    elseif T(i) >= 262
        A0 = A0gbsplus;
        Qcgbs = Qgbsplus;
    end
    Agbs(i) = A0*exp(-(Qcgbs./(R.*T(i))));
end

end
