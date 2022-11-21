function [Adis] = computeDislocationFlowRateParameter_TestParams(T,sr,A0displus,A0disminus,Qdisplus,Qdisminus)

R = 8.314; % J/mol K

% Dislocation Creep
Adis = zeros(size(T)); % MPa^-4 s^-1 
for i=1:length(T)
    if T(i) < 262
        A0 = A0disminus;
        Qcdis = Qdisminus;
    elseif T(i) >= 262
        A0 = A0displus;
        Qcdis = Qdisplus;
    end
    Adis(i) = A0.*exp(-(Qcdis./(R.*T(i))));
end

%f = computeFabricSoftening(sr);
f = 1;
Adis = f.*Adis;

end

