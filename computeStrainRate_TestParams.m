function [sr_gk,sr_dislocation,sr_gbs] = computeStrainRate_TestParams(T,grainsize,tau,ndis,ngbs,prevsr,A0displus,A0disminus,Qdisplus,Qdisminus,A0gbsplus,A0gbsminus,Qgbsplus,Qgbsminus)
% this function computes strain rate in a Goldsby-Kohlstedt rheology, from
% a combination of grain boundary sliding and dislocation creep

sr_dislocation = computeDislocationCreepSR_TestParams(tau,ndis,T,prevsr,A0displus,A0disminus,Qdisplus,Qdisminus);
sr_gbs = computeGBSCreepSR_TestParams(tau,ngbs,grainsize,T,A0gbsplus,A0gbsminus,Qgbsplus,Qgbsminus);

sr_gk = sr_dislocation + sr_gbs;

end


