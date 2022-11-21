function [T,grainsize,tau,frac_dis,sr_dislocation,sr_gbs,sr_gk] = computePartitioning_TestParams(tau_last,strainrate,ngbs,ndis,nglen,D,p,theta,smb,z,H,A0displus,A0disminus,Qdisplus,Qdisminus,A0gbsplus,A0gbsminus,Qgbsplus,Qgbsminus,Ts)

options = optimoptions('fsolve','Diagnostics','off','Display','off','MaxIterations',500,'FunctionTolerance',1e-10);
fun = @(x)findStressDiffFromSRGK_TestParams(x,strainrate,ngbs,ndis,nglen,D,p,theta,smb,z,H,A0displus,A0disminus,Qdisplus,Qdisminus,A0gbsplus,A0gbsminus,Qgbsplus,Qgbsminus,Ts);
x0 = tau_last;
[tau] = fsolve(fun,x0,options);
tau = abs(tau);

[T] = computeIceTemperature(theta,nglen,z,H,smb,strainrate,tau.*1e6,Ts);

[grainsize] = computeGrainSize(T,theta,D,p,tau.*1e6,strainrate);

[sr_gk,sr_dislocation,sr_gbs] = computeStrainRate_TestParams(T,grainsize,tau.*1e6,ndis,ngbs,strainrate,A0displus,A0disminus,Qdisplus,Qdisminus,A0gbsplus,A0gbsminus,Qgbsplus,Qgbsminus);

frac_dis = sr_dislocation./(sr_dislocation+sr_gbs);

end

