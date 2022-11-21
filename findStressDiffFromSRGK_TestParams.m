function [diff] = findStressDiffFromSRGK_TestParams(x,strainrate,ngbs,ndis,nglen,D,p,theta,smb,z,H,A0displus,A0disminus,Qdisplus,Qdisminus,A0gbsplus,A0gbsminus,Qgbsplus,Qgbsminus,Ts)

[T] = computeIceTemperature(theta,nglen,z,H,smb,strainrate,x.*1e6,Ts);

Adis = computeDislocationFlowRateParameter_TestParams(T,strainrate,A0displus,A0disminus,Qdisplus,Qdisminus);
Agbs = computeGBSFlowRateParameter_TestParams(T,A0gbsplus,A0gbsminus,Qgbsplus,Qgbsminus);

[d] = computeGrainSize(T,theta,D,p,x.*1e6,strainrate);

strainratest = Adis.*x.^ndis+Agbs.*(d./1e3).^(-1.4).*x.^ngbs;
%[strainratest,sr_dislocation,sr_gbs] = computeStrainRate(T,d./1e3,x.*1e6,ndis,ngbs,strainrate);

diff = abs(strainrate-strainratest).*1e20;

end
