% generate initial ensemble
function [Ens,mu,cov] = generateEnsemble(num_vars,num_members,mean,shift,sd)
cov = sd.^2*eye(num_vars);
mu = mean + randn(num_vars,1)*shift;
[ux,sx,vx] = svd(cov);
Ens = mu + ux*sqrt(sx)*randn(num_vars,num_members);
end