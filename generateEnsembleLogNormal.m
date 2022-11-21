% generate initial ensemble
function [Ens] = generateEnsembleLogNormal(num_members,mean,sd)

pd = makedist('Lognormal','mu',log(mean),'sigma',sd);
Ens = random(pd,num_members,1);

end