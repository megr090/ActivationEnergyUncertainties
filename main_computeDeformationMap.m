%% Varying Stress/Temperature, grain size computed from temperature
% In this section, we compute alpha for varying stress and temperature, for
% a grain size computed from temperature and stress.

clear all;

% define prefactors and activation energies
A0displus = 6.96e23;
A0disminus = 5e5;
Qdisplus = 155e3;
Qdisminus = 64e3;
A0gbsplus = 8.5e37;
A0gbsminus = 1.1e2;
Qgbsplus = 250e3;
Qgbsminus = 70e3;

% define other parameters
D = 0.03; % characteristic length scale for grain size model
p = 9; % grain growth exponent for grain size model
ndis = 4; % dislocation creep stress exponent
ngbs = 1.8; % grain boundary sliding stress exponent
nglen = 3; % glen's flow law stress exponent
theta = 0.99; % energy partitioning between thermal and stored energy
dep = 0.001; % range of strain rates

% range of strain rates and temperatures to compute gamma
strainrate = logspace(-13,-6,100);
temperature = linspace(240,273,100);

% initialize
frac_dis = zeros(length(strainrate),length(temperature));

% compute fraction of dislocation creep
for i=1:length(strainrate)
    tau_last = 0.1;
    tau_last_min = tau_last;
    tau_last_max = tau_last;
    for j=1:length(temperature)
        
        options = optimoptions('fsolve','Diagnostics','off','Display','off','MaxIterations',500,'FunctionTolerance',1e-10);
        fun = @(x)findStressDiffFromSRGK_constantT(x,theta,strainrate(i),ngbs,ndis,A0displus,A0disminus,Qdisplus,Qdisminus,A0gbsplus,A0gbsminus,Qgbsplus,Qgbsminus,temperature(j),D,p);
        x0 = 0.1;
        [tau,fval,exitflag] = fsolve(fun,x0,options);
        if exitflag < 1
            x0 = tau_last;
            [tau,fval,exitflag] = fsolve(fun,x0,options);
        end
        tau_last = tau;
        
        [d] = computeGrainSize(temperature(j),theta,D,p,tau.*1e6,strainrate(i));
        
        [sr_gk,sr_dislocation,sr_gbs] = computeStrainRate_SmoothTransition(temperature(j),d,tau.*1e6,ndis,ngbs,strainrate(i),A0displus,A0disminus,Qdisplus,Qdisminus,A0gbsplus,A0gbsminus,Qgbsplus,Qgbsminus);
        
        frac_dis(i,j) = sr_dislocation./(sr_dislocation+sr_gbs);

    end
    fprintf('Iteration %d of %d done \n',i,length(strainrate))
end

% plot
frac_dis = real(frac_dis);
frac_dis_flip = frac_dis(end:-1:1,:);

figure;
imagesc(frac_dis_flip)
hold on
cbar = colorbar;
colormap(colorcet('l17','reverse',0))
caxis([0 1])
set(gca,'FontSize',18,'FontWeight','b','GridColor','r');
ylabel('Strain Rate (s^{-1})')
xlabel('Temperature (K)')
% yticks([1 26 51 76 100])
% yticklabels({'10^{-7}','10^{-8}','10^{-9}','10^{-10}','10^{-11}'})
yticks([1 15 29 43 57 71 85 100])
yticklabels({'10^{-6}','10^{-7}','10^{-8}','10^{-9}','10^{-10}','10^{-11}','10^{-12}','10^{-13}'})
xticks([1 16 31 46 61 76 91])
xticklabels({'240','245','250','255','260','265','270'})
title('$$\alpha$$: Deformation Map, $$d_0(T)$$','Interpreter','Latex')

title_string = sprintf('deformationmap_varyingstrainratetemp_d(T)_An_intermediaten_smoothtransition_tanh_p%d_dep%d_Qdisminus%d_Qgbsminus%d_Qdisplus%d_Qgbsplus%d.mat',p,dep,Qdisminus,Qgbsminus,Qdisplus,Qgbsplus);
save(title_string,'frac_dis');