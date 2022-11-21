%% Compute fraction of dislocation creep over the Antarctic ice sheet

clear all;

% activation energy values
Qdisplus = 155e3;
Qdisminus = 64e3;
Qgbsplus = 250e3;
Qgbsminus = 70e3;

% other parameter values
p = 9; % grain size exponent from grain size model
D = 0.03; % grain length scale from grain size model
dep = 0.001; % range of strain rates
ngbs = 1.8; % grain-boundary sliding stress exponent
ndis = 4; % dislocation creep stress exponent
nglen = 3; % glen's flow law exponent

% load ice temperature estimates (computed from
% main_ComputeIceTemperature_Antarctica.m)
load('temperatureest_antarctica_theta099_intermediaten.mat');

% load in rest of data needed (strain rates, ice thickness, surface velocity, surface mass
% balance)
[SRmat,Tmat,Vmat,SMBmat,Xfix,Yfix] = readData_Antarctica();

% load deformation map that defines the fractions of dislocation creep for
% given strain rate, temperature
title_string = sprintf('deformationmap_varyingstrainratetemp_d(T)_An_intermediaten_smoothtransition_tanh_p%d_dep%d_Qdisminus%d_Qgbsminus%d_Qdisplus%d_Qgbsplus%d.mat',p,dep,Qdisminus,Qgbsminus,Qdisplus,Qgbsplus);
V = load(title_string); % deterministic version

% interpolate
strainrate = logspace(-13,-6,100);
temperature = linspace(240,273,100);
[X,Y] = meshgrid(strainrate,temperature);
Vq = interp2(X,Y,V.frac_dis,SRmat,temp_map);
fracdismat_lookup = Vq;
fracdismat_lookup(abs(Vmat)>5000) = -9999;
fracdismat_lookup(Vmat < 30) = -9999;

% save
title_save = sprintf('partitioningest_antarctica_fracdis_Qstdev%d_p%d_D%d_theta099_intermediaten_incldepthavgtemp_Qdisminus%d_Qgbsminus%d_Qdisplus%d_Qgbsplus%d.mat',p,D,Qdisminus,Qgbsminus,Qdisplus,Qgbsplus);
save(title_save,'temp_map','fracdismat_lookup');

% plot
al = ones(size(fracdismat_lookup));
al(Vmat<30) = 0;
al(abs(Vmat)>5000) = 0;
[Tssing,R_V] = geotiffread('/Users/meghanaranganathan/Documents/MIT/Research/Data/SameExtent/racmo_tskin_fixedextent.tif');

xaxis = linspace(R_V.XWorldLimits(1)./1000, R_V.XWorldLimits(2)./1000, size(SRmat,2));
yaxis = linspace(R_V.YWorldLimits(1)./1000,R_V.YWorldLimits(2)./1000,size(SRmat,1));
load('groundingline.mat')

% plot background
Vmat2 = Vmat;
Vmat2(abs(Vmat)>5000) = -9999;
al2 = ones(size(Vmat2));
al2(Vmat2==-9999) = 0;
Vmat2(Vmat2 >= 0) = 1;

figure;
ax1 = axes;
imagesc(ax1,xaxis,yaxis,Vmat2);
alpha(ax1,al2);
caxis([0 1.2])
colormap(ax1,colorcet('l1','reverse',0))
ax2 = axes;
imagesc(xaxis,yaxis,real(fracdismat_lookup))
alpha(ax2,al)
colormap(ax2,colorcet('l17','reverse',0));
caxis([0 1])
hold on
for i=1:676
    plot(ax2,S(i).X/1000,-S(i).Y/1000,'Color',[0 0 0],'LineWidth',1)
    hold on;
end
linkaxes([ax1,ax2])
%cbar = colorbar;
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
ax1.Visible = 'off';
ax1.XTick = [];
ax1.YTick = [];
set(ax1,'FontSize',24,'FontWeight','b','GridColor','r');
set(ax2,'FontSize',24,'FontWeight','b','GridColor','r');
title(ax1,'\gamma') 
axis equal

