%% Create plots of fraction of dislocation creep in probabilistic framework
% Find and plot fractio nof dislocation creep in idealized setup, subject to strain rate, and uncertainties in A0, Q 
% This code plots the full distributions

clear all;

% set standard deviations of A0, Q and strain rate
A0stdev = 1; % f_A0, the scaling on the standard deviation of the prefactor ensembles
Qstdev = 1; % f_Q, the scaling on the standard deviation of the activation energy ensembles
strainrate = 1e-8; % s^-1

% activation energy and prefactor parameters (means of distributions)
A0displus_exp = 6.96e23;
A0disminus_exp = 5e5;
Qdisplus_exp = 155e3;
Qdisminus_exp = 64e3;
A0gbsplus_exp = 8.5e37;
A0gbsminus_exp = 1.1e2;
Qgbsplus_exp = 250e3;
Qgbsminus_exp = 70e3;

% other parameters
p = 9; % grain growth exponent for grain size model
D = 0.03; % grain growth length scale for grain size model
ngbs = 1.8; % n value for grain-boundary sliding
ndis = 4; % n value for dislocation creep
nglen = 3; % n value in glen's flow law
H = 1000; % ice thickness
smb = 60; % surface mass balance
theta = 0.99; % energy partitioning between thermal and surface/strain energy
Ts = 248; % surface temperature
z = linspace(0,H,100); % depth profile

% generate ensembles
num_members = 1000; % number of ensemble members
num_vars = 1;
shift = 0;

% generate dislocation + ensemble
mean_qdisplus = Qdisplus_exp;
stdev_qdisplus = 2e4.*Qstdev;
mean_A0displus = A0displus_exp;
stdev_A0displus = A0stdev; % orders of magnitude
[Ens_qdisplus,mu,cov] = generateEnsemble(num_vars,num_members,mean_qdisplus,shift,stdev_qdisplus);
[Ens_A0displus] = generateEnsembleLogNormal(num_members,mean_A0displus,stdev_A0displus);
Qdisplus = Ens_qdisplus;
A0displus = Ens_A0displus;

% generate dislocation - ensemble
mean_qdisminus = Qdisminus_exp;
stdev_qdisminus = 1e4.*Qstdev;
mean_A0disminus = A0disminus_exp;
stdev_A0disminus = A0stdev;
[Ens_qdisminus,mu,cov] = generateEnsemble(num_vars,num_members,mean_qdisminus,shift,stdev_qdisminus);
[Ens_A0disminus] = generateEnsembleLogNormal(num_members,mean_A0disminus,stdev_A0disminus);
Qdisminus = Ens_qdisminus;
A0disminus = Ens_A0disminus;

% generate gbs + ensemble
mean_qgbsplus = Qgbsplus_exp;
stdev_qgbsplus = 8e4.*Qstdev;
mean_A0gbsplus = A0gbsplus_exp;
stdev_A0gbsplus = A0stdev;
[Ens_qgbsplus,mu,cov] = generateEnsemble(num_vars,num_members,mean_qgbsplus,shift,stdev_qgbsplus);
[Ens_A0gbsplus] = generateEnsembleLogNormal(num_members,mean_A0gbsplus,stdev_A0gbsplus);
Qgbsplus = Ens_qgbsplus;
A0gbsplus = Ens_A0gbsplus;

% generate gbs - ensemble
mean_qgbsminus = Qgbsminus_exp;
stdev_qgbsminus = 1e4.*Qstdev;
mean_A0gbsminus = A0gbsminus_exp;
stdev_A0gbsminus = A0stdev;
[Ens_qgbsminus,mu,cov] = generateEnsemble(num_vars,num_members,mean_qgbsminus,shift,stdev_qgbsminus);
[Ens_A0gbsminus] = generateEnsembleLogNormal(num_members,mean_A0gbsminus,stdev_A0gbsminus);
Qgbsminus = Ens_qgbsminus;
A0gbsminus = Ens_A0gbsminus;

% initialize matrices
Temp = zeros(length(num_members));
gsize = zeros(length(num_members));
fracdis = zeros(length(num_members));
Adis = zeros(length(num_members));
Agbs = zeros(length(num_members));
tau = zeros(length(num_members));
sr_gk = zeros(length(num_members));
sr_dislocation = zeros(length(num_members));
sr_gbs = zeros(length(num_members));

% compute partitioning of deformation mechanisms
for i=1:num_members
    if i==1
        tau_last = 0.1;
    end
    [Temp(i),gsize(i),tau(i),fracdis(i),sr_dislocation(i),sr_gbs(i),sr_gk(i)] = computePartitioning_TestParams(tau_last,strainrate,ngbs,ndis,nglen,D,p,theta,smb,z,H,A0displus(i),A0disminus(i),Qdisplus(i),Qdisminus(i),A0gbsplus(i),A0gbsminus(i),Qgbsplus(i),Qgbsminus(i),Ts);
    tau_last = tau(i);
    [Adis(i)] = computeDislocationFlowRateParameter_TestParams(Temp(i),strainrate,A0displus(i),A0disminus(i),Qdisplus(i),Qdisminus(i));
    [Agbs(i)] = computeGBSFlowRateParameter_TestParams(Temp(i),A0gbsplus(i),A0gbsminus(i),Qgbsplus(i),Qgbsminus(i));
end

% take only the members which estimate a sr_gk near inputted strain rate
idx = strainrate/10 < sr_gk & strainrate*10 > sr_gk;

% plot fraction of dislocation creep ensemble
figure;
h3 = histogram(real(fracdis(idx)));
set(gca,'FontSize',14,'FontWeight','b','GridColor','r');
xlabel('Fraction of Dislocation Creep') 
grid on
h3.Normalization = 'probability';
ylabel('Probability')
ylim([0 1])
%title('Strain Rate = 1*10^{-8} s^{-1}')
h3.FaceColor = [0.2 0.2 0.2];
h3.EdgeColor = [0.1 0.1 0.1];
h3.NumBins = 30;

% plot ice temperature ensemble 
figure;
h4 = histogram(real(Temp(idx)));
set(gca,'FontSize',14,'FontWeight','b','GridColor','r');
xlabel('Temperature (K)') 
grid on
h4.Normalization = 'probability';
ylabel('Probability')
%xlim([0 1e20])
h4.FaceColor = [0.2 0.2 0.2];
h4.EdgeColor = [0.1 0.1 0.1];
h4.NumBins = 30;

% plot grain size ensemble
figure;
h4 = histogram(real(gsize(idx)));
set(gca,'FontSize',14,'FontWeight','b','GridColor','r');
xlabel('Grain Size (mm)') 
grid on
h4.Normalization = 'probability';
ylabel('Probability')
xlim([0 100])
h4.FaceColor = [0.2 0.2 0.2];
h4.EdgeColor = [0.1 0.1 0.1];
h4.NumBins = 30;

% plot flow-rate parameter ensembles
figure;
h6 = histogram(log10(Adis(idx)));
hold on
h7 = histogram(log10(Agbs.*(gsize./1e3).^(-1.4)));
set(gca,'FontSize',14,'FontWeight','b','GridColor','r');
xlabel('A') 
grid on
h6.Normalization = 'probability';
ylabel('Probability')
h6.BinWidth = 1;
h7.BinWidth = 1;
h7.Normalization = 'probability';
ylabel('Probability')
legend('A_{dis}(T)','A_{gbs}(T,d)','FontSize',22)
xlim([-50 20])
xticks([-40 -20 0 20])
xticklabels({'10^{-40}','10^{-20}','10^{0}','10^{20}'})
