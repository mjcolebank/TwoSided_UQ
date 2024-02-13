% This is a driver file that will run the fluids model from c++ and
% fortran by passing parameter values needed by the model. This code is a
% working file and comes with NO GUARANTEES.
%
% Authors: MU Qureshi, MJ Colebank, M Bartolo, MS Olufsen
%
% Last edited: 2/14/2020, MJC
%%
clear; clc; close all;
%% Ensure the make file is compiled
! chmod +x treebranch
! make clean
! make veryclean
! make
if exist('alpha_beta')~=7
   mkdir('alpha_beta'); 
end
%

load micro_test.mat
r_vals = [0.757, 0.514,0.433, 0.293, 0.829, 0.562, 0.460, 0.610];
%% Decide which artery and vein to include for the analysis
ves_id = 2; 
which_samp = 10;

% Also include which admittance file we want to load
r = r_vals(ves_id);
num_art_offset = 7;
num_ven_offset = 19; 
pA_term = squeeze(p_test(:,num_art_offset+ves_id,which_samp));
qA_term = squeeze(q_test(:,num_art_offset+ves_id,which_samp));
pV_term = squeeze(p_test(:,num_ven_offset+ves_id,which_samp));
qV_term = -squeeze(q_test(:,num_ven_offset+ves_id,which_samp));
par_samp = par_test(which_samp,:);

%% Increase size of inputs if necessary
new_tmstps = 512;
pA_term = interp1(linspace(0,1,512),pA_term,linspace(0,1,new_tmstps))';
pV_term = interp1(linspace(0,1,512),pV_term,linspace(0,1,new_tmstps))';
qA_term = interp1(linspace(0,1,512),qA_term,linspace(0,1,new_tmstps))';
qV_term = interp1(linspace(0,1,512),qV_term,linspace(0,1,new_tmstps))';

%%
qp_art = [qA_term pA_term];
qp_ven = [qV_term pV_term];

%% Define parameters

fa1 = 0;
fa2 = 0;
fa3 = par_samp(3);
fv1 = 0;
fv2 = 0;
fv3 = par_samp(3);
alpha_b = par_samp(4);
beta_b  =  par_samp(5);
lrrA    = par_samp(6);
lrrV    = par_samp(7);
r_min   = par_samp(8);
r_root  = r; %terminal artery/vein radius, largest radius
maxgen  = 100; %maximum # of generations
tmstps = length(qA_term); %# of time points
period = 0.85; %length of heartbeat
mu     = 0.032; %large vessel viscosity
rho    = 1.055; %large vessel density
%% For the time being, calculate the maximum alpha and beta 
alpha = alpha_b;
beta  = beta_b;

alpha_vec = r_root.*alpha.^[0:100];
beta_vec  = r_root.*beta.^[0:100];
alpha_max = 3%find(alpha_vec<r_min,1);%-1;
beta_max  = 3%find(beta_vec<r_min,1);%-1;

%%
pars = [period mu rho fa1 fa2 fa3 fv1 fv2 fv3 r_root ...
        r_min  lrrA lrrV alpha_b beta_b maxgen tmstps alpha_max beta_max ves_id]';
    
% Write to file
dlmwrite('parameters_MJC.dat',pars); 

%% Find the file in the previous folder

% load('../output_artery.mat');
% load('../output_vein.mat');
% copyfile('../Admit_1_alpha.out');
% copyfile('../Admit_1_beta.out');



qp_art = [qA_term pA_term];
qp_ven = [qV_term pV_term];

%%
q_art_root = qp_art(1:tmstps,1);
p_art_root = qp_art(1:tmstps,2);
q_ven_root = qp_ven(1:tmstps,1);
p_ven_root = qp_ven(1:tmstps,2);

%% save
dlmwrite('qp_art.dat',qp_art);
dlmwrite('qp_ven.dat',qp_ven);

%% run the fortran file
flag = unix(sprintf('treebranch.exe'));
% flag = unix(sprintf('treebranch_mean.exe'));


%% All files are printed to the folder "alpha_beta". Read in those results 
% Now loop until you reach alpha or beta max
alpha_data = zeros(tmstps,alpha_max,4); 
beta_data  = zeros(tmstps,beta_max,4);
cd('alpha_beta');
for i=1:alpha_max
    name = strcat('alpha_',num2str(i-1),'.2d');
    data = dlmread(name);
    alpha_data(:,i,:) = data;
end
for i=1:beta_max
    name = strcat('beta_',num2str(i-1),'.2d');
    data = dlmread(name);
    beta_data(:,i,:) = data;
end
cd('../');

%% MJC plot pressure and flow
t = linspace(0,1,tmstps);
alpha_max = size(alpha_data,2);
beta_max  = size(beta_data,2);
P_art_alpha = alpha_data(:,:,1); P_ven_alpha = alpha_data(:,:,2);
Q_art_alpha = alpha_data(:,:,3); Q_ven_alpha = -1.*alpha_data(:,:,4);

P_art_beta = beta_data(:,:,1); P_ven_beta = beta_data(:,:,2);
Q_art_beta = beta_data(:,:,3); Q_ven_beta = -1.*beta_data(:,:,4);


color_fade_red = fliplr([linspace(0, 1,alpha_max); linspace(0.05, 1,alpha_max)*0; linspace(0.05, 1,alpha_max)*0]);
color_fade_pink = [linspace(0.9, 0,beta_max); linspace(0.3, 0,beta_max); linspace(0, 0,beta_max)];

color_fade_blue = [linspace(0.1, 0,beta_max); linspace(0.4, 0,beta_max); linspace(1, 0,beta_max)];
color_fade_lightblue = [linspace(0.5, 0,alpha_max); linspace(0.8, 0,alpha_max); linspace(1, 0,alpha_max)];

root_data = [p_art_root p_ven_root q_art_root q_ven_root];

figure(1);clf; 
subplot(2,2,1);
hold on;
for i=1:alpha_max
plot(t,P_art_alpha(:,i),'Color', color_fade_red(:,i),'LineWidth',5*exp(-i/7));
end
plot(t,root_data(:,1),'k','LineWidth',3);
set(gca,'FontSize',20); 
ylabel('Pressure [mmHg]');
xlabel('Time [s]');
title('Arteries (\alpha Branch)')

figure(1);
subplot(2,2,2);
hold on;
for i=1:beta_max
plot(t,P_art_beta(:,i),'Color', color_fade_pink(:,i),'LineWidth',5*exp(-i/7));
end
plot(t,root_data(:,1),'k','LineWidth',3);
set(gca,'FontSize',20); 
ylabel('Pressure [mmHg]');
xlabel('Time [s]');
title('Arteries (\beta Branch)')

figure(1);
subplot(2,2,3);
hold on;
for i=1:alpha_max
plot(t,P_ven_alpha(:,i),'Color', color_fade_lightblue(:,i),'LineWidth',5*exp(-i/7));
end
plot(t,root_data(:,2),'k','LineWidth',3);
set(gca,'FontSize',20); 
ylabel('Pressure [mmHg]');
xlabel('Time [s]');
title('Veins (\alpha Branch)')

figure(1);
subplot(2,2,4);
hold on;
for i=1:beta_max
plot(t,P_ven_beta(:,i),'Color', color_fade_blue(:,i),'LineWidth',5*exp(-i/7));
end
plot(t,root_data(:,2),'k','LineWidth',3);
set(gca,'FontSize',20); 
ylabel('Pressure [mmHg]');
xlabel('Time [s]');
title('Veins (\beta Branch)')


figure(2);clf; 
subplot(2,2,1);
hold on;
for i=1:alpha_max
plot(t,Q_art_alpha(:,i),'Color', color_fade_red(:,i),'LineWidth',5*exp(-i/7));
end
plot(t,root_data(:,3),'k','LineWidth',3);
set(gca,'FontSize',20); 
ylabel('Flow [ml/s]');
xlabel('Time [s]');
title('Arteries (\alpha Branch)')


subplot(2,2,2);
hold on;
for i=1:beta_max
plot(t,Q_art_beta(:,i),'Color', color_fade_pink(:,i),'LineWidth',5*exp(-i/7));
end
plot(t,root_data(:,3),'k','LineWidth',3);
set(gca,'FontSize',20); 
ylabel('Flow [ml/s]');
xlabel('Time [s]');
title('Arteries (\beta Branch)')

subplot(2,2,3);
hold on;
for i=1:alpha_max
plot(t,Q_ven_alpha(:,i),'Color', color_fade_lightblue(:,i),'LineWidth',5*exp(-i/7));
end
plot(t,root_data(:,4),'k','LineWidth',3);
set(gca,'FontSize',20); 
ylabel('Flow [ml/s]');
xlabel('Time [s]');
title('Veins (\alpha Branch)')

subplot(2,2,4);
hold on;
for i=1:beta_max
plot(t,Q_ven_beta(:,i),'Color', color_fade_blue(:,i),'LineWidth',5*exp(-i/7));
end
plot(t,root_data(:,4),'k','LineWidth',3);
set(gca,'FontSize',20); 
ylabel('Flow [ml/s]');
xlabel('Time [s]');
title('Veins (\beta Branch)')
%% Pass the alpha and beta data into the shear stress function
return
pars_shear = [alpha beta period mu rho fa1 fa2 fa3 fv1 fv2 fv3 r_root ...
        r_min alpha_b beta_b lrrA lrrV maxgen length(p_art_root)]';
get_shear_micro(alpha_data,beta_data,pars_shear,root_data);

