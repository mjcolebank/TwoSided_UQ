% This function will calculate the shear stress in the microcirculation
% using the alpha and beta branch data provided. The alpha and beta data
% structures contain arterial pressure, venous pressure, arterial flow,
% venous flow. THERE ARE NO OPTIONAL ARGUMENTS
% This code is a working file and comes with NO GUARANTEES.
%
% Authors: MU Qureshi, MJ Colebank, M Bartolo, MS Olufsen
%
% Last edited: 1/27/2020, MJC

function A_art_alpha = get_shear_micro(alpha_data,beta_data,pars,root_data)
alpha_max = size(alpha_data,2);
beta_max  = size(beta_data,2);
%% Define all the variables passed in
alpha = pars(1); beta = pars(2); T = pars(3);
mu = pars(4); rho = pars(5);
fa1 = pars(6); fa2 = pars(7); fa3 = pars(8);
fv1 = pars(9); fv2 = pars(10); fv3 = pars(11);
r_root = pars(12); r_min = pars(13);
lrrA = pars(16);
lrrV = pars(17);

tmstps = pars(19);
t = linspace(0,T,tmstps);
dt = diff(t(1:2));
conv = 1333.22; % Conversion factor to go from mmHg to SI units
%% Rename the pressures and flows
P_art_alpha = alpha_data(:,:,1); P_ven_alpha = alpha_data(:,:,2);
Q_art_alpha = alpha_data(:,:,3); Q_ven_alpha = -1.*alpha_data(:,:,4);

P_art_beta = beta_data(:,:,1); P_ven_beta = beta_data(:,:,2);
Q_art_beta = beta_data(:,:,3); Q_ven_beta = -1.*beta_data(:,:,4);

%% First, calculate the radii, areas, and lengths for the alpha and beta branches
r_alpha = r_root.*alpha.^[0:alpha_max-1];
r_beta  = r_root.*beta.^[0:beta_max-1];
A0_alpha = pi.*r_alpha.^2;
A0_beta  = pi.*r_beta.^2;

% NOTE: The lrr relation will change in future versions. These current
% functions are from literature (???).
% for i = 1:length(r_alpha)
%     if(r_alpha(i)<0.005)
%       L_art_alpha(i) = 1.88*r_alpha(i)^0.47;
%     else
%       L_art_alpha(i) = lrrA*r_alpha(i)^1.1;
%     end
% end
% 
% for i = 1:length(r_beta)
%     if(r_beta(i)<0.005)
%       L_art_beta(i) = 1.88*r_beta(i)^0.47;
%     else
%       L_art_beta(i) = lrrA*r_beta(i)^1.1;
%     end
% end
L_ven_alpha = lrrA.*r_alpha;
L_ven_beta  = lrrA.*r_beta;
L_ven_alpha = lrrV.*r_alpha;
L_ven_beta  = lrrV.*r_beta;
%% Calculate the stiffness
stiff_art_alpha = (4/3).*(fa1*exp(fa2*r_alpha) + fa3);
stiff_art_beta  = (4/3).*(fa1*exp(fa2*r_beta) + fa3);

stiff_ven_alpha = (4/3).*(fv1*exp(fv2*r_alpha) + fv3);
stiff_ven_beta  = (4/3).*(fv1*exp(fv2*r_beta) + fv3);

%% Define the viscosity, which has a dynamic dependence on radius
% The function is originally defined in microns
D = 2*r_alpha*10^4;
eta = 6*exp(-0.085*2*D)+3.2-2.44*exp(-0.06*(D.^0.645));
chi = (D./(D- 1.1)).^2;
mu_eff  = (1 + (eta - 1).*chi).*chi;
mu_alpha = mu*mu_eff/3.2;

D = 2*r_beta*10^4;
eta = 6*exp(-0.085*2*D)+3.2-2.44*exp(-0.06*(D.^0.645));
chi = (D./(D- 1.1)).^2;
mu_eff  = (1 + (eta - 1).*chi).*chi;
mu_beta = mu*mu_eff/3.2;

%% Calculate the dynamic area of the vessel using the wall model defintion

% Linear B
% A_art_alpha = A0_alpha.*(1+P_art_alpha*conv./stiff_art_alpha).^2 ;%
% A_art_beta  = A0_beta.*(1+P_art_beta*conv./stiff_art_beta).^2 ;%  

A_art_alpha = A0_alpha./((1-P_art_alpha*conv ./ stiff_art_alpha).^2);
A_art_beta = A0_beta./ (1-P_art_beta*conv  ./ stiff_art_beta).^2;

% A_ven_alpha = A0_alpha.*(1+P_ven_alpha*conv./stiff_ven_alpha).^2 ;
% A_ven_beta  = A0_beta.*(1+P_ven_beta*conv./stiff_ven_beta).^2 ;%

A_ven_alpha = A0_alpha./(1-P_ven_alpha*conv./stiff_ven_alpha).^2;
A_ven_beta = A0_beta./(1-P_ven_beta*conv./stiff_ven_beta).^2;

figure(100)
plot(A_art_alpha, 'r','linewidth',2)
set(gca,'fontsize',16)
axis tight
hold on
plot(A_art_beta, 'b','linewidth',2)
title('Artery \alpha and \beta area')
set(gca,'fontsize',16)
axis tight

figure(102)
plot(A_art_beta, 'r','linewidth',2)
title('Arteries Area (\beta)')
set(gca,'fontsize',16)
axis tight

figure(103)
plot(A_ven_beta, 'b','linewidth',2)
title('Veins Area (\beta)')
set(gca,'fontsize',16)
axis tight

%% Calculate the cyclic stretch in the arteries & veins
r_art_alpha = sqrt(A_art_alpha./pi); r_art_beta = sqrt(A_art_beta./pi);
r_ven_alpha = sqrt(A_ven_alpha./pi); r_ven_beta = sqrt(A_ven_beta./pi);

CS_art_alpha = 100*(max(r_art_alpha)-min(r_art_alpha))./min(r_art_alpha);
CS_ven_alpha = 100*(max(r_ven_alpha)-min(r_ven_alpha))./min(r_ven_alpha);

CS_art_beta = 100*(max(r_art_beta)-min(r_art_beta))./min(r_art_beta);
CS_ven_beta = 100*(max(r_ven_beta)-min(r_ven_beta))./min(r_ven_beta);

% Mean shear
TM_art_alpha = 4.*mu_alpha.*mean(Q_art_alpha)./(pi.*mean(r_art_alpha).^3);
TM_art_beta = 4.*mu_beta.*mean(Q_art_beta)./(pi.*mean(r_art_beta).^3);
TM_ven_alpha = 4.*mu_alpha.*mean(Q_ven_alpha)./(pi.*mean(r_ven_alpha).^3);
TM_ven_beta = 4.*mu_beta.*mean(Q_ven_beta)./(pi.*mean(r_ven_beta).^3);

%% Calculate the shear stress as done in Zamir 
% Alpha Branch
omega = 2*pi;
rho = 1.055;
Omega = sqrt((rho.*omega)./mu_alpha).*r_alpha;
lambda = ((1i-1)./sqrt(2)).*Omega;

%Arteries - alpha
T_art_alpha = ((lambda.* mu_alpha) ./ (pi .* r_alpha.^3)).*...
        (besselj(1,lambda)./((besselj(0,lambda)-(2./lambda).*besselj(1,lambda)))).*...
        Q_art_alpha;
T_art_alpha = abs(T_art_alpha);

%Veins - alpha 
T_ven_alpha = ((lambda.* mu_alpha) ./ (pi .* r_alpha.^3)).*...
        (besselj(1,lambda)./((besselj(0,lambda)-(2./lambda).*besselj(1,lambda)))).*...
        Q_ven_alpha;
T_ven_alpha = abs(T_ven_alpha);

%Beta branch
Omega = sqrt((rho.*omega)./mu_beta).*r_beta;
lambda = ((1i-1)./sqrt(2)).*Omega;

%Arteries - beta
T_art_beta = ((lambda.* mu_beta) ./ (pi .* r_beta.^3)).*...
        (besselj(1,lambda)./((besselj(0,lambda)-(2./lambda).*besselj(1,lambda)))).*...
        Q_art_beta;
T_art_beta = abs(T_art_beta);

%Veins - beta
T_ven_beta = ((lambda.* mu_beta) ./ (pi .* r_beta.^3)).*...
        (besselj(1,lambda)./((besselj(0,lambda)-(2./lambda).*besselj(1,lambda)))).*...
        Q_ven_beta;
T_ven_beta = abs(T_ven_beta);
%% Now begin plotting all the above results:
color_fade_red = fliplr([linspace(0, 1,alpha_max); linspace(0.05, 1,alpha_max)*0; linspace(0.05, 1,alpha_max)*0]);
color_fade_pink = [linspace(0.9, 0,beta_max); linspace(0.3, 0,beta_max); linspace(0, 0,beta_max)];

color_fade_blue = [linspace(0.1, 0,beta_max); linspace(0.4, 0,beta_max); linspace(1, 0,beta_max)];
color_fade_lightblue = [linspace(0.5, 0,alpha_max); linspace(0.8, 0,alpha_max); linspace(1, 0,alpha_max)];

%% Pressure
%{
figure; 
for i=1:alpha_max
 hold on;
h1=plot(t,P_art_alpha(:,1),'Color', color_fade_lightblue(:,1),'LineWidth',5*exp(-1/7));
h2=plot(t,P_art_alpha(:,5),'Color', color_fade_lightblue(:,5),'LineWidth',5*exp(-5/7));
h3=plot(t,P_art_alpha(:,10),'Color', color_fade_lightblue(:,10),'LineWidth',5*exp(-10/7));
h4=plot(t,P_art_alpha(:,15),'Color', color_fade_lightblue(:,15),'LineWidth',5*exp(-15/7));
h5=plot(t,P_art_alpha(:,alpha_max),'Color', color_fade_lightblue(:,alpha_max),'LineWidth',5*exp(-alpha_max/7));

plot(t,P_art_alpha(:,i),'Color', color_fade_lightblue(:,i),'LineWidth',5*exp(-i/7));
end
%plot(t,root_data(:,1),'-.c','LineWidth',3);
set(gca,'FontSize',20); 
ylabel('Pressure [mmHg]');
xlabel('Time [s]');
title('Arteries (\alpha Branch)')
leg=legend([h1,h2,h3,h4,h5],'Largest Vessel','','','','Smallest Vessel');
leg.FontSize = 15;

figure;
for i=1:alpha_max
 hold on;
h1=plot(t,P_ven_alpha(:,1),'Color', color_fade_red(:,1),'LineWidth',5*exp(-1/7));
h2=plot(t,P_ven_alpha(:,5),'Color', color_fade_red(:,5),'LineWidth',5*exp(-5/7));
h3=plot(t,P_ven_alpha(:,10),'Color', color_fade_red(:,10),'LineWidth',5*exp(-10/7));
h4=plot(t,P_ven_alpha(:,15),'Color', color_fade_red(:,15),'LineWidth',5*exp(-15/7));
h5=plot(t,P_ven_alpha(:,alpha_max),'Color', color_fade_red(:,alpha_max),'LineWidth',5*exp(-alpha_max/7));

plot(t,P_ven_alpha(:,i),'Color',color_fade_red(:,i),'LineWidth',5*exp(-i/7));
end
%plot(t,root_data(:,2),'-.c','LineWidth',3);
set(gca,'FontSize',20); 
ylabel('Pressure [mmHg]');
xlabel('Time [s]');
title('Veins (\alpha Branch)')
leg=legend([h5,h4,h3,h2,h1],'Smallest Vessel','','','','Largest Vessel');
leg.FontSize = 15;
ylim([0 35])

figure;
for i=1:beta_max
 hold on;
h1=plot(t,P_art_beta(:,1),'Color', color_fade_blue(:,1),'LineWidth',5*exp(-1/7));
h2=plot(t,P_art_beta(:,4),'Color', color_fade_blue(:,4),'LineWidth',5*exp(-4/7));
h3=plot(t,P_art_beta(:,7),'Color', color_fade_blue(:,7),'LineWidth',5*exp(-7/7));
h4=plot(t,P_art_beta(:,11),'Color', color_fade_blue(:,11),'LineWidth',5*exp(-11/7));
h5=plot(t,P_art_beta(:,beta_max),'Color', color_fade_blue(:,beta_max),'LineWidth',5*exp(-beta_max/7));

plot(t,P_art_beta(:,i),'Color',color_fade_blue(:,i),'LineWidth',5*exp(-i/7));
end
%plot(t,root_data(:,1),'-.c','LineWidth',3);
leg=legend([h1,h2,h3,h4,h5],'Largest Vessel','','','','Smallest Vessel');
leg.FontSize = 15;
set(gca,'FontSize',20); 
ylabel('Pressure [mmHg]');
xlabel('Time [s]');
title('Arteries (\beta Branch)')

figure;
for i=1:beta_max
 hold on;
h1=plot(t,P_ven_beta(:,1),'Color',color_fade_pink(:,1),'LineWidth',5*exp(-1/7));
h2=plot(t,P_ven_beta(:,4),'Color', color_fade_pink(:,4),'LineWidth',5*exp(-4/7));
h3=plot(t,P_ven_beta(:,7),'Color', color_fade_pink(:,7),'LineWidth',5*exp(-7/7));
h4=plot(t,P_ven_beta(:,11),'Color', color_fade_pink(:,11),'LineWidth',5*exp(-11/7));
h5=plot(t,P_ven_beta(:,beta_max),'Color', color_fade_pink(:,beta_max),'LineWidth',5*exp(-beta_max/7));

plot(t,P_ven_beta(:,i),'Color',color_fade_pink(:,i),'LineWidth',5*exp(-i/7));
end
%plot(t,root_data(:,2),'-.c','LineWidth',3);
set(gca,'FontSize',20); 
ylabel('Pressure [mmHg]');
xlabel('Time [s]');
leg=legend([h5,h4,h3,h2,h1],'Smallest Vessel','','','','Largest Vessel');
leg.FontSize = 15;%ylim([0 20])
title('Veins (\beta Branch)')
ylim([0 35])

figure; clf; 
%subplot(1,2,1);
plot(-r_alpha,mean(P_art_alpha),'-r.','LineWidth',3,'MarkerSize',36);
hold on; grid on;
plot(-r_beta,mean(P_art_beta) ,'-bo','LineWidth',3,'MarkerSize',10);
set(gca,'FontSize',20);
ylabel('Mean Pressure (mmHg)');
%set(gca,'XDir','reverse');
%subplot(1,2,2);
plot(r_alpha,mean(P_ven_alpha),'-r.','LineWidth',3,'MarkerSize',36);
hold on; grid on;
plot(r_beta,mean(P_ven_beta),'-bo','LineWidth',3,'MarkerSize',10);
set(gca,'FontSize',20);
xlabel('Radius (mm)');
xlim([-0.5,0.5])
set(gca, 'XTick', -0.5:0.25:0.5)
legend('\alpha Branch','\beta Branch', 'location','southwest')
axis tight
grid on
%}
%% Flow
%{
figure;
for i=1:alpha_max
 hold on;
 h1=plot(t,Q_art_alpha(:,1),'Color',color_fade_lightblue(:,1),'LineWidth',5*exp(-1/7));
 h2=plot(t,Q_art_alpha(:,5),'Color',color_fade_lightblue(:,5),'LineWidth',5*exp(-5/7));
 h3=plot(t,Q_art_alpha(:,10),'Color',color_fade_lightblue(:,10),'LineWidth',5*exp(-10/7));
 h4=plot(t,Q_art_alpha(:,15),'Color',color_fade_lightblue(:,15),'LineWidth',5*exp(-15/7));
 h5=plot(t,Q_art_alpha(:,22),'Color',color_fade_lightblue(:,22),'LineWidth',5*exp(-22/7));
 plot(t,Q_art_alpha(:,i),'Color',color_fade_lightblue(:,i),'LineWidth',5*exp(-i/7));
end
%plot(t,root_data(:,3),'-.c','LineWidth',3);
set(gca,'FontSize',20); 
title('Arteries (\alpha Branch)')
ylabel('Flow [mL/s]');
xlabel('Time [s]');
leg=legend([h1,h2,h3,h4,h5],'Largest Vessel','','','','Smallest Vessel','location','northeast');
leg.FontSize = 15;
ylim([0 50])

figure;
for i=1:alpha_max
hold on;
h1=plot(t,Q_ven_alpha(:,1),'Color',color_fade_red(:,1),'LineWidth',5*exp(-1/7));
h2=plot(t,Q_ven_alpha(:,5),'Color',color_fade_red(:,5),'LineWidth',5*exp(-5/7));
h3=plot(t,Q_ven_alpha(:,8),'Color',color_fade_red(:,10),'LineWidth',5*exp(-10/7));
h4=plot(t,Q_ven_alpha(:,15),'Color',color_fade_red(:,15),'LineWidth',5*exp(-15/7));
h5=plot(t,Q_ven_alpha(:,22),'Color',color_fade_red(:,22),'LineWidth',5*exp(-22/7));

plot(t,Q_ven_alpha(:,i),'Color',color_fade_red(:,i),'LineWidth',5*exp(-i/7));
end
%plot(t,root_data(:,4),'-.c','LineWidth',3);
set(gca,'FontSize',20); 
title('Veins (\alpha Branch)')
ylabel('Flow [mL/s]');
xlabel('Time [s]');
leg=legend([h1,h2,h3,h4,h5],'Largest Vessel','','','','Smallest Vessel','location','northeast');
leg.FontSize = 15;
ylim([0 50])

figure; 
for i=1:beta_max
 hold on;
h1=plot(t,Q_art_beta(:,1),'Color',color_fade_blue(:,1),'LineWidth',5*exp(-1/7));
h2=plot(t,Q_art_beta(:,4),'Color',color_fade_blue(:,4),'LineWidth',5*exp(-4/7));
h3=plot(t,Q_art_beta(:,7),'Color',color_fade_blue(:,7),'LineWidth',5*exp(-7/7));
h4=plot(t,Q_art_beta(:,11),'Color',color_fade_blue(:,11),'LineWidth',5*exp(-11/7));
h5=plot(t,Q_art_beta(:,beta_max),'Color',color_fade_blue(:,beta_max),'LineWidth',5*exp(-14/7));

plot(t,Q_art_beta(:,i),'Color',color_fade_blue(:,i),'LineWidth',5*exp(-i/7));
end
%plot(t,root_data(:,3),'-.c','LineWidth',3);
title('Arteries (\beta Branch)')
ylabel('Flow [mL/s]');
set(gca,'FontSize',20); 
leg=legend([h1,h2,h3,h4,h5],'Largest Vessel','','','','Smallest Vessel','location','northeast');
xlabel('Time [s]');
leg.FontSize = 15;
ylim([0 50])

figure;
for i=1:beta_max
hold on;
h1=plot(t,Q_ven_beta(:,1),'Color',color_fade_pink(:,1),'LineWidth',5*exp(-1/7));
h2=plot(t,Q_ven_beta(:,4),'Color',color_fade_pink(:,4),'LineWidth',5*exp(-4/7));
h3=plot(t,Q_ven_beta(:,7),'Color',color_fade_pink(:,7),'LineWidth',5*exp(-7/7));
h4=plot(t,Q_ven_beta(:,11),'Color',color_fade_pink(:,11),'LineWidth',5*exp(-11/7));
h5=plot(t,Q_ven_beta(:,beta_max),'Color',color_fade_pink(:,beta_max),'LineWidth',5*exp(-14/7));
plot(t,Q_ven_beta(:,i),'Color',color_fade_pink(:,i),'LineWidth',5*exp(-i/7));
end
%plot(t,root_data(:,4),'-.c','LineWidth',3);
title('Veins (\beta Branch)')
ylabel('Flow [mL/s]');
set(gca,'FontSize',20); 
leg=legend([h1,h2,h3,h4,h5],'Largest Vessel','','','','Smallest Vessel','location','northeast');
xlabel('Time [s]');
leg.FontSize = 15;
ylim([0 50])

figure; clf; 
plot([-r_alpha fliplr(r_alpha)],abs([mean(Q_art_alpha) fliplr(mean(Q_ven_alpha))]),...
    '-r.','LineWidth',3,'MarkerSize',36);
hold on;
plot([-r_beta fliplr(r_beta)],abs([mean(Q_art_beta) fliplr(mean(Q_ven_beta))]),...
    '-bo','LineWidth',3,'MarkerSize',10);
set(gca,'FontSize',20);
set(gca,'FontSize',20);
ylabel('Mean Flow (ml/s)');
xlabel('Radius (mm)');
legend('\alpha Branch','\beta Branch', 'location','north')
grid on
axis tight

%}
%% Shear Stress

figure; 
for i=1:alpha_max
 hold on;
 h1=plot(t,T_art_alpha(:,1),'Color',color_fade_lightblue(:,1),'LineWidth',5*exp(-1/7));
plot(t,T_art_alpha(:,i),'Color',color_fade_lightblue(:,i),'LineWidth',5*exp(-i/7));
end
set(gca,'FontSize',20); 
title('Arteries (\alpha)')
ylabel('\tau (dyn/mm^2)');
xlabel('Time (s)');
%ylim([-50 400])
grid on 

figure;
for i=2:alpha_max
hold on;
h2=plot(t,T_ven_alpha(:,1),'Color',color_fade_red(:,1),'LineWidth',5*exp(-1/7));
plot(t,T_ven_alpha(:,i),'Color',color_fade_red(:,i),'LineWidth',5*exp(-i/7));
end
set(gca,'FontSize',20); 
title('\alpha Branches')
legend([h1,h2],'Artery','Vein')
ylabel('\tau (dyn/mm^2)');
xlabel('Time (s)');
%ylim([-50 400])
grid on

figure;
for i=1:beta_max
hold on;
% h1=plot(t,T_art_beta(:,1),'Color',color_fade_blue(:,1),'LineWidth',5*exp(-1/7));
% h2=plot(t,T_ven_beta(:,1),'Color',color_fade_blue(:,4),'LineWidth',5*exp(-4/7));
% h3=plot(t,T_ven_beta(:,1),'Color',color_fade_blue(:,8),'LineWidth',5*exp(-8/7));
% h4=plot(t,T_ven_beta(:,1),'Color',color_fade_blue(:,12),'LineWidth',5*exp(-12/7));
% h5=plot(t,T_ven_beta(:,1),'Color',color_fade_blue(:,14),'LineWidth',5*exp(-14/7));
plot(t,T_art_beta(:,i),'Color',color_fade_blue(:,i),'LineWidth',5*exp(-i/7));
end
title('Arteries (\beta Branch)')
xlabel('Time (s)');
% legend([h1,h2,h3,h4,h5],'Largest','','','','Smallest')
ylabel('\tau (dyn/mm^2)');

figure;
for i=1:beta_max
hold on;
% h1=plot(t,T_ven_beta(:,1),'Color',color_fade_pink(:,1),'LineWidth',5*exp(-1/7));
% h2=plot(t,T_ven_beta(:,1),'Color',color_fade_pink(:,4),'LineWidth',5*exp(-4/7));
% h3=plot(t,T_ven_beta(:,1),'Color',color_fade_pink(:,8),'LineWidth',5*exp(-8/7));
% h4=plot(t,T_ven_beta(:,1),'Color',color_fade_pink(:,12),'LineWidth',5*exp(-12/7));
% h5=plot(t,T_ven_beta(:,1),'Color',color_fade_pink(:,14),'LineWidth',5*exp(-14/7));
plot(t,T_ven_beta(:,i),'Color',color_fade_pink(:,i),'LineWidth',5*exp(-i/7));
end
set(gca,'FontSize',20); 
title('Veins (\beta Branch)')
xlabel('Time (s)');
% legend([h1,h2,h3,h4,h5],'Largest','','','','Smallest')
ylabel('\tau (dyn/mm^2)');
%ylim([-50 400])
grid on


figure(23); clf; 
semilogx(r_alpha,mean(T_art_alpha),'-r.','LineWidth',3,'MarkerSize',36);
hold on; grid on;
semilogx(r_beta,mean(T_art_beta),'-bo','LineWidth',3,'MarkerSize',10);
set(gca,'FontSize',20);
ylabel('Mean \tau (dyn/mm^2)');
%set(gca,'XDir','reverse');
semilogx(r_alpha,mean(T_ven_alpha),'-r.','LineWidth',3,'MarkerSize',36);
hold on; grid on;
semilogx(r_beta,mean(T_ven_beta),'-bo','LineWidth',3,'MarkerSize',10);
set(gca,'FontSize',20);
xlabel('Radius (mm)');
legend('\alpha Branch','\beta Branch', 'location','northeast')
axis tight
grid on



%% Cyclic Stretch
%{
figure(30); clf; 
h1=plot(-r_alpha, CS_art_alpha,'-r','LineWidth',3,'MarkerSize',36);
hold on
plot(r_alpha,CS_ven_alpha, '-b','LineWidth',3,'MarkerSize',36);
%h2=plot(-r_beta,CS_art_beta, '-bo','LineWidth',3,'MarkerSize',10);
%plot(r_beta,CS_ven_beta, '-bo','LineWidth',3,'MarkerSize',10);
set(gca,'FontSize',20);
ylabel('CS (%)');
xlabel('r (mm)');
xlim([-0.05 0.05])
%ylim([0 250])
Ax = gca;
%set(gca, 'YTick', 0:10:100)
%set(gca, 'XTick', -0.005:0.001:0.005)
%legend([h1,h2],'\alpha Branch','\beta Branch', 'location','northwest')
set(gca,'fontsize',30)
grid on
%}
 save('save_1s.mat', 'P_art_alpha', 'P_art_beta', 'Q_art_alpha','Q_art_beta',...
    'CS_art_alpha', 'CS_art_beta', 'T_art_alpha', 'T_art_beta',...
    'P_ven_alpha', 'P_ven_beta', 'Q_ven_alpha','Q_ven_beta',...
    'CS_ven_alpha', 'CS_ven_beta', 'T_ven_alpha', 'T_ven_beta', ...
    'r_alpha','mu_alpha','r_beta','mu_beta','t')

end
