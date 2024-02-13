% This function will calculate the shear stress in the microcirculation
% using the alpha and beta branch data provided. The alpha and beta data
% structures contain arterial pressure, venous pressure, arterial flow,
% venous flow. THERE ARE NO OPTIONAL ARGUMENTS
% This code is a working file and comes with NO GUARANTEES.
%
% Authors: MU Qureshi, MJ Colebank, M Bartolo, MS Olufsen
%
% Last edited: 7/11/2023, MJC

function [WSS,CS,rads] = get_shear_CS(alpha_data,beta_data,pars)
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


%% Calculate the cyclic stretch in the arteries & veins
r_art_alpha = sqrt(A_art_alpha./pi); r_art_beta = sqrt(A_art_beta./pi);
r_ven_alpha = sqrt(A_ven_alpha./pi); r_ven_beta = sqrt(A_ven_beta./pi);

CS_art_alpha = 100*(max(r_art_alpha)-min(r_art_alpha))./min(r_art_alpha);
CS_ven_alpha = 100*(max(r_ven_alpha)-min(r_ven_alpha))./min(r_ven_alpha);

CS_art_beta = 100*(max(r_art_beta)-min(r_art_beta))./min(r_art_beta);
CS_ven_beta = 100*(max(r_ven_beta)-min(r_ven_beta))./min(r_ven_beta);

% Mean shear
% WSS_art_alpha = 4.*mu_alpha.*mean(Q_art_alpha)./(pi.*mean(r_art_alpha).^3);
% WSS_art_beta = 4.*mu_beta.*mean(Q_art_beta)./(pi.*mean(r_art_beta).^3);
% WSS_ven_alpha = 4.*mu_alpha.*mean(Q_ven_alpha)./(pi.*mean(r_ven_alpha).^3);
% WSS_ven_beta = 4.*mu_beta.*mean(Q_ven_beta)./(pi.*mean(r_ven_beta).^3);

WSS_art_alpha = 4.*mu_alpha.*(dt.*trapz(Q_art_alpha)./T)./(pi.*mean(r_art_alpha).^3);
WSS_art_beta = 4.*mu_beta.*(dt.*trapz(Q_art_beta)./T)./(pi.*mean(r_art_beta).^3);
WSS_ven_alpha = 4.*mu_alpha.*(dt.*trapz(Q_ven_alpha)./T)./(pi.*mean(r_ven_alpha).^3);
WSS_ven_beta = 4.*mu_beta.*(dt.*trapz(Q_ven_beta)./T)./(pi.*mean(r_ven_beta).^3);


rads.alpha = [r_art_alpha;r_ven_alpha];
rads.beta = [r_art_beta;r_ven_beta];
CS.alpha = [CS_art_alpha;CS_ven_alpha];
CS.beta =  [CS_art_beta;CS_ven_beta];
WSS.alpha = [WSS_art_alpha;WSS_ven_alpha];
WSS.beta =  [WSS_art_beta;WSS_ven_beta];
end
