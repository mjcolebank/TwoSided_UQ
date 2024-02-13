close all; clear all; clc;

addpath '/Users/michellebartolo/Documents/NCSU/Research/ST_2Sided_MJC'

%Large vessels
save_normal =load('save_normal','-mat');
save2 =load('save_halfgen','-mat');
save3 =load('save_rmin','-mat');
save4 =load('save_rmin2','-mat');
% save5 =load('save8','-mat');


%Small Vessels
save_normal_s =load('save_normal_s','-mat');
save2_s =load('save_half_s','-mat');
save3_s =load('save_rmin_s','-mat');
save4_s =load('save_rmin_2s','-mat');
% save5_s =load('save8_s','-mat');


%MPA
%Flow
figure;
subplot(1,3,1)
plot(save_normal.model.q_art(:,1),'k','linewidth',3)
hold on
plot(save2.model.q_art(:,1),'r','linewidth',3)
plot(save3.model.q_art(:,1),'b','linewidth',3)
plot(save4.model.q_art(:,1),'m','linewidth',3)
% plot(save5.model.q_art(:,1),'g','linewidth',3)

title('Flow - MPA')
set(gca,'fontsize',20)

%Pressure
subplot(1,3,2)
plot(save_normal.model.p_art(:,1),'k','linewidth',3)
hold on
plot(save2.model.p_art(:,1),'r','linewidth',3)
plot(save3.model.p_art(:,1),'b','linewidth',3)
plot(save4.model.p_art(:,1),'m','linewidth',3)
% plot(save5.model.p_art(:,1),'g','linewidth',3)
title('Pressure - MPA')
set(gca,'fontsize',20)

%Shear Stress
subplot(1,3,3)
plot(save_normal.shear.tau_art_exp(:,1),'k','linewidth',3)
plot(save2.shear.tau_art_exp(:,1),'r','linewidth',3)
plot(save3.shear.tau_art_exp(:,1),'b','linewidth',3)
plot(save4.shear.tau_art_exp(:,1),'m','linewidth',3)
% plot(save5.shear.tau_art_exp(:,1),'g','linewidth',3)
title('Shear Stress - MPA')
set(gca,'fontsize',20)


%RPA
%Flow
figure;
subplot(1,3,1)
plot(save_normal.model.q_art(:,2),'k','linewidth',3)
hold on
plot(save2.model.q_art(:,2),'r','linewidth',3)
plot(save3.model.q_art(:,2),'b','linewidth',3)
plot(save4.model.q_art(:,2),'m','linewidth',3)
% plot(save5.model.q_art(:,2),'g','linewidth',3)

title('Flow - RPA')
set(gca,'fontsize',20)

%Pressure
subplot(1,3,2)
plot(save_normal.model.p_art(:,2),'k','linewidth',3)
hold on
plot(save2.model.p_art(:,2),'r','linewidth',3)
plot(save3.model.p_art(:,2),'b','linewidth',3)
plot(save4.model.p_art(:,2),'m','linewidth',3)
% plot(save5.model.p_art(:,2),'g','linewidth',3)
title('Pressure - RPA')
set(gca,'fontsize',20)

%Shear Stress
subplot(1,3,3)
plot(save_normal.shear.tau_art_exp(:,2),'k','linewidth',3)
hold on
plot(save2.shear.tau_art_exp(:,2),'r','linewidth',3)
plot(save3.shear.tau_art_exp(:,2),'b','linewidth',3)
plot(save4.shear.tau_art_exp(:,2),'m','linewidth',3)
% plot(save5.shear.tau_art_exp(:,2),'g','linewidth',3)
title('Shear Stress - RPA')
set(gca,'fontsize',20)




%LPA
%Flow
figure;
subplot(1,3,1)
plot(save_normal.model.q_art(:,3),'k','linewidth',3)
hold on
plot(save2.model.q_art(:,3),'r','linewidth',3)
plot(save3.model.q_art(:,3),'b','linewidth',3)
plot(save4.model.q_art(:,3),'m','linewidth',3)
% plot(save5.model.q_art(:,3),'g','linewidth',3)
title('Flow - LPA')
set(gca,'fontsize',20)

%Pressure
subplot(1,3,2)
plot(save_normal.model.p_art(:,3),'k','linewidth',3)
hold on
plot(save2.model.p_art(:,3),'r','linewidth',3)
plot(save3.model.p_art(:,3),'b','linewidth',3)
plot(save4.model.p_art(:,3),'m','linewidth',3)
% plot(save5.model.p_art(:,3),'g','linewidth',3)
title('Pressure - LPA')
set(gca,'fontsize',20)

%Shear Stress
subplot(1,3,3)
plot(save_normal.shear.tau_art_exp(:,3),'k','linewidth',3)
hold on
plot(save2.shear.tau_art_exp(:,3),'r','linewidth',3)
plot(save3.shear.tau_art_exp(:,3),'b','linewidth',3)
plot(save4.shear.tau_art_exp(:,3),'m','linewidth',3)
% plot(save5.shear.tau_art_exp(:,3),'g','linewidth',3)
title('Shear Stress - LPA')
set(gca,'fontsize',20)



%Large Veins
%Flow
figure;
subplot(1,3,1)
plot(save_normal.model.q_ven(:,1),'k','linewidth',3)
hold on
plot(save2.model.q_ven(:,1),'r','linewidth',3)
plot(save3.model.q_ven(:,1),'b','linewidth',3)
plot(save4.model.q_ven(:,1),'m','linewidth',3)
% plot(save5.model.q_ven(:,1),'g','linewidth',3)
title('Flow - RIV')
set(gca,'fontsize',20)

%Pressure
subplot(1,3,2)
plot(save_normal.model.p_ven(:,1),'k','linewidth',3)
hold on
plot(save2.model.p_ven(:,1),'r','linewidth',3)
plot(save3.model.p_ven(:,1),'b','linewidth',3)
plot(save4.model.p_ven(:,1),'m','linewidth',3)
% plot(save5.model.p_ven(:,1),'g','linewidth',3)
title('Pressure - RIV')
set(gca,'fontsize',20)

%Shear Stress
subplot(1,3,3)
plot(save_normal.shear.tau_ven_exp(:,1),'k','linewidth',3)
hold on
plot(save2.shear.tau_ven_exp(:,1),'r','linewidth',3)
plot(save3.shear.tau_ven_exp(:,1),'b','linewidth',3)
plot(save4.shear.tau_ven_exp(:,1),'m','linewidth',3)
% plot(save5.shear.tau_ven_exp(:,1),'g','linewidth',3)
title('Shear Stress - RIV')
set(gca,'fontsize',20)





%Small arteries/veins
%Flow
figure;

subplot(1,3,1)
plot(-save_normal_s.r_alpha,mean(save_normal_s.Q_art_alpha),'k','linewidth',3)
hold on
plot(-save2_s.r_alpha,mean(save2_s.Q_art_alpha),'r','linewidth',3)
plot(-save3_s.r_alpha,mean(save3_s.Q_art_alpha),'b','linewidth',3)
plot(-save4_s.r_alpha,mean(save4_s.Q_art_alpha),'m','linewidth',3)
% plot(-save5_s.r_alpha,mean(save5_s.Q_art_alpha),'g','linewidth',3)

plot(save_normal_s.r_alpha,mean(save_normal_s.Q_ven_alpha),'k','linewidth',3)
plot(save2_s.r_alpha,mean(save2_s.Q_ven_alpha),'r','linewidth',3)
plot(save3_s.r_alpha,mean(save3_s.Q_ven_alpha),'b','linewidth',3)
plot(save4_s.r_alpha,mean(save4_s.Q_ven_alpha),'m','linewidth',3)
% plot(save5_s.r_alpha,mean(save5_s.Q_ven_alpha),'g','linewidth',3)

title('Flow')
set(gca,'fontsize',20)

%Pressure
subplot(1,3,2)
plot(-save_normal_s.r_alpha,mean(save_normal_s.P_art_alpha),'k','linewidth',3)
hold on
plot(-save2_s.r_alpha,mean(save2_s.P_art_alpha),'r','linewidth',3)
plot(-save3_s.r_alpha,mean(save3_s.P_art_alpha),'b','linewidth',3)
plot(-save4_s.r_alpha,mean(save4_s.P_art_alpha),'m','linewidth',3)
% plot(-save5_s.r_alpha,mean(save5_s.P_art_alpha),'g','linewidth',3)

plot(save_normal_s.r_alpha,mean(save_normal_s.P_ven_alpha),'k','linewidth',3)
plot(save2_s.r_alpha,mean(save2_s.P_ven_alpha),'r','linewidth',3)
plot(save3_s.r_alpha,mean(save3_s.P_ven_alpha),'b','linewidth',3)
plot(save4_s.r_alpha,mean(save4_s.P_ven_alpha),'m','linewidth',3)
% plot(save5_s.r_alpha,mean(save5_s.P_ven_alpha),'g','linewidth',3)

title('Pressure')
set(gca,'fontsize',20)

%Shear Stress
subplot(1,3,3)
plot(-save_normal_s.r_alpha,mean(save_normal_s.T_art_alpha),'k','linewidth',3)
hold on
plot(-save2_s.r_alpha,mean(save2_s.T_art_alpha),'r','linewidth',3)
plot(-save3_s.r_alpha,mean(save3_s.T_art_alpha),'b','linewidth',3)
plot(-save4_s.r_alpha,mean(save4_s.T_art_alpha),'m','linewidth',3)
% plot(-save5_s.r_alpha,mean(save5_s.T_art_alpha),'g','linewidth',3)

plot(save_normal_s.r_alpha,mean(save_normal_s.T_ven_alpha),'k','linewidth',3)
plot(save2_s.r_alpha,mean(save2_s.T_ven_alpha),'r','linewidth',3)
plot(save3_s.r_alpha,mean(save3_s.T_ven_alpha),'b','linewidth',3)
plot(save4_s.r_alpha,mean(save4_s.T_ven_alpha),'m','linewidth',3)
% plot(save5_s.r_alpha,mean(save5_s.T_ven_alpha),'g','linewidth',3)

title('Shear Stress')
set(gca,'fontsize',20)


%Shear Stress
figure;
sgt=sgtitle('Shear Stress over Time');
sgt.FontSize = 20;

subplot(1,2,1)
plot(save_normal_s.t,save_normal_s.T_art_alpha,'k','linewidth',3)
hold on
plot(save2_s.t,save2_s.T_art_alpha,'r','linewidth',3)
plot(save3_s.t,save3_s.T_art_alpha,'b','linewidth',3)
plot(save4_s.t,save4_s.T_art_alpha,'b','linewidth',3)
title('Arteries')
set(gca,'fontsize',18)

subplot(1,2,2)
plot(save_normal_s.t,save_normal_s.T_ven_alpha,'k','linewidth',3)
hold on
plot(save2_s.t,save2_s.T_ven_alpha,'r','linewidth',3)
plot(save3_s.t,save3_s.T_ven_alpha,'b','linewidth',3)
plot(save4_s.t,save4_s.T_ven_alpha,'m','linewidth',3)
title('Veins')
set(gca,'fontsize',18)

%Cyclic Stretch
figure;
plot(-save_normal_s.r_alpha,save_normal_s.CS_art_alpha,'k','linewidth',3)
hold on
plot(-save2_s.r_alpha,save2_s.CS_art_alpha,'r','linewidth',3)
plot(-save3_s.r_alpha,save3_s.CS_art_alpha,'b','linewidth',3)
plot(-save4_s.r_alpha,save4_s.CS_art_alpha,'m','linewidth',3)
% plot(-save5_s.r_alpha,save5_s.CS_art_alpha,'g','linewidth',3)

plot(save_normal_s.r_alpha,save_normal_s.CS_ven_alpha,'k','linewidth',3)
plot(save2_s.r_alpha,save2_s.CS_ven_alpha,'r','linewidth',3)
plot(save3_s.r_alpha,save3_s.CS_ven_alpha,'b','linewidth',3)
plot(save4_s.r_alpha,save4_s.CS_ven_alpha,'m','linewidth',3)
% plot(save5_s.r_alpha,save5_s.CS_ven_alpha,'g','linewidth',3)

title('Cyclic Stress')
set(gca,'fontsize',20)



figure;
plot(save_normal_s.t,save_normal_s.Q_art_alpha,'k')
hold on
plot(save2_s.t,save2_s.Q_art_alpha,'r')
plot(save3_s.t,save3_s.Q_art_alpha,'b')
plot(save4_s.t,save4_s.Q_art_alpha,'b')