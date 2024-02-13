%% Try to calculate the correlation coefficient for the differenct PCE indicies

close all; clc; clear;
%%
load new_parsonly\test_data.mat
which_ves = 17;
vel_power = 9;
visc = 0.032; % 1 g/cm/s is 1 dyne s /cm^2


p_data = squeeze(p_test_data(1:8:end,which_ves,:))';
if which_ves>15
    q_data = -squeeze(q_test_data(1:8:end,which_ves,:))';
else
    q_data = squeeze(q_test_data(1:8:end,which_ves,:))';
end
A_data = squeeze(A_test_data(1:8:end,which_ves,:))';
r_data = sqrt(A_data./pi);
CS_data = (max(r_data,[],2)-min(r_data,[],2))./min(r_data,[],2);
U_data = q_data./A_data;
WSS_data = visc.*(vel_power+2).* U_data./r_data;

%%
load UQ_macro_new\ves_pce4_17_subsamp.mat

res_P = (p_data-uq_eval_uq_metamodel(PCE_P,par_test));
res_Q = (q_data-uq_eval_uq_metamodel(PCE_Q,par_test));
res_CS = (CS_data-uq_eval_uq_metamodel(PCE_CS,par_test));
res_WSS = (WSS_data-uq_eval_uq_metamodel(PCE_WSS,par_test));

R2_P = 1.0 - sum(res_P.^2)./(100.*var(p_data));
R2_Q = 1.0 - sum(res_Q.^2)./(100.*var(q_data));
R2_CS = 1.0 - sum(res_CS.^2)./(100.*var(CS_data));
R2_WSS = 1.0 - sum(res_WSS.^2)./(100.*var(WSS_data));

RMS_P = sum(sqrt(((res_P'./p_data').^2)));