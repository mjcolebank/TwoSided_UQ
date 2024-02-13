% This is a driver file that will run the fluids model from c++ and
% fortran by passing parameter values needed by the model. This code is a
% working file and comes with NO GUARANTEES.
%
% Authors: MU Qureshi, MJ Colebank, M Bartolo, MS Olufsen
%
% Last edited: 2/14/2020, MJC
%%
% clear; clc; close all;
%% Ensure the make file is compiled
! chmod +x treebranch
! make clean
! make veryclean
! make
if exist('alpha_beta')~=7
    mkdir('alpha_beta');
end
%

% load micro_test.mat
% load ../new_parsonly/pq_TwoSided_parsonly_terminal.mat p_PCE q_PCE par_sample
load ../new_BC_run/pq_Terminal_BCs.mat p_PCE q_PCE par_sample
%%
par_test = par_sample;
p_test = p_PCE;
q_test = q_PCE;
r_vals = [0.757, 0.514,0.433, 0.293, 0.829, 0.562, 0.460, 0.610];
num_term = length(r_vals);
num_samp = length(par_test);


for which_ves=2:num_term
    P_ST   = cell(4,num_samp);
    Q_ST   = cell(4,num_samp);
    CS_ST  = cell(4,num_samp);
    WSS_ST = cell(4,num_samp);
    r_ST   = cell(4,num_samp);
    for which_samp = 500%1:num_samp
        %% Decide which artery and vein to include for the analysis
        ves_id = which_ves;

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
        alpha_max = find(alpha_vec<r_min,1);%-1;
        beta_max  = find(beta_vec<r_min,1);%-1;

        %%
        pars = [period mu rho fa1 fa2 fa3 fv1 fv2 fv3 r_root ...
            r_min  lrrA lrrV alpha_b beta_b maxgen tmstps alpha_max beta_max ves_id]';

        % Write to file
        dlmwrite('parameters_MJC.dat',pars);

        %% Find the file in the previous folder
        qp_art = [qA_term pA_term];
        qp_ven = [qV_term pV_term];

        %% save
        dlmwrite('qp_art.dat',qp_art);
        dlmwrite('qp_ven.dat',qp_ven);

        %% run the fortran file
        flag = unix(sprintf('treebranch.exe'));

        %% All files are printed to the folder "alpha_beta". Read in those results
        % Now loop until you reach alpha or beta max
        P_art_alpha = zeros(tmstps,alpha_max); P_ven_alpha = P_art_alpha;
        Q_art_alpha = P_art_alpha; Q_ven_alpha = P_art_alpha;

        P_art_beta = zeros(tmstps,beta_max); P_ven_beta = P_art_beta;
        Q_art_beta = P_art_beta; Q_ven_beta = P_art_beta;
        alpha_data = zeros(tmstps,alpha_max,4); 
        beta_data  = zeros(tmstps,beta_max,4);
        cd('alpha_beta');
        for i=1:alpha_max
            name = strcat('alpha_',num2str(i-1),'.2d');
            data = dlmread(name);
            P_art_alpha(:,i) = data(:,1); P_ven_alpha(:,i) = data(:,2);
            Q_art_alpha(:,i) = data(:,3); Q_ven_alpha(:,i) = -1.*data(:,4);
            alpha_data(:,i,:) = data;

        end
        for i=1:beta_max
            name = strcat('beta_',num2str(i-1),'.2d');
            data = dlmread(name);
            P_art_beta(:,i) = data(:,1); P_ven_beta(:,i) = data(:,2);
            Q_art_beta(:,i) = data(:,3); Q_ven_beta(:,i) = -1.*data(:,4);
            beta_data(:,i,:) = data;
        end
        cd('../');

        P_ST{1,which_samp} = P_art_alpha; 
        P_ST{2,which_samp} = P_art_beta;
        P_ST{3,which_samp} = P_ven_alpha; 
        P_ST{4,which_samp} = P_ven_beta;

        Q_ST{1,which_samp} = Q_art_alpha; 
        Q_ST{2,which_samp} = Q_art_beta;
        Q_ST{3,which_samp} = Q_ven_alpha; 
        Q_ST{4,which_samp} = Q_ven_beta;
        
        pars_shear = [alpha beta period mu rho fa1 fa2 fa3 fv1 fv2 fv3 r_root ...
        r_min alpha_b beta_b lrrA lrrV maxgen new_tmstps]';
        [WSS,CS,rads] = get_shear_CS(alpha_data,beta_data,pars_shear);

        CS_ST{1,which_samp} = CS.alpha(1,:); 
        CS_ST{2,which_samp} = CS.beta(1,:);
        CS_ST{3,which_samp} = CS.alpha(2,:);
        CS_ST{4,which_samp} = CS.beta(2,:);

        WSS_ST{1,which_samp} = WSS.alpha(1,:); 
        WSS_ST{2,which_samp} = WSS.beta(1,:);
        WSS_ST{3,which_samp} = WSS.alpha(2,:);
        WSS_ST{4,which_samp} = WSS.beta(2,:);

        r_ST{1,which_samp} = rads.alpha(1,:); 
        r_ST{2,which_samp} = rads.beta(1,:);
        r_ST{3,which_samp} = rads.alpha(2,:);
        r_ST{4,which_samp} = rads.beta(2,:);
    end
    fname = strcat('../UQ_micro_BC/p_micro_BC_',num2str(which_ves));
    save(fname,'P_ST','Q_ST','WSS_ST','CS_ST','r_ST');
end

% %% Now save all the data
% % 2300 --> just parameters
% % 2381 --> BCs included as well
% n_vals = 2381; %2300
% for i=1%:8
%     P_ST = cell(4,n_vals);
%     Q_ST = cell(4,n_vals);
%     WSS_ST = cell(4,n_vals);
%     CS_ST  = cell(4,n_vals);
%     r_ST   = cell(4,n_vals);
%     for j=1:n_vals
%         P_ST{1,j} = P_ST{1,j,i}; P_ST{2,j} = P_ST{2,j,i};
%         P_ST{3,j} = P_ST{3,j,i}; P_ST{4,j} = P_ST{4,j,i};
% 
%         Q_ST{1,j} = Q_ST{1,j,i}; Q_ST{2,j} = Q_ST{2,j,i};
%         Q_ST{3,j} = Q_ST{3,j,i}; Q_ST{4,j} = Q_ST{4,j,i};
% 
%         WSS_ST{1,j} = WSS_ST{1,j,i}; WSS_ST{2,j} = WSS_ST{2,j,i};
%         WSS_ST{3,j} = WSS_ST{3,j,i}; WSS_ST{4,j} = WSS_ST{4,j,i};
% 
%         CS_ST{1,j} = CS_ST{1,j,i}; CS_ST{2,j} = CS_ST{2,j,i};
%         CS_ST{3,j} = CS_ST{3,j,i}; CS_ST{4,j} = CS_ST{4,j,i};
% 
%         r_ST{1,j} = r_ST{1,j,i}; r_ST{2,j} = r_ST{2,j,i};
%         r_ST{3,j} = r_ST{3,j,i}; r_ST{4,j} = r_ST{4,j,i};
%     end
% end
