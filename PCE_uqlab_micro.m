%% PCE METAMODELING: MULTIPLE OUTPUTS
%
% This example showcases an application of polynomial chaos expansion
% (PCE) to the metamodeling of a simply supported beam model
% with multiple outputs.
% The model computes the deflections at several points along the length
% of the beam subjected to a uniform random load.

%% 1 - INITIALIZE UQLAB
%
% Clear all variables from the workspace, set the random number generator
% for reproducible results, and initialize the UQLab framework:
clearvars
rng(100,'twister')
uqlab

fname = {'UQ_micro_new/Art_Alpha_',...
         'UQ_micro_new/Art_Beta_',...
         'UQ_micro_new/Ven_Alpha_',...
         'UQ_micro_new/Ven_Beta_'};

% fname = {'UQ_micro_BC/Art_Alpha_',...
%          'UQ_micro_BC/Art_Beta_',...
%          'UQ_micro_BC/Ven_Alpha_',...
%          'UQ_micro_BC/Ven_Beta_'};
%% 2- Names
Names = {'f3_{A}','f3_{MV}','f3_{V}',...
    '\alpha','\beta','lrr_A','lrr_V','rm'};

% Names = {'f3_{A}','f3_{MV}','f3_{V}',...
%     '\alpha','\beta','lrr_A','lrr_V','rm',...
%          'QA0','QA1','QA2',...
%          'LA0','LA1','LA2','LA5',...
%          '\phi_1','\phi_2', '\phi_5'};

%% 3 - PROBABILISTIC INPUT MODEL
%
% The simply supported beam model has five inputs,
% modeled by independent lognormal random variables.
% The detailed model is given in the following table:
f_fload = 'micro/p_micro_';
counter = 1;
load new_parsonly\pq_TwoSided_parsonly_4_6_23.mat upp low par_sample chi

% load new_BC_run\pq_Terminal_BCs.mat chi upp low par_sample
MetaOpts.Display = 'quiet';
[num_samp,num_par] = size(par_sample);
% For BCs
% num_samp = 2381;

Input = [];
% Define an INPUT object with the following marginals:
for i=1:num_par
    Input.Marginals(i).Name = Names{i};  % beam width
    Input.Marginals(i).Type = 'Uniform';
    Input.Marginals(i).Parameters = [low(i) upp(i)];  % (m)
end
myInput = uq_createInput(Input);



% Data
pce_degree = 3;
N_train = min(nchoosek(num_par+pce_degree,pce_degree).*2,length(par_sample));
samps = 100:100+N_train-1;
num_samp = N_train;
% Data
par_sample = par_sample(samps,:);
X_train   = par_sample;


% Initialize arrays
num_pts_outs = 50;
QoI_P    = zeros(num_pts_outs,num_samp);
QoI_Q    = zeros(num_pts_outs,num_samp);
QoI_CS   = zeros(num_pts_outs,num_samp);
QoI_WSS  = zeros(num_pts_outs,num_samp);
x_interp = linspace(0,1,num_pts_outs);

%%
MetaOpts.Type = 'Metamodel';
MetaOpts.MetaType = 'PCE';
MetaOpts.Degree = 3;
MetaOpts.Method = 'OLS';
MetaOpts.ExpDesign.X = X_train;

PCESobol.Type = 'Sensitivity';
PCESobol.Method = 'Sobol';
PCESobol.Sobol.Order = 3;

for which_data = 1:8
    
    load(strcat(f_fload,num2str(which_data)));

    for which_ves = 1:4
        for k=1:num_samp
            n_ST_ves = size(P_ST{which_ves,samps(k)},2);
            QoI_P(:,k)  = interp1(linspace(0,1,n_ST_ves),mean(P_ST{which_ves,samps(k)}),x_interp);
            QoI_Q(:,k)  = interp1(linspace(0,1,n_ST_ves),mean(Q_ST{which_ves,samps(k)}),x_interp);
            QoI_CS(:,k) = interp1(linspace(0,1,n_ST_ves),CS_ST{which_ves,samps(k)},x_interp);
            QoI_WSS(:,k) = interp1(linspace(0,1,n_ST_ves),WSS_ST{which_ves,samps(k)},x_interp);
        end

        if which_ves>2
            QoI_P   = flipud(QoI_P);
            QoI_Q   = flipud(QoI_Q);
            QoI_CS  = flipud(QoI_CS);
            QoI_WSS = flipud(QoI_WSS);            
        end
        
        
        %%
        % Create the PCE metamodels:
        MetaOpts.ExpDesign.Y = QoI_P';
        PCE_P = uq_createModel(MetaOpts);
        Sobol_P = uq_createAnalysis(PCESobol);

        MetaOpts.ExpDesign.Y = QoI_Q';
        PCE_Q = uq_createModel(MetaOpts);
        Sobol_Q = uq_createAnalysis(PCESobol);

        MetaOpts.ExpDesign.Y = QoI_CS';
        PCE_CS = uq_createModel(MetaOpts);
        Sobol_CS = uq_createAnalysis(PCESobol);

        MetaOpts.ExpDesign.Y = QoI_WSS';
        PCE_WSS = uq_createModel(MetaOpts);
        Sobol_WSS = uq_createAnalysis(PCESobol);

        
        fname_save = strcat(fname{which_ves},num2str(which_data));
        save(fname_save,'PCE_P','PCE_Q','PCE_CS','PCE_WSS',...
              'Sobol_P','Sobol_Q','Sobol_CS','Sobol_WSS')
        
    end
end


return;
% %% Use bootstrapping?
% [YPCval,YPC_var,YPCval_Bootstrap] = uq_evalModel(myPCE_OLS,X_test);
% which_test = 50;
% uq_figure
% p(1) = uq_plot(Y_test(which_test,:), 'k');
% hold on
% cmap = get(gca,'ColorOrder');
% p(2) = uq_plot(YPCval(which_test,:));
% p(3) = uq_plot(...
%     squeeze(quantile(YPCval_Bootstrap(which_test,:,:),0.025,2)), '--',...
%     'Color', cmap(1,:));
% p(4) = uq_plot(...
%     squeeze(quantile(YPCval_Bootstrap(which_test,:,:),0.975,2)), '--',...
%     'Color', cmap(1,:));
% pb = uq_plot(...
%     squeeze(YPCval_Bootstrap(which_test,:,:))',...
%     'LineWidth', 1.5, 'Color', cmap(2,:));
% h = get(gca, 'Children');
% set(gca, 'Children', h(end:-1:1))  % Reorder the lines in the plot
% hold off
% % Add labels and a legend
% xlabel('$\mathrm{X}$')
% ylabel('$\mathcal{M}(X)$')
% uq_legend(...
%     [p([1:3]) pb(1)],...
%     {'True', 'PCE', '95\% Confidence Bounds', 'Replications'})
% 

%%
%% Now perform Sensitivity analysis
% PCESobol.Type = 'Sensitivity';
% PCESobol.Method = 'Sobol';
% PCESobol.Sobol.Order = 3;
% PCESobolAnalysis = uq_createAnalysis(PCESobol);
% 
% %% Now plot
% % Create the plot
% uq_figure('Name', 'Total Sobol'' Indices')
% barWidth = 1;
% uq_bar(1:num_par, PCESobolAnalysis.Results.Total(:,1), barWidth)
% % Set axes limits
% ylim([0 1])
% xlim([0 num_par])
% % Set labels
% xlabel('Variable name');
% xticklabels(Names)
% ylabel('Total Sobol'' indices')
