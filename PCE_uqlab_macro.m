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
% clearvars
rng(100,'twister')
uqlab

fname = 'UQ_macro_new/ves_smooth_final_';
% fname = 'UQ_testmacro/pce3'
% load new_parsonly\pq_TwoSided_parsonly_4_6_23.mat

% fname = 'UQ_macro_BC/ves_';
% load new_BC_run\pq_newest_BCs.mat

% Radius Values

rvals = [1.35 0.9 1.1 0.842 0.481 0.922 0.755 ...
        0.757 0.514 0.433 0.293 0.829 0.562 0.460 0.610 ...
        0.641 0.716 0.864 0.824 ...
        0.757 0.514 0.433 0.293 0.829 0.562 0.460 0.610];
vel_power = 9;
visc = 0.032; % 1 g/cm/s is 1 dyne s /cm^2
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
MetaOpts.Display = 'quiet';
[num_samp,num_par] = size(par_sample);

Input = [];
% Define an INPUT object with the following marginals:
for i=1:num_par
    Input.Marginals(i).Name = Names{i};  % beam width
    Input.Marginals(i).Type = 'Uniform';
    Input.Marginals(i).Parameters = [low(i) upp(i)];  % (m)
end
myInput = uq_createInput(Input);



% Initialize arrays
subsamp = 4;
num_pts_outs = 512/subsamp;
QoI_P    = zeros(num_pts_outs,num_samp);
QoI_Q    = zeros(num_pts_outs,num_samp);
QoI_CS   = zeros(num_pts_outs,num_samp);
QoI_WSS  = zeros(num_pts_outs,num_samp);
x_interp = linspace(0,1,num_pts_outs);

%%
pce_deg = 4;
N_train = 1900;%min(nchoosek(num_par+pce_deg,pce_deg).*2,length(par_sample));
% samps = 100:100+N_train-1;
samps = 1:N_train;
% Data
par_sample = par_sample(samps,:);
X_train   = par_sample;

MetaOpts.Type = 'Metamodel';
MetaOpts.MetaType = 'PCE';
MetaOpts.Degree = pce_deg;
MetaOpts.Method = 'OLS';
MetaOpts.ExpDesign.X = X_train;

PCESobol.Type = 'Sensitivity';
PCESobol.Method = 'Sobol';
PCESobol.Sobol.Order = 3;
t = linspace(0,0.85,num_pts_outs);

%%
for which_ves = 1:27%[1 2 3 16 17 18 19]
        %Look at the QoI
        
    QoI_P = squeeze(p_PCE(1:subsamp:end,which_ves,samps));
    QoI_Q = squeeze(q_PCE(1:subsamp:end,which_ves,samps));
    par_stiff = par_sample(:,1);

    if which_ves>15
        QoI_Q = -QoI_Q;
        par_stiff = par_sample(:,2);
    end
    QoI_A = squeeze(A_PCE(1:subsamp:end,which_ves,samps));
    QoI_r = sqrt(QoI_A./pi);
    QoI_CS = (max(QoI_r)-min(QoI_r))./min(QoI_r);
    QoI_U = QoI_Q./QoI_A;
    QoI_WSS = visc.*(vel_power+2).* QoI_U./QoI_r;

    %% Average first few points of pressure
%     QoI_P(1,:) = QoI_P(5,:);
%     QoI_P(2,:) = QoI_P(5,:);
%     QoI_P(3,:) = QoI_P(5,:);
%     QoI_P(4,:) = QoI_P(5,:);
% QoI_P = movmean(QoI_P,5,1);

%     [fwd_comp, fwd_exp, bwd_comp, bwd_exp] = WIA(t,QoI_P,QoI_Q,QoI_A,rvals(which_ves),par_stiff',1);

    % Two cycles
%     P_WIA = [QoI_P;QoI_P]; t_WIA = [t t]'; 
%     Q_WIA = [QoI_Q; QoI_Q]; A_WIA = [QoI_A; QoI_A];
    % Less temporal resolution
    if subsamp==1
        which_pts = 1:512;%6:512;
    else
        which_pts = 1:num_pts_outs;
    end
    P_WIA = QoI_P(which_pts,:); 
    Q_WIA = QoI_Q(which_pts,:);
    t_WIA = t(which_pts); 
    A_WIA = QoI_A(which_pts,:);

    [fwd_comp, fwd_exp, bwd_comp, bwd_exp] = WIA(t_WIA,P_WIA,Q_WIA,A_WIA,rvals(which_ves),par_stiff',1);

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

    %% WIA PCEs
    MetaOpts.ExpDesign.Y = fwd_comp';
    PCE_fwdC = uq_createModel(MetaOpts);
    Sobol_fwdC = uq_createAnalysis(PCESobol);

    MetaOpts.ExpDesign.Y = fwd_exp';
    PCE_fwdE = uq_createModel(MetaOpts);
    Sobol_fwdE = uq_createAnalysis(PCESobol);

    MetaOpts.ExpDesign.Y = bwd_comp';
    PCE_bwdC = uq_createModel(MetaOpts);
    Sobol_bwdC = uq_createAnalysis(PCESobol);

    MetaOpts.ExpDesign.Y = bwd_exp';
    PCE_bwdE = uq_createModel(MetaOpts);
    Sobol_bwdE = uq_createAnalysis(PCESobol);

            fname_save = strcat(fname,num2str(which_ves),'_128subsamp');
            save(fname_save,'PCE_P','PCE_Q','PCE_CS','PCE_WSS',...
                'PCE_fwdC','PCE_fwdE','PCE_bwdC','PCE_bwdE',...
                  'Sobol_P','Sobol_Q','Sobol_CS','Sobol_WSS',...
                  'Sobol_fwdC','Sobol_fwdE','Sobol_bwdC','Sobol_bwdE')

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
