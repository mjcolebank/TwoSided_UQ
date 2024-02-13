fname = {'UQ_micro_new/Art_Alpha_',...
    'UQ_micro_new/Art_Beta_',...
    'UQ_micro_new/Ven_Alpha_',...
    'UQ_micro_new/Ven_Beta_'};
fname2 = {'UQ_micro_new/figures/Art_Alpha_',...
    'UQ_micro_new/figures/Art_Beta_',...
    'UQ_micro_new/figures/Ven_Alpha_',...
    'UQ_micro_new/figures/Ven_Beta_'};
num_par = 8;
% With BCs
% fname = {'UQ_micro_BC/Art_Alpha_',...
%     'UQ_micro_BC/Art_Beta_',...
%     'UQ_micro_BC/Ven_Alpha_',...
%     'UQ_micro_BC/Ven_Beta_'};
% fname2 = {'UQ_micro_BC/figures/Art_Alpha_',...
%     'UQ_micro_BC/figures/Art_Beta_',...
%     'UQ_micro_BC/figures/Ven_Alpha_',...
%     'UQ_micro_BC/figures/Ven_Beta_'};
% num_par = 18;
%% 2- Names
% Only parameters in the system
Names = {'f3_{A}','f3_{MV}','f3_{V}',...
    '\alpha','\beta','lrr_A','lrr_V','rm'};%,...
% cmap = get_color_map(num_par);

% Include the BCs
% Names = {'f3_{A}','f3_{MV}','f3_{V}',...
%     '\alpha','\beta','lrr_A','lrr_V','rm',...
%          'QA0','QA1','QA2',...
%          'LA0','LA1','LA2','LA5',...
%          '\phi_1','\phi_2', '\phi_5'};
cmap = get_color_map(num_par);

num_pts_outs = 50;
x_interp = linspace(0,1,num_pts_outs);
% x_interp_A = linspace(0,1,num_pts_outs);
% x_interp_V = x_interp_A+0.1;
mu_P  = zeros(num_pts_outs,1); std_P = zeros(num_pts_outs,1);
mu_Q  = zeros(num_pts_outs,1); std_Q = zeros(num_pts_outs,1);
mu_CS  = zeros(num_pts_outs,1); std_CS = zeros(num_pts_outs,1);
mu_WSS  = zeros(num_pts_outs,1); std_WSS = zeros(num_pts_outs,1);
for which_data = 1:8
    for which_ves = 1:4
        fname_save = strcat(fname{which_ves},num2str(which_data));
        fname_plot = strcat(fname2{which_ves},num2str(which_data));
        load(fname_save);
%         if which_ves<3
%             x_interp = x_interp_A;
%         else
%             x_interp = x_interp_V;
%         end
        %% Get Mean & STD
        for i=1:num_pts_outs
            mu_P(i)   = PCE_P.PCE(i).Moments.Mean;
            mu_Q(i)   = PCE_Q.PCE(i).Moments.Mean;
            mu_CS(i)  = PCE_CS.PCE(i).Moments.Mean;
            mu_WSS(i) = PCE_WSS.PCE(i).Moments.Mean;

            std_P(i)   = sqrt(PCE_P.PCE(i).Moments.Var);
            std_Q(i)   = sqrt(PCE_Q.PCE(i).Moments.Var);
            std_CS(i)  = sqrt(PCE_CS.PCE(i).Moments.Var);
            std_WSS(i) = sqrt(PCE_WSS.PCE(i).Moments.Var);
        end

        % Pressure
        figure(10.*which_ves); clf; hold on;
        fill([x_interp'; fliplr(x_interp)'],[mu_P+std_P; flipud(mu_P-std_P)],[0.9 0.6 0.6],'EdgeAlpha',0,'FaceAlpha',0.4);
        plot(x_interp,mu_P,'k','LineWidth',3);
        set(gca,'FontSize',30); grid on;
        ylabel('Pressure (mmHg)');
        xticks([0 0.5 1]);
        ylim([0 60]);
        print(strcat(fname_plot,'_P_UQ'),'-dpng');

        % Flow
        figure(10.*which_ves+1); clf; hold on;
        fill([x_interp'; fliplr(x_interp)'],[mu_Q+std_Q; flipud(mu_Q-std_Q)],[0.9 0.6 0.6],'EdgeAlpha',0,'FaceAlpha',0.4);
        plot(x_interp,mu_Q,'k','LineWidth',3);
        set(gca,'FontSize',30); grid on;
        ylabel('Flow (mL/s)');
        xticks([0 0.5 1]);
        ylim([0 26])
        print(strcat(fname_plot,'_Q_UQ'),'-dpng');

        % CS
        figure(10.*which_ves+2); clf; hold on;
        fill([x_interp'; fliplr(x_interp)'],[mu_CS+std_CS; flipud(mu_CS-std_CS)],[0.9 0.6 0.6],'EdgeAlpha',0,'FaceAlpha',0.4);
        plot(x_interp,mu_CS,'k','LineWidth',3);
        set(gca,'FontSize',30); grid on;
        ylabel('Cyclic Stretch (%)');
        xticks([0 0.5 1]);
        ylim([0 8]);
        print(strcat(fname_plot,'_CS_UQ'),'-dpng');

        % WSS
        figure(10.*which_ves+3); clf; hold on;
        fill([x_interp'; fliplr(x_interp)'],[mu_WSS+std_WSS; flipud(mu_WSS-std_WSS)],[0.9 0.6 0.6],'EdgeAlpha',0,'FaceAlpha',0.4);
        plot(x_interp,mu_WSS,'k','LineWidth',3);
        set(gca,'FontSize',30); grid on;
        ylabel('WSS (dyne/cm^2)');
        xticks([0 0.5 1]);
        ylim([-25 175]);
        print(strcat(fname_plot,'_WSS_UQ'),'-dpng');

        %% Now perform Sensitivity analysis
        % Create the plot
        %         barWidth = 1;
        %         bar(1:num_par, PCESobolAnalysis.Results.Total(:,1), barWidth)
% 
%         figure(100.*which_ves); clf; hold on;
%         for jj=1:num_par
%             semilogy(x_interp,Sobol_P.Results.Total(jj,:),'LineWidth',3,'Color',cmap(jj,:));
%         end
%         ylim([0 1])
%         ylabel('Total Sobol'' indices');
%         grid on; title('Pressure')
%         set(gca,'FontSize',20);
%         print(strcat(fname_plot,'_P_ST'),'-dpng');
% 
%         figure(100.*which_ves+1); clf; hold on;
%         for jj=1:num_par
%             semilogy(x_interp,Sobol_Q.Results.Total(jj,:),'LineWidth',3,'Color',cmap(jj,:));
%         end
%         ylim([0 1])
%         ylabel('Total Sobol'' indices');
%         grid on; title('Flow')
%         set(gca,'FontSize',20);
%         print(strcat(fname_plot,'_Q_ST'),'-dpng');
% 
%         figure(100.*which_ves+2); clf; hold on;
%         for jj=1:num_par
%             semilogy(x_interp,Sobol_CS.Results.Total(jj,:),'LineWidth',3,'Color',cmap(jj,:));
%         end
%         ylim([0 1])
%         ylabel('Total Sobol'' indices');
%         grid on; title('Cyclic Stretch')
%         set(gca,'FontSize',20);
%         print(strcat(fname_plot,'_CS_ST'),'-dpng');
% 
%         figure(100.*which_ves+3); clf; hold on;
%         for jj=1:num_par
%             semilogy(x_interp,Sobol_WSS.Results.Total(jj,:),'LineWidth',3,'Color',cmap(jj,:));
%         end
%         ylim([0 1])
%         ylabel('Total Sobol'' indices');
%         grid on; title('WSS')
%         set(gca,'FontSize',20);
%         print(strcat(fname_plot,'_WSS_ST'),'-dpng');
%         %% First order
%         figure(1000.*which_ves); clf; hold on;
%         for jj=1:num_par
%             semilogy(x_interp,Sobol_P.Results.FirstOrder(jj,:),'LineWidth',3,'Color',cmap(jj,:));
%         end
%         ylim([0 1])
%         ylabel('First order Sobol'' indices');
%         grid on; title('Pressure')
%         set(gca,'FontSize',20);
%         print(strcat(fname_plot,'_P_Si'),'-dpng');
% 
%         figure(1000.*which_ves+1); clf; hold on;
%         for jj=1:num_par
%             semilogy(x_interp,Sobol_Q.Results.FirstOrder(jj,:),'LineWidth',3,'Color',cmap(jj,:));
%         end
%         ylim([0 1])
%         ylabel('First order Sobol'' indices');
%         grid on; title('Flow')
%         set(gca,'FontSize',20);
%         print(strcat(fname_plot,'_Q_Si'),'-dpng');
% 
%         figure(1000.*which_ves+2); clf; hold on;
%         for jj=1:num_par
%             semilogy(x_interp,Sobol_CS.Results.FirstOrder(jj,:),'LineWidth',3,'Color',cmap(jj,:));
%         end
%         ylim([0 1])
%         ylabel('First order Sobol'' indices');
%         grid on; title('Cyclic Stretch')
%         set(gca,'FontSize',20);
%         print(strcat(fname_plot,'_CS_Si'),'-dpng');
% 
%         figure(1000.*which_ves+3); clf; hold on;
%         for jj=1:num_par
%             semilogy(x_interp,Sobol_WSS.Results.FirstOrder(jj,:),'LineWidth',3,'Color',cmap(jj,:));
%         end
%         ylim([0 1])
%         ylabel('First order Sobol'' indices');
%         grid on; title('WSS')
%         set(gca,'FontSize',20);
%         print(strcat(fname_plot,'_WSS_Si'),'-dpng');
        %% Lastly, Calculate the generalized Sobol indices based on Alexanderian et al.
        %         intD = trapz(ti, Di);
        %         w = Di ./ intD;
        %         gS(n) = trapz(ti, w .* Si);
        GPi   = zeros(num_par,num_pts_outs);  GPT = zeros(num_par,num_pts_outs);
        GQi   = zeros(num_par,num_pts_outs);  GQT = zeros(num_par,num_pts_outs);
        GCSi  = zeros(num_par,num_pts_outs);  GCST = zeros(num_par,num_pts_outs);
        GWSSi = zeros(num_par,num_pts_outs);  GWSST = zeros(num_par,num_pts_outs);

        for i=2:50
            intV_P = trapz(std_P(1:i).^2);
            GPi(:,i)    = trapz(std_P(1:i).^2 .* Sobol_P.Results.FirstOrder(:,1:i)'./intV_P);
            GPT(:,i)    = trapz(std_P(1:i).^2 .* Sobol_P.Results.Total(:,1:i)'./intV_P);

            intV_Q = trapz(std_Q(1:i).^2);
            GQi(:,i)    = trapz(std_Q(1:i).^2 .* Sobol_Q.Results.FirstOrder(:,1:i)'./intV_Q);
            GQT(:,i)    = trapz(std_Q(1:i).^2 .* Sobol_Q.Results.Total(:,1:i)'./intV_Q);

            intV_CS = trapz(std_CS(1:i).^2);
            GCSi(:,i)    = trapz(std_CS(1:i).^2 .* Sobol_CS.Results.FirstOrder(:,1:i)'./intV_CS);
            GCST(:,i)    = trapz(std_CS(1:i).^2 .* Sobol_CS.Results.Total(:,1:i)'./intV_CS);

            intV_WSS = trapz(std_WSS(1:i).^2);
            GWSSi(:,i)    = trapz(std_WSS(1:i).^2 .* Sobol_WSS.Results.FirstOrder(:,1:i)'./intV_WSS);
            GWSST(:,i)    = trapz(std_WSS(1:i).^2 .* Sobol_WSS.Results.Total(:,1:i)'./intV_WSS);
        end

        %         %% Generalized Sobol
                figure(10000.*which_ves);
                plot(x_interp,GPT','LineWidth',3);
                ylim([0 1])
                ylabel('First order Sobol'' indices');
                grid on; title('Pressure')
                set(gca,'FontSize',20);
                legend(Names);
        %
        %         figure(10000.*which_ves+1);
        %         semilogy(x_interp,GQT','LineWidth',3);
        %         ylim([0 1])
        %         ylabel('First order Sobol'' indices');
        %         grid on; title('Flow')
        %         set(gca,'FontSize',20);
        %         legend(Names);
        %
        %         figure(10000.*which_ves+2);
        %         semilogy(x_interp,GCST','LineWidth',3);
        %         ylim([0 1])
        %         ylabel('First order Sobol'' indices');
        %         grid on; title('Cyclic Stretch')
        %         set(gca,'FontSize',20);
        %         legend(Names);
        %
        %         figure(10000.*which_ves+3);
        %         semilogy(x_interp,GWSST','LineWidth',3);
        %         ylim([0 1])
        %         ylabel('First order Sobol'' indices');
        %         grid on; title('WSS')
        %         set(gca,'FontSize',20);
        %         legend(Names);



        %         figure(999); hold on;



    end
    pause(2); close all;

end
