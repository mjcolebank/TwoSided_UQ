fname = 'UQ_macro_new/ves_smooth_final_';
num_par = 8;
Names = {'f3_{A}','f3_{MV}','f3_{V}',...
    '\alpha','\beta','lrr_A','lrr_V','rm'};%,...
% nsamp = 331;%2300;
nsamp = 1900%990;
% fname = 'UQ_macro_BC/ves_';
% num_par = 18;
% Names = {'f3_{A}','f3_{MV}','f3_{V}',...
%     '\alpha','\beta','lrr_A','lrr_V','rm',...
%          'QA0','QA1','QA2',...
%          'LA0','LA1','LA2','LA5',...
%          '\phi_1','\phi_2', '\phi_5'};
% num_samp = 2381;
cmap = get_color_map(num_par);

num_pts_outs = 512/4;
num_pts_WIA = num_pts_outs-1;
t = linspace(0,0.85,num_pts_outs);
% Full output
% tW = linspace(0,0.85,num_pts_outs-6);
% For BC
tW = linspace(0,0.85,num_pts_outs-1);
num_WIA = length(tW);


mu_P  = zeros(num_pts_outs,1); std_P = zeros(num_pts_outs,1);
mu_Q  = zeros(num_pts_outs,1); std_Q = zeros(num_pts_outs,1);
mu_CS  = zeros(num_pts_outs,1); std_CS = zeros(num_pts_outs,1);
mu_WSS  = zeros(num_pts_outs,1); std_WSS = zeros(num_pts_outs,1);

mu_fwdC  = zeros(num_WIA,1); std_fwdC = zeros(num_WIA,1);
mu_fwdE  = zeros(num_WIA,1); std_fwdE = zeros(num_WIA,1);
mu_bwdC  = zeros(num_WIA,1); std_bwdC = zeros(num_WIA,1);
mu_bwdE  = zeros(num_WIA,1); std_bwdE = zeros(num_WIA,1);

%%
    for which_ves = [18]%[1 2 3 16 17 18 19]
        if which_ves<16
            shade_color = [0.6 0.6 0.9];
        else
            shade_color = [0.9 0.6 0.6];
        end
%         fname_save = strcat(fname,num2str(which_ves));
        % BC
%         fname_save = strcat(fname,num2str(which_ves),'_subsamp');
        fname_save = strcat(fname,num2str(which_ves),'_128subsamp');
        load(fname_save);
        %% Get Mean & STD
        for i=1:num_pts_outs
            mu_P(i)   = PCE_P.PCE(i).Moments.Mean;
            mu_Q(i)   = PCE_Q.PCE(i).Moments.Mean;
%             mu_CS(i)  = PCE_CS.PCE(i).Moments.Mean;
            mu_WSS(i) = PCE_WSS.PCE(i).Moments.Mean;

            std_P(i)   = sqrt(PCE_P.PCE(i).Moments.Var);
            std_Q(i)   = sqrt(PCE_Q.PCE(i).Moments.Var);
%             std_CS(i)  = sqrt(PCE_CS.PCE(i).Moments.Var);
            std_WSS(i) = sqrt(PCE_WSS.PCE(i).Moments.Var);
            
        end
        for i=1:length(tW)%num_pts_outs-1
            mu_fwdC(i) = PCE_fwdC.PCE(i).Moments.Mean;
            mu_fwdE(i) = PCE_fwdE.PCE(i).Moments.Mean;
            mu_bwdC(i) = PCE_bwdC.PCE(i).Moments.Mean;
            mu_bwdE(i) = PCE_bwdE.PCE(i).Moments.Mean;

            std_fwdC(i)   = sqrt(PCE_fwdC.PCE(i).Moments.Var);
            std_fwdE(i)   = sqrt(PCE_fwdE.PCE(i).Moments.Var);
            std_bwdC(i)  = sqrt(PCE_bwdC.PCE(i).Moments.Var);
            std_bwdE(i) = sqrt(PCE_bwdE.PCE(i).Moments.Var);
            
        end

        % Pressure
        figure(10.*which_ves); clf; hold on;
        fill([t'; fliplr(t)'],[mu_P+std_P; flipud(mu_P-std_P)],shade_color,'EdgeAlpha',0,'FaceAlpha',0.4);
        plot(t,mu_P,'k','LineWidth',3);
        set(gca,'FontSize',30); grid on;
        ylabel('Pressure (mmHg)');
        xticks([0 0.5 1]);
        xlim([t(1) t(end)])
        print(strcat(fname_save,'_P_UQ'),'-dpng');

        % Flow
        figure(10.*which_ves+1); clf; hold on;
        fill([t'; fliplr(t)'],[mu_Q+std_Q; flipud(mu_Q-std_Q)],shade_color,'EdgeAlpha',0,'FaceAlpha',0.4);
        plot(t,mu_Q,'k','LineWidth',3);
        set(gca,'FontSize',30); grid on;
        ylabel('Flow (mL/s)');
        xticks([0 0.5 1]);
        if which_ves<15
        ylim([-50 450])
        else
            ylim([-60 220])
        end
        xlim([t(1) t(end)])
        print(strcat(fname_save,'_Q_UQ'),'-dpng');

        % CS
%         figure(10.*which_ves+2); clf; hold on;
%         fill([t'; fliplr(t)'],[mu_CS+std_CS; flipud(mu_CS-std_CS)],[0.9 0.6 0.6],'EdgeAlpha',0,'FaceAlpha',0.4);
%         plot(t,mu_CS,'k','LineWidth',3);
%         set(gca,'FontSize',30); grid on;
%         ylabel('Cyclic Stretch (%)');
%         xticks([0 0.5 1]);
%         print(strcat(fname_save,'_CS_UQ'),'-dpng');

        % WSS
        figure(10.*which_ves+3); clf; hold on;
        fill([t'; fliplr(t)'],[mu_WSS+std_WSS; flipud(mu_WSS-std_WSS)],shade_color,'EdgeAlpha',0,'FaceAlpha',0.4);
        plot(t,mu_WSS,'k','LineWidth',3);
        set(gca,'FontSize',30); grid on;
        ylabel('WSS (dyne/cm^2)');
        xticks([0 0.5 1]);
        if which_ves<15
        ylim([-5 25])
        else
%         axis tight
            ylim([-20 72])
        end
        xlim([t(1) t(end)])
        print(strcat(fname_save,'_WSS_UQ'),'-dpng');

        %% WIA
        figure(10.*which_ves+4); clf; hold on;

%         plot(tW,mu_fwdC+std_fwdC,'--r','LineWidth',3);
%         plot(tW,mu_fwdC-std_fwdC,'--r','LineWidth',3);
        fill([tW'; fliplr(tW)'],[mu_fwdC+std_fwdC; flipud(mu_fwdC-std_fwdC)],[1 0 0],'EdgeAlpha',0,'FaceAlpha',0.2);
        plot(tW,mu_fwdC,'r','LineWidth',3);

%         plot(tW,mu_fwdE+std_fwdE,'--c','LineWidth',3);
%         plot(tW,mu_fwdE-std_fwdE,'--c','LineWidth',3);
        fill([tW'; fliplr(tW)'],[mu_fwdE+std_fwdE; flipud(mu_fwdE-std_fwdE)],[0 1 1],'EdgeAlpha',0,'FaceAlpha',0.3);
        plot(tW,mu_fwdE,'c','LineWidth',3);

%         plot(tW,mu_bwdC+std_bwdC,'--b','LineWidth',3);
%         plot(tW,mu_bwdC-std_bwdC,'--b','LineWidth',3);
        fill([tW'; fliplr(tW)'],[mu_bwdC+std_bwdC; flipud(mu_bwdC-std_bwdC)],[0 0 1],'EdgeAlpha',0,'FaceAlpha',0.1);
        plot(tW,mu_bwdC,'b','LineWidth',3);

%         plot(tW,mu_bwdE+std_bwdE,'--m','LineWidth',3);
%         plot(tW,mu_bwdE-std_bwdE,'--m','LineWidth',3);
        fill([tW'; fliplr(tW)'],[mu_bwdE+std_bwdE; flipud(mu_bwdE-std_bwdE)],[1 0 1],'EdgeAlpha',0,'FaceAlpha',0.2);
        plot(tW,mu_bwdE,'m','LineWidth',3);

        set(gca,'FontSize',30); grid on;
        ylabel('WI ');
        if which_ves<15
            ylim([-1e5 2e5]);
        else
            axis tight
%             ylim([-1.7e5 1e5]);
        end
        xlim([tW(1) tW(end)])
        print(strcat(fname_save,'_WIA_UQ'),'-dpng');



        close all;
    end


