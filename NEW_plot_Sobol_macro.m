fname = 'UQ_macro_new/ves_pce4_';
num_par = 8;
Names = {'f3_{A}','f3_{MV}','f3_{V}',...
    '\alpha','\beta','lrr_A','lrr_V','rm'};%,...
% nsamp = 331;%2300;
nsamp = 1900;%990;%
% fname = 'UQ_macro_BC/ves_';
% num_par = 18;
% Names = {'f3_{A}','f3_{MV}','f3_{V}',...
%     '\alpha','\beta','lrr_A','lrr_V','rm',...
%          'QA0','QA1','QA2',...
%          'LA0','LA1','LA2','LA5',...
%          '\phi_1','\phi_2', '\phi_5'};
% num_samp = 2381;
cmap = get_color_map(num_par);

num_pts_outs = 512/8;

t = linspace(0,0.85,num_pts_outs);
% Full output
% tW = linspace(0,0.85,num_pts_outs-6);
% For BC
% tW = linspace(0,0.85,num_pts_outs-2);
tW = linspace(0,0.85,num_pts_outs-1);

num_WIA = length(tW);

grey_color_A = [linspace(0,0.9,15); linspace(0,0.9,15); linspace(0,0.9,15)];
grey_color_V = [linspace(0,0.9,12); linspace(0,0.9,12); linspace(0,0.9,12)];


mu_P  = zeros(num_pts_outs,1); std_P = zeros(num_pts_outs,1);
mu_Q  = zeros(num_pts_outs,1); std_Q = zeros(num_pts_outs,1);
mu_CS  = zeros(num_pts_outs,1); std_CS = zeros(num_pts_outs,1);
mu_WSS  = zeros(num_pts_outs,1); std_WSS = zeros(num_pts_outs,1);

mu_fwdC  = zeros(num_WIA,1); std_fwdC = zeros(num_WIA,1);
mu_fwdE  = zeros(num_WIA,1); std_fwdE = zeros(num_WIA,1);
mu_bwdC  = zeros(num_WIA,1); std_bwdC = zeros(num_WIA,1);
mu_bwdE  = zeros(num_WIA,1); std_bwdE = zeros(num_WIA,1);


Si_P_comp = zeros(27,num_par); ST_P_comp = zeros(27,num_par);
Si_Q_comp = zeros(27,num_par); ST_Q_comp = zeros(27,num_par);
Si_WSS_comp = zeros(27,num_par); ST_WSS_comp = zeros(27,num_par);
Si_CS_comp = zeros(27,num_par); ST_CS_comp = zeros(27,num_par);

Si_fwdC_comp = zeros(27,num_par); ST_fwdC_comp = zeros(27,num_par);
Si_bwdC_comp = zeros(27,num_par); ST_bwdC_comp = zeros(27,num_par);
Si_fwdE_comp = zeros(27,num_par); ST_fwdE_comp = zeros(27,num_par);
Si_bwdE_comp = zeros(27,num_par); ST_bwdE_comp = zeros(27,num_par);

TAWSS = zeros(27,1);


%%
    for which_ves = [8:15 20:27]%1:27%[1 2 3 16 17 18 19]%[1 2 3]
        
        fname_save = strcat(fname,num2str(which_ves),'_subsamp');
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

        figure(which_ves);
        subplot(2,2,1); plot(std_P);
        subplot(2,2,2); plot(std_Q);
        subplot(2,2,3); plot(sqrt(PCE_CS.PCE.Moments.Var),'o');
        subplot(2,2,4); plot(std_WSS);
        
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

%         figure(5000+which_ves); hold on;
%         plot(mu_WSS,'k','LineWidth',2);
%         plot(mu_WSS+std_WSS,'--r','LineWidth',1.5);
%         plot(mu_WSS-std_WSS,'--r','LineWidth',1.5);

%         TAWSS(which_ves) = mean(mu_WSS);
        


       
        %% Lastly, Calculate the generalized Sobol indices based on Alexanderian et al.
%         intD = trapz(ti, Di);
%         w = Di ./ intD;
%         gS(n) = trapz(ti, w .* Si);
        GPi   = zeros(num_par,num_pts_outs);  GPT = zeros(num_par,num_pts_outs);
        GQi   = zeros(num_par,num_pts_outs);  GQT = zeros(num_par,num_pts_outs);
        GWSSi = zeros(num_par,num_pts_outs);  GWSST = zeros(num_par,num_pts_outs);

        GfwdCi = zeros(num_par,num_WIA);  GfwdCT = zeros(num_par,num_WIA);
        GbwdCi = zeros(num_par,num_WIA);  GbwdCT = zeros(num_par,num_WIA);

        GfwdEi = zeros(num_par,num_WIA);  GfwdET = zeros(num_par,num_WIA);
        GbwdEi = zeros(num_par,num_WIA);  GbwdET = zeros(num_par,num_WIA);

        for i=2:num_pts_outs
            intV_P = trapz(std_P(1:i).^2);
            GPi(:,i)    = trapz(std_P(1:i).^2 .* Sobol_P.Results.FirstOrder(:,1:i)'./intV_P);
            GPT(:,i)    = trapz(std_P(1:i).^2 .* Sobol_P.Results.Total(:,1:i)'./intV_P);

            intV_Q = trapz(std_Q(1:i).^2);
            GQi(:,i)    = trapz(std_Q(1:i).^2 .* Sobol_Q.Results.FirstOrder(:,1:i)'./intV_Q);
            GQT(:,i)    = trapz(std_Q(1:i).^2 .* Sobol_Q.Results.Total(:,1:i)'./intV_Q);

%             intV_CS = trapz(std_CS(1:i).^2);
%             GCSi(:,i)    = trapz(std_CS(1:i).^2 .* Sobol_CS.Results.FirstOrder(:,1:i)'./intV_CS);
%             GCST(:,i)    = trapz(std_CS(1:i).^2 .* Sobol_CS.Results.Total(:,1:i)'./intV_CS);

            intV_WSS = trapz(std_WSS(1:i).^2);
            GWSSi(:,i)    = trapz(std_WSS(1:i).^2 .* Sobol_WSS.Results.FirstOrder(:,1:i)'./intV_WSS);
            GWSST(:,i)    = trapz(std_WSS(1:i).^2 .* Sobol_WSS.Results.Total(:,1:i)'./intV_WSS);
        end
        for i=2:num_WIA
            % Wave intensity - Compression
            intV_fwdC = trapz(std_fwdC(1:i).^2);
            GfwdCi(:,i)    = trapz(std_fwdC(1:i).^2 .* Sobol_fwdC.Results.FirstOrder(:,1:i)'./intV_fwdC);
            GfwdCT(:,i)    = trapz(std_fwdC(1:i).^2 .* Sobol_fwdC.Results.Total(:,1:i)'./intV_fwdC);

            intV_bwdC = trapz(std_bwdC(1:i).^2);
            GbwdCi(:,i)    = trapz(std_bwdC(1:i).^2 .* Sobol_bwdC.Results.FirstOrder(:,1:i)'./intV_bwdC);
            GbwdCT(:,i)    = trapz(std_bwdC(1:i).^2 .* Sobol_bwdC.Results.Total(:,1:i)'./intV_bwdC);
            
            % Wave intensity - Expansion
            intV_fwdE = trapz(std_fwdE(1:i).^2);
            GfwdEi(:,i)    = trapz(std_fwdE(1:i).^2 .* Sobol_fwdE.Results.FirstOrder(:,1:i)'./intV_fwdE);
            GfwdET(:,i)    = trapz(std_fwdE(1:i).^2 .* Sobol_fwdE.Results.Total(:,1:i)'./intV_fwdE);

            intV_bwdE = trapz(std_bwdE(1:i).^2);
            GbwdEi(:,i)    = trapz(std_bwdE(1:i).^2 .* Sobol_bwdE.Results.FirstOrder(:,1:i)'./intV_bwdE);
            GbwdET(:,i)    = trapz(std_bwdE(1:i).^2 .* Sobol_bwdE.Results.Total(:,1:i)'./intV_bwdE);
        end

        GCSi    = Sobol_CS.Results.FirstOrder;
        GCST    = Sobol_CS.Results.Total;



        Si_P_comp(which_ves,:)   = GPi(:,end); 
        ST_P_comp(which_ves,:)   = GPT(:,end); 
        Si_Q_comp(which_ves,:)   = GQi(:,end); 
        ST_Q_comp(which_ves,:)   = GQT(:,end); 
        Si_WSS_comp(which_ves,:) = GWSSi(:,end); 
        ST_WSS_comp(which_ves,:) = GWSST(:,end); 
        Si_CS_comp(which_ves,:)  = GCSi(:,end); 
        ST_CS_comp(which_ves,:)  = GCST(:,end); 

        %% WIA sobols
        Si_fwdC_comp(which_ves,:) = GfwdCi(:,end);
        ST_fwdC_comp(which_ves,:) = GfwdCT(:,end);
        Si_bwdC_comp(which_ves,:) = GbwdCi(:,end);
        ST_bwdC_comp(which_ves,:) = GbwdCT(:,end);

        Si_fwdE_comp(which_ves,:) = GfwdEi(:,end);
        ST_fwdE_comp(which_ves,:) = GfwdET(:,end);
        Si_bwdE_comp(which_ves,:) = GbwdEi(:,end);
        ST_bwdE_comp(which_ves,:) = GbwdET(:,end);

        
%% For Si and St next to each other
%% For two horizontal plots
        Si_vals = [GPi(:,end) GQi(:,end) GWSSi(:,end) GCSi(:,end)];
        ST_vals = [GPT(:,end) GQT(:,end) GWSST(:,end) GCST(:,end)];
 
%         h = figure(100+which_ves);
%         b = bar([Si_vals(:,1) ST_vals(:,1)]); ylim([0 1]);
%             b(1).FaceColor = [0.7 0.7 0.7];
%             b(2).FaceColor = [0.4 0.4 0.4];            
%         grid on; set(gca,'FontSize',20);
%         h.Position=[-1280 1310 1792 420];
%         print(strcat(fname_save,'_SobolP_bar'),'-dpng');
        
%         h = figure(200+which_ves);
%         b = bar([Si_vals(:,2) ST_vals(:,2)]); ylim([0 1]);
%             b(1).FaceColor = [0.7 0.7 0.7];
%             b(2).FaceColor = [0.4 0.4 0.4];            
%         grid on; set(gca,'FontSize',20);
%         h.Position=[-1280 1310 1792 420];
%         print(strcat(fname_save,'_SobolQ_bar'),'-dpng');

%          h = figure(300+which_ves);
%         b = bar([Si_vals(:,3) ST_vals(:,3)]); ylim([0 1]);
%             b(1).FaceColor = [0.7 0.7 0.7];
%             b(2).FaceColor = [0.4 0.4 0.4];            
%         grid on; set(gca,'FontSize',20);
%         h.Position=[-1280 1310 1792 420];
%         print(strcat(fname_save,'_SobolWSS_bar'),'-dpng');
        
%          h = figure(400+which_ves);
%         b = bar([Si_vals(:,4) ST_vals(:,4)]); ylim([0 1]);
%             b(1).FaceColor = [0.7 0.7 0.7];
%             b(2).FaceColor = [0.4 0.4 0.4];            
%         grid on; set(gca,'FontSize',20);
%         h.Position=[-1280 1310 1792 420];
%         print(strcat(fname_save,'_SobolCS_bar'),'-dpng');


% WIA plots
% figure(1000+which_ves);clf;hold on;
% plot(mu_fwdC);
% plot(mu_fwdE);
% plot(mu_bwdC);
% plot(mu_bwdE);
% legend('FCW','FEW','BCW','BEW')
% 
% if which_ves<16
%     figure(2000);
%     subplot(2,2,1); hold on; plot(mu_fwdC,'LineWidth',2,'Color',grey_color_A(:,which_ves));
%     subplot(2,2,2); hold on; plot(mu_fwdE,'LineWidth',2,'Color',grey_color_A(:,which_ves));
%     subplot(2,2,3); hold on; plot(mu_bwdC,'LineWidth',2,'Color',grey_color_A(:,which_ves));
%     subplot(2,2,4); hold on; plot(mu_bwdE,'LineWidth',2,'Color',grey_color_A(:,which_ves));
% 
%     figure(2001);
%     subplot(2,2,1); hold on; plot(std_fwdC,'LineWidth',2,'Color',grey_color_A(:,which_ves));
%     subplot(2,2,2); hold on; plot(std_fwdE,'LineWidth',2,'Color',grey_color_A(:,which_ves));
%     subplot(2,2,3); hold on; plot(std_bwdC,'LineWidth',2,'Color',grey_color_A(:,which_ves));
%     subplot(2,2,4); hold on; plot(std_bwdE,'LineWidth',2,'Color',grey_color_A(:,which_ves));
% else
%      figure(2002);
%     subplot(2,2,1); hold on; plot(mu_fwdC,'LineWidth',2,'Color',grey_color_V(:,which_ves-15));
%     subplot(2,2,2); hold on; plot(mu_fwdE,'LineWidth',2,'Color',grey_color_V(:,which_ves-15));
%     subplot(2,2,3); hold on; plot(mu_bwdC,'LineWidth',2,'Color',grey_color_V(:,which_ves-15));
%     subplot(2,2,4); hold on; plot(mu_bwdE,'LineWidth',2,'Color',grey_color_V(:,which_ves-15));
% 
%     figure(2003);
%     subplot(2,2,1); hold on; plot(std_fwdC,'LineWidth',2,'Color',grey_color_V(:,which_ves-15));
%     subplot(2,2,2); hold on; plot(std_fwdE,'LineWidth',2,'Color',grey_color_V(:,which_ves-15));
%     subplot(2,2,3); hold on; plot(std_bwdC,'LineWidth',2,'Color',grey_color_V(:,which_ves-15));
%     subplot(2,2,4); hold on; plot(std_bwdE,'LineWidth',2,'Color',grey_color_V(:,which_ves-15));
% 
% end
% 
% figure(1000+which_ves);clf;
% subplot(2,2,1); hold on; plot(PCE_fwdC.ExpDesign.Y','r'); plot(mu_fwdC,'--k','LineWidth',2);
% grid on; set(gca,'FontSize',20);
% subplot(2,2,2); hold on; plot(PCE_fwdE.ExpDesign.Y','c'); plot(mu_fwdE,'--k','LineWidth',2);
% grid on; set(gca,'FontSize',20);
% subplot(2,2,3); hold on; plot(PCE_bwdC.ExpDesign.Y','b'); plot(mu_bwdC,'--k','LineWidth',2);
% grid on; set(gca,'FontSize',20);
% subplot(2,2,4); hold on; plot(PCE_bwdE.ExpDesign.Y','m'); plot(mu_bwdE,'--k','LineWidth',2);
% grid on; set(gca,'FontSize',20);
% 
% set(gcf,'Position',[1.4500e+02   1.1370e+03   9.6400e+02   6.9920e+02]);
% print(strcat(fname_save,'_WIA_all'),'-dpng');
%         close all;
    end

%%
si_x = 0.7:2:num_par.*2;
st_x = 1.3:2:num_par.*2;
j=1:15;
    for jj=1:4
        if jj==1
            Si = squeeze(Si_P_comp(j,:));
            St = squeeze(ST_P_comp(j,:));
        elseif jj==2
            Si = squeeze(Si_Q_comp(j,:));
            St = squeeze(ST_Q_comp(j,:));
        elseif jj==3
            Si = squeeze(Si_WSS_comp(j,:));
            St = squeeze(ST_WSS_comp(j,:));
        else
            Si = squeeze(Si_CS_comp(j,:));
            St = squeeze(ST_CS_comp(j,:));
        end
        mu_si = median(Si,1);
        sd_si = std(Si,[],1);
        mu_st = median(St,1);
        sd_st = std(St,[],1);

        min_si = mu_si - min(Si);
        min_st = mu_st - min(St);

        max_si = -mu_si + max(Si);
        max_st = -mu_st + max(St);


        figure(1000);
        subplot(4,1,jj); hold on;
        bar(si_x,mu_si,'BarWidth',0.3,'FaceColor',[0.7 0.7 0.7])
        bar(st_x,mu_st,'BarWidth',0.3,'FaceColor',[0.4 0.4 0.4])
        errorbar(si_x,mu_si,min_si,max_si,'Color',[0 0 0],'LineStyle','none');
        errorbar(st_x,mu_st,min_st,max_st,'Color',[0 0 0],'LineStyle','none');
        axis tight;
        ylim([0 1]);
        grid on; set(gca,'FontSize',20);
        xticklabels('')
    end

    j=16:27;
    for jj=1:4
        if jj==1
            Si = squeeze(Si_P_comp(j,:));
            St = squeeze(ST_P_comp(j,:));
        elseif jj==2
            Si = squeeze(Si_Q_comp(j,:));
            St = squeeze(ST_Q_comp(j,:));
        elseif jj==3
            Si = squeeze(Si_WSS_comp(j,:));
            St = squeeze(ST_WSS_comp(j,:));
        else
            Si = squeeze(Si_CS_comp(j,:));
            St = squeeze(ST_CS_comp(j,:));
        end
        mu_si = median(Si,1);
        sd_si = std(Si,[],1);
        mu_st = median(St,1);
        sd_st = std(St,[],1);

        min_si = mu_si - min(Si);
        min_st = mu_st - min(St);

        max_si = -mu_si + max(Si);
        max_st = -mu_st + max(St);


        figure(2000);
        subplot(4,1,jj); hold on;
        bar(si_x,mu_si,'BarWidth',0.3,'FaceColor',[0.7 0.7 0.7])
        bar(st_x,mu_st,'BarWidth',0.3,'FaceColor',[0.4 0.4 0.4])
        errorbar(si_x,mu_si,min_si,max_si,'Color',[0 0 0],'LineStyle','none');
        errorbar(st_x,mu_st,min_st,max_st,'Color',[0 0 0],'LineStyle','none');
        axis tight;
        ylim([0 1]);
        grid on; set(gca,'FontSize',20);
        xticklabels('')
    end

    %% WIA
    %%

j=1:15;
    for jj=1:4
        if jj==1
            Si = squeeze(Si_fwdC_comp(j,:));
            St = squeeze(ST_fwdC_comp(j,:));
        elseif jj==2
            Si = squeeze(Si_fwdE_comp(j,:));
            St = squeeze(ST_fwdE_comp(j,:));
        elseif jj==3
            Si = squeeze(Si_bwdC_comp(j,:));
            St = squeeze(ST_bwdC_comp(j,:));
        else
            Si = squeeze(Si_bwdE_comp(j,:));
            St = squeeze(ST_bwdE_comp(j,:));
        end
        mu_si = median(Si,1);
        sd_si = std(Si,[],1);
        mu_st = median(St,1);
        sd_st = std(St,[],1);

        min_si = mu_si - min(Si);
        min_st = mu_st - min(St);

        max_si = -mu_si + max(Si);
        max_st = -mu_st + max(St);


        figure(1100);
        subplot(4,1,jj); hold on;
        bar(si_x,mu_si,'BarWidth',0.3,'FaceColor',[0.7 0.7 0.7])
        bar(st_x,mu_st,'BarWidth',0.3,'FaceColor',[0.4 0.4 0.4])
        errorbar(si_x,mu_si,min_si,max_si,'Color',[0 0 0],'LineStyle','none');
        errorbar(st_x,mu_st,min_st,max_st,'Color',[0 0 0],'LineStyle','none');
        axis tight;
        ylim([0 1]);
        grid on; set(gca,'FontSize',20);
        xticklabels('')
    end

    j=16:27;
    for jj=1:4
        if jj==1
            Si = squeeze(Si_fwdC_comp(j,:));
            St = squeeze(ST_fwdC_comp(j,:));
        elseif jj==2
            Si = squeeze(Si_fwdE_comp(j,:));
            St = squeeze(ST_fwdE_comp(j,:));
        elseif jj==3
            Si = squeeze(Si_bwdC_comp(j,:));
            St = squeeze(ST_bwdC_comp(j,:));
        else
            Si = squeeze(Si_bwdE_comp(j,:));
            St = squeeze(ST_bwdE_comp(j,:));
        end
        mu_si = median(Si,1);
        sd_si = std(Si,[],1);
        mu_st = median(St,1);
        sd_st = std(St,[],1);

        min_si = mu_si - min(Si);
        min_st = mu_st - min(St);

        max_si = -mu_si + max(Si);
        max_st = -mu_st + max(St);


        figure(2100);
        subplot(4,1,jj); hold on;
        bar(si_x,mu_si,'BarWidth',0.3,'FaceColor',[0.7 0.7 0.7])
        bar(st_x,mu_st,'BarWidth',0.3,'FaceColor',[0.4 0.4 0.4])
        errorbar(si_x,mu_si,min_si,max_si,'Color',[0 0 0],'LineStyle','none');
        errorbar(st_x,mu_st,min_st,max_st,'Color',[0 0 0],'LineStyle','none');
        axis tight;
        ylim([0 1]);
        grid on; set(gca,'FontSize',20);
        xticklabels('')
    end
