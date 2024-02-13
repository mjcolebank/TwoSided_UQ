fname = {'UQ_micro/Art_Alpha_',...
    'UQ_micro/Art_Beta_',...
    'UQ_micro/Ven_Alpha_',...
    'UQ_micro/Ven_Beta_'};
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
% x_interp = linspace(0,1,num_pts_outs);
x_interp_A = linspace(0,1,num_pts_outs);
x_interp_V = x_interp_A+1;
mu_P  = zeros(num_pts_outs,1); std_P = zeros(num_pts_outs,1);
mu_Q  = zeros(num_pts_outs,1); std_Q = zeros(num_pts_outs,1);
mu_CS  = zeros(num_pts_outs,1); std_CS = zeros(num_pts_outs,1);
mu_WSS  = zeros(num_pts_outs,1); std_WSS = zeros(num_pts_outs,1);

Si_P_comp = zeros(4,8,8); ST_P_comp = zeros(4,8,8);
Si_Q_comp = zeros(4,8,8); ST_Q_comp = zeros(4,8,8);
Si_WSS_comp = zeros(4,8,8); ST_WSS_comp = zeros(4,8,8);
Si_CS_comp = zeros(4,8,8); ST_CS_comp = zeros(4,8,8);
for which_data = 1:8
    for which_ves = 1:4
        fname_save = strcat(fname{which_ves},num2str(which_data));
        fname_plot = strcat(fname2{which_ves},num2str(which_data));
        load(fname_save);
        if which_ves<3
            x_interp = x_interp_A;
            shade_color = [0.6 0.6 0.9];
        else
            x_interp = x_interp_V;
            shade_color = [0.9 0.6 0.6];
        end
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

         

            intV_WSS = trapz(std_WSS(1:i).^2);
            GWSSi(:,i)    = trapz(std_WSS(1:i).^2 .* Sobol_WSS.Results.FirstOrder(:,1:i)'./intV_WSS);
            GWSST(:,i)    = trapz(std_WSS(1:i).^2 .* Sobol_WSS.Results.Total(:,1:i)'./intV_WSS);
        end
        GCSi    = Sobol_CS.Results.FirstOrder;
        GCST    = Sobol_CS.Results.Total;

        %% For Si and St next to each other

        Si_P_comp(which_ves,which_data,:)   = GPi(:,end); 
        ST_P_comp(which_ves,which_data,:)   = GPT(:,end); 
        Si_Q_comp(which_ves,which_data,:)   = GQi(:,end); 
        ST_Q_comp(which_ves,which_data,:)   = GQT(:,end); 
        Si_WSS_comp(which_ves,which_data,:) = GWSSi(:,end); 
        ST_WSS_comp(which_ves,which_data,:) = GWSST(:,end); 
        Si_CS_comp(which_ves,which_data,:)  = GCSi(:,end); 
        ST_CS_comp(which_ves,which_data,:)  = GCST(:,end); 
%% For two horizontal plots
        Si_vals = [GPi(:,end) GQi(:,end) GWSSi(:,end) GCSi(:,end)];
        ST_vals = [GPT(:,end) GQT(:,end) GWSST(:,end) GCST(:,end)];
 
%         h = figure(100+which_data);
%         subplot(1,4,which_ves);
%         b = bar([Si_vals(:,1) ST_vals(:,1)]); ylim([0 1]);
%             b(1).FaceColor = [0.7 0.7 0.7];
%             b(2).FaceColor = [0.4 0.4 0.4];            
%         grid on; set(gca,'FontSize',20);
%         h.Position=[-1280 1310 1792 420];
%         print(strcat(fname_save,'_SobolP_bar'),'-dpng');
        
%         h = figure(200+which_data);
%         subplot(1,4,which_ves);
%         b = bar([Si_vals(:,2) ST_vals(:,2)]); ylim([0 1]);
%             b(1).FaceColor = [0.7 0.7 0.7];
%             b(2).FaceColor = [0.4 0.4 0.4];            
%         grid on; set(gca,'FontSize',20);
%         h.Position=[-1280 1310 1792 420];
%         print(strcat(fname_save,'_SobolQ_bar'),'-dpng');

%          h = figure(300+which_data);
%          subplot(1,4,which_ves);
%         b = bar([Si_vals(:,3) ST_vals(:,3)]); ylim([0 1]);
%             b(1).FaceColor = [0.7 0.7 0.7];
%             b(2).FaceColor = [0.4 0.4 0.4];            
%         grid on; set(gca,'FontSize',20);
%         h.Position=[-1280 1310 1792 420];
%         print(strcat(fname_save,'_SobolWSS_bar'),'-dpng');
        
%          h = figure(400+which_data);
%          subplot(1,4,which_ves);
%         b = bar([Si_vals(:,4) ST_vals(:,4)]); ylim([0 1]);
%             b(1).FaceColor = [0.7 0.7 0.7];
%             b(2).FaceColor = [0.4 0.4 0.4];            
%         grid on; set(gca,'FontSize',20);
%         h.Position=[-1280 1310 1792 420];
%         print(strcat(fname_save,'_SobolCS_bar'),'-dpng');

    end

end
%%
si_x = 0.7:2:16;
st_x = 1.3:2:16;
for j=1:4
    for jj=1:4
        if jj==1
            Si = squeeze(Si_P_comp(j,:,:));
            St = squeeze(ST_P_comp(j,:,:));
        elseif jj==2
            Si = squeeze(Si_Q_comp(j,:,:));
            St = squeeze(ST_Q_comp(j,:,:));
        elseif jj==3
            Si = squeeze(Si_WSS_comp(j,:,:));
            St = squeeze(ST_WSS_comp(j,:,:));
        else
            Si = squeeze(Si_CS_comp(j,:,:));
            St = squeeze(ST_CS_comp(j,:,:));
        end
        mu_si = median(Si,1);
        sd_si = std(Si,[],1);
        mu_st = median(St,1);
        sd_st = std(St,[],1);

        min_si = mu_si - min(Si);
        min_st = mu_st - min(St);

        max_si = -mu_si + max(Si);
        max_st = -mu_st + max(St);


        figure(1000.*jj);
        subplot(1,4,j); hold on;
        bar(si_x,mu_si,'BarWidth',0.3,'FaceColor',[0.7 0.7 0.7])
        bar(st_x,mu_st,'BarWidth',0.3,'FaceColor',[0.4 0.4 0.4])
        errorbar(si_x,mu_si,min_si,max_si,'Color',[0 0 0],'LineStyle','none');
        errorbar(st_x,mu_st,min_st,max_st,'Color',[0 0 0],'LineStyle','none');
        axis tight;
        ylim([0 1]);
        grid on; set(gca,'FontSize',20);
        xticklabels('')
    end
end

% bar(x,data)                
% hold on
% er = errorbar(x,data,errlow,errhigh);    
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% hold off

%%

% % Boxplots for Si and ST values
% for j=1:4
%     figure(1000); 
%     subplot(2,2,j); hold on;
%     boxplot(squeeze(Si_P_comp(j,:,:)));
% 
%     figure(1001); 
%     subplot(2,2,j); hold on;
%     boxplot(squeeze(ST_P_comp(j,:,:)));
% 
%     figure(2000); 
%     subplot(2,2,j); hold on;
%     boxplot(squeeze(Si_Q_comp(j,:,:)));
% 
%     figure(2001); 
%     subplot(2,2,j); hold on;
%     boxplot(squeeze(ST_Q_comp(j,:,:)));
% 
%     figure(3000); 
%     subplot(2,2,j); hold on;
%     boxplot(squeeze(Si_WSS_comp(j,:,:)));
% 
%     figure(3001); 
%     subplot(2,2,j); hold on;
%     boxplot(squeeze(ST_WSS_comp(j,:,:)));
% 
%     figure(4000); 
%     subplot(2,2,j); hold on;
%     boxplot(squeeze(Si_CS_comp(j,:,:)));
% 
%     figure(4001); 
%     subplot(2,2,j); hold on;
%     boxplot(squeeze(ST_CS_comp(j,:,:)));
% 
% end
