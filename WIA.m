function [fwd_comp, fwd_exp, bwd_comp, bwd_exp] = WIA(t,p,q,a,r,stiff,plot_ID)
p_og = p;
q_og = q;
p = [p;p];
q = [q;q];
a = [a;a];
for i=1:size(p,2)
    p(:,i) = smooth(p(:,i),20);
    q(:,i) = smooth(q(:,i),20);
    a(:,i) = smooth(a(:,i),20);
end
p = p((floor(end/2)+1):end,:);
q = q((floor(end/2)+1):end,:);
a = a((floor(end/2)+1):end,:);

% Define constants
rho = 1.055;
n = length(t);
p = p.*1333.22;
u = q./a;
dt = t(2)-t(1);

%% SUBSAMPLE
% p = p(1:8:end,:);
% u = u(1:8:end,:);
% a = a(1:8:end,:);
% t = t(1:8:end,:);

dp = p(2:end,:)-p(1:end-1,:);%diff(p);
du = diff(u);
tnew = t;%(1:end-1);


% NOTE: W/m^2 is 0.001 g/s^3
Wconv = 0.001;

% Fudge factor: C was not nondimensionalized correctly
c = wave_speed(a,r,stiff);

rho_cW = rho.*c(1:end-1,:);

dp_fwd = 0.5.*(dp+rho_cW.*du);
dp_bwd = 0.5.*(dp-rho_cW.*du);

P_fwd = fwd_bwd(dp_fwd,p);
P_bwd = fwd_bwd(dp_bwd,p);

du_fwd = 0.5.*(du+(dp./rho_cW));
du_bwd = 0.5.*(du-(dp./rho_cW));

U_fwd = fwd_bwd(du_fwd,u);
U_bwd = fwd_bwd(du_bwd,u);

WI_fwd = Wconv.*(dp_fwd.*du_fwd)./dt^2;
WI_bwd = Wconv.*(dp_bwd.*du_bwd)./dt^2;


%% Try plotting differently
ids = ones(size(dp_fwd)).*[1:(length(t)-1)]';
bool_fwd = dp_fwd>0;
bool_bwd = dp_bwd>0;

fwd_comp = WI_fwd.*bool_fwd;
fwd_exp  = WI_fwd.*~bool_fwd;
bwd_comp = WI_bwd.*bool_bwd;
bwd_exp  = WI_bwd.*~bool_bwd;
return
%% red, magenta, cyan, and blue
figure(plot_ID+1);hold on;
h = area(tnew,fwd_comp); 
h.FaceColor = 'r'; %h.LineWidth = 2;
h = area(tnew,fwd_exp); 
h.FaceColor = 'c'; %h.LineWidth = 2;
h = area(tnew,bwd_comp); 
h.FaceColor = 'b'; %h.LineWidth = 3; h.LineStyle = ':';
h = area(tnew,bwd_exp); 
h.FaceColor = 'm'; %h.LineWidth = 3; h.LineStyle = ':';
plot(tnew,0.*tnew,'k','LineWidth',3);


%% black and white
% figure(plot_ID+1);hold on;
% h = area(tnew(fwd_comp),WI_fwd(fwd_comp)); 
% h.FaceColor = [0.1 0.1 0.1]; %h.LineWidth = 2;
% h = area(tnew(fwd_exp),WI_fwd(fwd_exp)); 
% h.FaceColor = [0.4 0.4 0.4]; %h.LineWidth = 2;
% h = area(tnew(bwd_comp),WI_bwd(bwd_comp)); 
% h.FaceColor = [0.6 0.6 0.6]; h.LineWidth = 3; h.LineStyle = ':';
% h = area(tnew(bwd_exp),WI_bwd(bwd_exp)); 
% h.FaceColor = [0.9 0.9 0.9]; h.LineWidth = 3; h.LineStyle = ':';
% plot(tnew,0.*tnew,'k','LineWidth',3);

% figure(plot_ID+1);hold on;
% plot(fwd_comp,WI_fwd(fwd_comp),'or','LineWidth',2,'MarkerSize',10)
% plot(fwd_exp,WI_fwd(fwd_exp),'oc','LineWidth',2,'MarkerSize',10)
% plot(bwd_comp,WI_bwd(bwd_comp),'ob','LineWidth',2,'MarkerSize',10)
% plot(bwd_exp,WI_bwd(bwd_exp),'om','LineWidth',2,'MarkerSize',10)
ylim([-6.5e4 11.5e4]);
xlim([0 tnew(end)]);
set(gca,'FontSize',24);

disp([sum(P_bwd(fwd_comp))./sum(P_fwd(fwd_comp)) max(P_bwd)./max(P_fwd)]);
end

function c = wave_speed(a,r,fr)
rho = 1.055;
a0 = pi.*r.^2;
% c =  sqrt(0.5.*fr.*sqrt(a./a0)./rho);
c =  sqrt(2.*fr.*sqrt(a./a0)./rho./3);

% c =  0.5*fr(r)*sqrt(a/a0); % use if fr is a function
end

function out = fwd_bwd(dX,X)
[m,n] = size(dX);
out = zeros(m,n);
out(1,:) = dX(1,:);
for i=2:m
    out(i,:) = out(i-1,:)+dX(i,:);
end
out = out + X(1,:);
end

