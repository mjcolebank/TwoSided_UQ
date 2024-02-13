% Run the 1D ST code in the 21 vessel mouse geometry
! make -f Makefile
clear; clc; close all;
% Large artery stiffness 
fa3 = 8e5; % g/cm^2/s

% Small vessel stiffness 
fs3 = 2.5e5; % g/cm^2/s

% Large venous stiffness
fv3 = 8.5e5; % g/cm^2/s

% Structured tree parameters
alpha = 0.88; % Dimensionless
beta =  0.66; % Dimensionless
lrrA = 25;    % Dimensionless, arteriole length-to-radius ratio
lrrV = 25;    % Dimensionless, arteriole length-to-radius ratio
rm   = 0.001; % cm, minimum readius


pars = [fa3,fs3,fv3,alpha,beta,lrrA,lrrV,rm];
% pars_str = mat2str(pars);

%% Now look at PCE
num_par = 8+7; % all
par_id  = [1:8];
poly_order = 3;
size_DP = nchoosek((num_par+poly_order),poly_order);
N_s = 2.*size_DP;
N_t = 256;

% Get amplitude and phase for inlet/outlet boundary condition augmentation
Qdata = dlmread('Qin_orig.dat');
Qfreq = fft(Qdata);
Pdata = dlmread('LA_orig.dat');
Pfreq = fft(Pdata);

% Select some of the frequencies to perturb
Q_amp_ind = [1 2 3];
P_amp_ind = [1 2 3];
P_pha_ind = [3];

Q_amps = abs(real(Qfreq(Q_amp_ind)));
P_amps = abs(real(Pfreq(P_amp_ind)));
P_phas = abs(imag(Pfreq(P_pha_ind)));



% Input the parameter means and standard deviations.
par_nom = [pars(:); Q_amps(:); P_amps(:); P_phas(:)];
% Bounds for the model parameters and boundary conditions 3/9/23
upp = 1.30.*par_nom; % BIG
low = 0.7.*par_nom;

% Now, specify certain upper and lower bounds
% pars = [fa3,fs3,fv3,alpha,beta,lrrA,lrrV,rm];
% FROM OLD CODE
% upp = [0 0 1e6 0 0 3e5 0 0 1e6 0.89 0.68 50 50 0.005]; % BIG
% low = [0 0 6e5 0 0 9e4 0 0 6e5 0.84 0.65 20 20 0.0005]; % BIG
% upp(1:3) = 1e6; low(1:3) = 8e4;
upp(4) = 0.92; low(4) = 0.80;
upp(5) = 0.76; low(5) = 0.60;
upp(6:7) = 50; low(6:7) = 10;
upp(8) = 0.01; low(8) = 0.001;

mu = (upp+low)./2;
sig = sqrt((upp-low).^2./12);
%%

% X = sobolset(num_par, 'Skip',2e12,'Leap',0.1e13); % draw 4503 values

rng('shuffle');
X = unifrnd(0,1,N_s+100,num_par);
chi = (X(:,:)-0.5).*2;
par_sample = sqrt(3).*chi.*sig'+mu';

%%

% ind_mat = get_PSImap(num_par,poly_order);
% ind_rows = size(ind_mat,1);
% 
% reg_mat = zeros(N_s,ind_rows);


p_storage = cell(N_s,1);
q_storage = cell(N_s,1);
A_storage = cell(N_s,1);

format shorte
% cluster_P = gcp('nocreate');
% if isempty(cluster_P)
%         parpool(4);
% end
%delete(gcp('nocreate'));
% tic
parpool(2)

parfor i=1:2%N_s
    disp(i)
    Q = par_sample(i,:);

    % Get new inlet and outlet boundary condition
    get_inflow(Q(9:11),Q_amp_ind,i);
    get_outletpressure(Q(12:14),Q(15),P_amp_ind,P_pha_ind,i);
   
    par_in = pars;
    par_in(par_id) = Q(par_id);
    par_in(end+1) = i;
    par_in = round(par_in,5);
    pars_str = mat2str(par_in);
    %     clc;
    %     str_pars = mat2str(pars);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Run an initial simulation to get the BC parameters
    out = unix(sprintf('sor06.exe  %s',pars_str(2:end-1)));
%     pause(1);
	disp('Out complete');
    if out ==2
        fname = strcat('output_',num2str(i),'.2d');
        data = load(fname);
	disp('data loaded');
        [~,~,p,q,a,~] = gnuplot(data);
       %% figure; plot(p)
       disp('gnuplot finished');
        p_storage{i} = p;
        q_storage{i} = q;
        A_storage{i} = a;

	disp('stored solutions');
       
       	delete(fname);
	disp('solution deleted');

        delete(strcat('Qin_',num2str(i),'.dat'));
	disp('Qin deleted');

        delete(strcat('LAout_',num2str(i),'.dat'));
        disp('LA deleted');
    end

end
% return;
%% save all the files
out_size = round(N_s/10);
outs_1 = cell(out_size,3); outs_2 = cell(out_size,3); outs_3 = cell(out_size,3);
outs_4 = cell(out_size,3); outs_5 = cell(out_size,3); outs_6 = cell(out_size,3);
outs_7 = cell(out_size,3); outs_8 = cell(out_size,3); outs_9 = cell(out_size,3);
outs_10 = cell(N_s-9*out_size,3);

i = 1; counter_start = (i-1)*out_size+1; counter_end   = (i)*out_size;
counter_range = counter_start:counter_end;
for j=1:out_size
    outs_1{j,1} = p_storage{counter_range(j)};
    outs_1{j,2} = q_storage{counter_range(j)};
    outs_1{j,3} = A_storage{counter_range(j)};
end

i = i+1; counter_start = (i-1)*out_size+1; counter_end   = (i)*out_size;
counter_range = counter_start:counter_end;
for j=1:out_size
    outs_2{j,1} = p_storage{counter_range(j)};
    outs_2{j,2} = q_storage{counter_range(j)};
    outs_2{j,3} = A_storage{counter_range(j)};
end

i = i+1; counter_start = (i-1)*out_size+1; counter_end   = (i)*out_size;
counter_range = counter_start:counter_end;
for j=1:out_size
    outs_3{j,1} = p_storage{counter_range(j)};
    outs_3{j,2} = q_storage{counter_range(j)};
    outs_3{j,3} = A_storage{counter_range(j)};
end

i = i+1; counter_start = (i-1)*out_size+1; counter_end   = (i)*out_size;
counter_range = counter_start:counter_end;
for j=1:out_size
    outs_4{j,1} = p_storage{counter_range(j)};
    outs_4{j,2} = q_storage{counter_range(j)};
    outs_4{j,3} = A_storage{counter_range(j)};
end

i = i+1; counter_start = (i-1)*out_size+1; counter_end   = (i)*out_size;
counter_range = counter_start:counter_end;
for j=1:out_size
    outs_5{j,1} = p_storage{counter_range(j)};
    outs_5{j,2} = q_storage{counter_range(j)};
    outs_5{j,3} = A_storage{counter_range(j)};
end

i = i+1; counter_start = (i-1)*out_size+1; counter_end   = (i)*out_size;
counter_range = counter_start:counter_end;
for j=1:out_size
    outs_6{j,1} = p_storage{counter_range(j)};
    outs_6{j,2} = q_storage{counter_range(j)};
    outs_6{j,3} = A_storage{counter_range(j)};
end

i = i+1; counter_start = (i-1)*out_size+1; counter_end   = (i)*out_size;
counter_range = counter_start:counter_end;
for j=1:out_size
    outs_7{j,1} = p_storage{counter_range(j)};
    outs_7{j,2} = q_storage{counter_range(j)};
    outs_7{j,3} = A_storage{counter_range(j)};
end

i = i+1; counter_start = (i-1)*out_size+1; counter_end   = (i)*out_size;
counter_range = counter_start:counter_end;
for j=1:out_size
    outs_8{j,1} = p_storage{counter_range(j)};
    outs_8{j,2} = q_storage{counter_range(j)};
    outs_8{j,3} = A_storage{counter_range(j)};
end

i = i+1; counter_start = (i-1)*out_size+1; counter_end   = (i)*out_size;
counter_range = counter_start:counter_end;
for j=1:out_size
    outs_9{j,1} = p_storage{counter_range(j)};
    outs_9{j,2} = q_storage{counter_range(j)};
    outs_9{j,3} = A_storage{counter_range(j)};
end

i = i+1; counter_start = (i-1)*out_size+1; counter_end   = (i)*out_size;
counter_range = counter_start:N_s;
for j=1:out_size
    outs_10{j,1} = p_storage{counter_range(j)};
    outs_10{j,2} = q_storage{counter_range(j)};
    outs_10{j,3} = A_storage{counter_range(j)};
end

save('PCE_out_1',"outs_1");
save('PCE_out_2',"outs_2");
save('PCE_out_3',"outs_3");
save('PCE_out_4',"outs_4");
save('PCE_out_5',"outs_5");
save('PCE_out_6',"outs_6");
save('PCE_out_7',"outs_7");
save('PCE_out_8',"outs_8");
save('PCE_out_9',"outs_9");
save('PCE_out_10',"outs_10");

save('PCE_testBCs_3_10_15par',"upp","low",'par_sample','chi');
