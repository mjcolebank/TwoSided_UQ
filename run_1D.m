% This is a driver file that will run the fluids model from c++ and
% fortran by passing parameter values needed by the model. This code is a
% working file and comes with NO GUARANTEES.
%
% Authors: MJ Colebank, M Bartolo, MU Qureshi, MS Olufsen
%
% Last edited: 1/27/2020, MJC
%%
clear; clc; %close all;
%% Ensure the make file is compiled
! make clean -f Makefile
! make -f Makefile
! chmod +x sor06
%%
% Up: 1.0400e+06   3.2500e+05   1.1050e+06   9.2000e-01   7.6000e-01   5.0000e+01   5.0000e+01   1.0000e-02
% Low:5.6000e+05   1.7500e+05   5.9500e+05   8.0000e-01   6.0000e-01   1.0000e+01   1.0000e+01   1.0000e-03
% Large artery stiffness - exponential relationship
fa3 = 7.0e5;%8e5; % g/cm^2/s

% Small vessel stiffness (same on arteriole and venule side) - exponential relationship
fs3 = 2.0e5;%2.5e5;%8e4; % g/cm^2/s

% Large artery stiffness - exponential relationship
fv3 = 6e5;%8.5e5; % g/cm^2/s

% Structured tree parameters
alpha = 0.83;%0.885; % Dimensionless
beta =  0.67;%0.655; % Dimensionless
lrrA = 26;%35;    % Dimensionless, arteriole length-to-radius ratio
lrrV = 15;%25;    % Dimensionless, arteriole length-to-radius ratio
rm   = 0.0085;%0.002; % cm, minimum readius

%% Put the parameters into a vector and run the model
par_nom = [fa3,fs3,fv3,alpha,beta,lrrA,lrrV,rm,1];
param_str = mat2str(par_nom);
Qin = dlmread('Qin_orig.dat');
dlmwrite('Qin_1.dat',Qin);
LAP = dlmread('LA_orig.dat');
dlmwrite('LAout_1.dat',LAP);
% Run the model
% NOTE: Windows users need 'sor06.exe', Mac/Linux users need ./sor06
tic
out = unix(sprintf('sor06.exe %s',param_str(2:end-1)));
toc
%% Load all the model results
% NOTE: Results are stored such that the first 1:N entries are the PROXIMAL
% large artery (1-15) and large vein (16-27) predictions, N+1:2N are the
% MIDPOINT predictions in the same vessels, and 2N+1:3N are the DISTAL
% predictions.
% Vessel 1 - Main Pulmonary artery, Vessel 2/3 - Left/Right pulmonary artery
% Vessel 16:19 - Left inferior, left superior, right inferior, and right
% superior pulmonary veins (connected to left atrium)
% Vessels 8:15 - terminal arteries
% Vessels 20:27 - terminal veins
name = 'output_1.2d';
data = load(name);
%%
% p - pressure (mmHg)
% q - flow (cm^3/s)
% a - area (cm^2)
[~,~,p,q,a,~] = gnuplot(data);

% Time vector
t = linspace(0,0.11,512);
Num_ves = 27;   % Number of large vessels in the tree
art     = 1:15; % Arteries
ven     = 16:27;% Veins

% Plot midpoint  in the arteries and veins
figure;%(1);clf;
subplot(1,2,1);
plot(t,p(:,art+Num_ves),'LineWidth',3);
ylabel('Pressure (mmHg)');
xlabel('Time (s)');
grid on; set(gca,'FontSize',20);

subplot(1,2,2);
plot(t,p(:,ven+Num_ves),'LineWidth',3);
xlabel('Time (s)');
grid on; set(gca,'FontSize',20);

figure;%(2);clf;
subplot(1,2,1);
plot(t,q(:,art+Num_ves),'LineWidth',3);
ylabel('mL/s')
xlabel('Time (s)');
grid on; set(gca,'FontSize',20);

subplot(1,2,2);
plot(t,-q(:,ven+Num_ves),'LineWidth',3);
xlabel('Time (s)');
grid on; set(gca,'FontSize',20);