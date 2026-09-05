%% Neutronics data for the point kinetic model

Lambda = 1.05E-6;
beta_i = [7.91E-5, 7.03E-4, 5.04E-4, 1.17E-3, 4.57E-4, 1.1E-4];
beta = 3.02E-3;
lambda_i = [1.27E-2, 3E-2, 1.1E-1, 3.19E-1, 1.18, 7.02];
phi_tot = 3.61E21;
v_rec = 4.16E-9;
N_0 = 1.5E13;
Pi_0 = [1.58E-1, 6.42E-2, 6.98E-2, 7.1E-2, ...
    7.14E-2, 7.12E-2, 7.07E-2, 7E-2, ...
    6.91E-2, 6.76E-2, 6.49E-2, 1.52E-1];
%% Update 21-04-2024 ( normalizing  power distribution to sum(Pi_0)=1 )
Pi_0 = Pi_0 + (1-sum(Pi_0))*(1/12)*ones(1,12);
%%
P_0 = 1.E8;
tau_c = 5.43;
tau_e = 10;
lambda_c = 0.184;
lambda_e = 0.1;
alfa_f = -2.08E-4;
alfa_c = -8.28E-6;
rho_0 = 1.11E-3; %reactivity, in cents, reported in Liu et al 2022
%rho_0 = 204.2E-6; % Considering rho_0 reported in Wang, 2017
C_i0 = [4.29E14, 4.15E15, 3.83E15, 1.16E16, 5.7E15, 1.53E15]';
Ce_i0 = [7.01E14, 5.88E15, 3.36E15, 5.09E15, 8.20E14, 3.96E13]';

%% Parameters

L = 6;          %delayed neutronic groups
Nodes = 12;     %nodes for energy balance

%% Thermodynamic parameters
fuelpipes = 1027;
coolpipes = [2166, 2166*ones(1,10),0];
A_fw = [13.6, 9.0*ones(1,10), 13.6];
A_wc = [10.9, 10.3*ones(1,10), 10.9];
mdot_f = [327.7, 327.7*ones(1,10), 327.7];
mdot_c = [5550.5, 5550.5*ones(1,10), 5550.5];
M_f = [307.7, 116.4*ones(1,10), 307.7];
M_w = [39.3, 31.1*ones(1,10), 39.3];
M_c = [219.9, 770.5*ones(1,10), 219.9];
c_pf = [400, 400*ones(1,10), 400];
c_pw = [690, 690*ones(1,10), 690];
c_pc = [140.2, 140.2*ones(1,10), 140.2];
Heat = [18744000, 5970000*ones(1,10), 16956000];
h_fw = [6546.0, 4305.1*ones(1,10), 7629.5];
h_wc = [27405.4, 13833.4*ones(1,10), 32177.3];
Tin_f = 1300;
Tin_c = 973;

%% Assumptions
%Tbar_c0 = (973 + 1100)/2; %by inspection
%Tbar_c0 = 1050;    %by inspection
%Tbar_f0 = 1300;    %by inspection
%Tbar_c0 = 1043;     %ROM solution
%Tbar_c0 = 1060;
%Tbar_f0 = 1266;     %ROM solution
%Tf_0 = Tbar_f0*ones(12,1);
%Tc_0 = Tbar_c0*ones(12,1);
%Tw_0 = 1100*ones(12,1);

%% Update 15-12-22, failure injection method
%M_f_tot = sum(M_f);
%mdot_f_delta = 0;

%% Assumptions (vact. Sep 27. 2023)
%0. Total reactivity formula: modified 

%1. Tin_x and Tout_x reported in Liu et al. (2021)

%Defined by design
%Tin_f = 1300;
%Tout_f = 1300;

%Tin_c = 973;
%Tout_c = 1100;

%2. Tbar_x: key for stability at first seconds of simulation

%They do not have to be equal to Tout_x!

%Defined by inspection (in order to P_N(t=tsim) = 1)
%Tbar_f0 = 1288;
%Tbar_c0 = 1060;

%3. Initial conditions Tf_0, Tw_0, Tc_0, Te_f0, Te_c0

% Defined by inspection (check the simulink model initial conditions)
% Te_f0, Te_c0: based on design values. Only by inspection
% These values play a fundamental role for reactivity stability

%% Assumptions (vact. April 2024)
%Tbar_f0 = (Tin_f + Tout_f)/2;
%Tbar_c0 = (Tin_c + Tout_c)/2;

%Defined by inspection (in order to P_N(t=tsim) = 1)
%Tbar_f0 = 1291.9;
%Tbar_c0 = 1030;

%Defined by inspection (in order to P_N(t=tsim) = 1)
%Tbar_f0 = 1288;
%Tbar_c0 = 1060;

%% Temperature profile (01-07-2024)

% PKE + TH models

% Loading steady-state values for temperature profiles (steady-state model)
temp_vars = {'Tc', 'Tc_e', 'Tw', 'Tf', 'Tf_e'};

filename = "SA_DualFluid.mat";

load(filename, temp_vars{:});

%
Tc_0 = Tc; Te_c0 = Tc_e; Tf_0 = Tf; Te_f0 = Tf_e; Tw_0 = Tw;
clear(temp_vars{:});

Tin_f0 = 1300; Tin_c0 = 973;

Tbar_f0 = mean(Te_f0(12), Tin_f0);
Tbar_c0 = mean(Te_c0(12), Tin_c0);