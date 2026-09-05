% Definition of parameters an values
parameters_dualfluid_rom_v2

%% MC simulation - m
% HYPERPARAMETERS - SIMULATION
%nsim: taken from the HRM
nsim = 1;
tsim = 1500;
fault_time = 500;
step = 1;
signals_md = {'mdot_f_signal', 'mdot_c_signal', 'rhoex_signal', 'M_f_signal', 'M_c_signal', 'tf_signal', 'tw_signal', 'tc_signal', 'totalpower_signal'};
%mat_seeds: TAKEN FROM THE HRM
mat_seeds = {};
mat_seeds{1}= [4.000 0 0; -0.6470 -0.6215 -0.0001; -0.0016 -0.0016 -0.0000];
%mat_seeds = load('Casesim-HRM_TEST2_400-800.mat', 'mat_seeds').mat_seeds;

%Fault mode cases (case /in {0,..n} where 0 is default mode (no failure))
% Fuel mass flow rate (mdot_f)
%Fault case: taken from the HRM
fault_cases = {mat_seeds{1}(1,:)};
setup_mat = zeros(9,3);

% BIAS AND DRIFT - INITIAL VALUES FOR THE SIMULATION
mdot_f_value = mdot_f(1);    
mdot_f_delta = 0;           %mdot_f is a vector!
mdot_f_drift = 0;
mdot_f_randu = 0;
mdot_f_rands = (0.05)*0.001*M_f(1);
%!!!
M_f_faultcase = 1;

mdot_c_value = mdot_c(1);    %mdot_c is a vector!
mdot_c_delta = 0;
mdot_c_drift = 0;
mdot_c_randu = 0;
mdot_c_rands = (0.05)*0.001*M_c(1);
%!!!
M_c_faultcase = 0;

% external reactivity
rhoex_value = 0;    
rhoex_bias = 0;
rhoex_drift = 0;
rhoex_randu = 0;
rhoex_rands = (0.05)*rhoex_value;

%%
% MEMORY - FAULT GENERATION SEEDS

% random generation values - parameters
% nominal value of the forcing parameter
%nom_values = [0, 0, rho_0];
%rand_perc = 0.3;

% parameters for the generation of random numbers (uniform dist)
% bias case
%mdot_f_s = -(min(M_f)/(tsim - fault_time))*12;
%mdot_c_s = -(min(M_c)/(tsim - fault_time))*12;
%unif_dist_bias = [mdot_f_s, mdot_c_s, -0.1*rho_0; 0, 0, 0.1*rho_0];
% drift case
%mdot_drift_f_s = -abs(mdot_f_s/(tsim - fault_time));
%mdot_drift_c_s = -abs(mdot_c_s/(tsim - fault_time));
%unif_dist_drift = [mdot_drift_f_s, mdot_drift_c_s, -0.01*rho_0; 0, 0, 0.01*rho_0];

%data for ann
% t, mdot_f, mdot_c, rho_ex, M_f1, M_f2, M_c1, M_c2, tf, tw, tc, P_N
ann_data = {};
% t, mdot_f, mdot_c, rho_ex, M_f1, M_f2, M_c1, M_c2, t_f, t_w, t_c, P_N
ann_data_deltas = zeros(nsim, 12);
% storing the computing times
comp_times_hrm = zeros(nsim);

%data for PARAMETERS, i.e.: mdot_f, mdot_c, rho_ex, M_f1, M_f2, M_c1, M_c2
par_signals = {};
for k = 1:7
    par_signals{k} = {};
end

%setup matrix
%setup_mat(1,:) = fault_cases{1}';
%setup_mat(2,:) = nom_values';
%setup_mat(3:end,1) = [unif_dist_bias(:,1)', unif_dist_drift(:,1)', mdot_f_randu, mdot_f_rands, M_f_faultcase]';
%setup_mat(3:end,2) = [unif_dist_bias(:,2)', unif_dist_drift(:,2)', mdot_c_randu, mdot_c_rands, M_c_faultcase]';
%setup_mat(3:end,3) = [unif_dist_bias(:,3)', unif_dist_drift(:,3)', rhoex_randu, rhoex_rands, 0]';

%% Run ROM simulation parameters

% Case: loss-of-flow accident (fuel circuit)
% Option 1 (by defect). Fault injection: drift with random noise
% mdot_f_drift = -(unif_dist_drift(1,1) / 5);
% mdot_f_drift_sigma = abs(mdot_f_drift) * 0.7;
% caso = fault_cases{1};
% mdot_f_faultcase = caso(1);
% mdot_c_faultcase = caso(2);
% rhoex_faultcase = caso(3);
% Case: loss-of-flow accident (coolant circuit)
% Option 2. Fault injection: drift with random noise
%mdot_f_drift = (unif_dist_drift(1,1) / 5);
%mdot_f_drift_sigma = abs(mdot_f_drift) * 0.3;
%mdot_c_drift = (unif_dist_drift(1,1) / 5);
%mdot_c_drift_sigma = abs(mdot_c_drift) * 0.3;
%caso = fault_cases{1};
%mdot_f_faultcase = 0; M_f_faultcase = 0;
%mdot_c_faultcase = 4; M_c_faultcase = 1;
%rhoex_faultcase = 0;

% Option 2. Fault injection: bias with random noise
% mdot_f_randu = mdot_f_s(1,1) / 10;
% mdot_f_delta = mdot_f_randu;
% mdot_f_rands = abs(mdot_f_randu) * 0.7;
% Option 3 (experiment with step reactivity insertion 100 pcm)
% caso = [0 0 1];
% mdot_f_faultcase = caso(1);
% mdot_c_faultcase = caso(2);
% rhoex_faultcase = caso(3);
% rhoex_value = 0;
% rhoex_bias = 100E-5;
% mdot_f_drift_sigma = 0;

%% Sampling leakage rate (May 2024)
% LOCA: step for leakage rate, connected with coolant mass
%mdot_drift_f_s = abs(mdot_f(1)/(fault_time));
%mdot_drift_f_s = 0;
%mdot_drift_c_s = abs(mdot_c(1)/(1500 - fault_time)); 
unif_dist_drift = [0, 8, -0.01*rho_0; 0, 10, 0.01*rho_0];

mdot_c_drift = 10;
%mdot_c_drift_sigma = mdot_c_drift*0.3;

vel_ACS_ub = 0.09;
vel_ACS_lb = 0.11;
vel_ACS = 0.1;

%% MC Dynamic Event Tree parameters - Auxiliary Cooling System
%   - Set simulation clock
%   - Set event list
%   obs: consider logic rules to open the branches
Ncomps = 5;
% Nstates_max: max number of states, per component
Nstates = [2 2 2 2 2];
% mat_states_index: Ncomps x Nstates
mat_states_index = [1 2; 3 4; 5 6; 7 8; 9 10];

% lambda [occ/year]; mu [year]
year_to_second = 365*24*3600;

%% Importing data for FT and MC direct sampling
load('MCDET_conferencepaper_database.mat')

%%

% Lambda_s/Mu_s: distributional parameters for (IMPORTANCE) SAMPLING. PDF
%Lambda_s = Lambda_t*(1E5);
%Mu_s = Mu_t*(1E4);

% Vact 14-05 - works
Lambda_s = Lambda_t;%*(1E3);
Mu_s = Mu_t;%*(5E3);

% Vact 14-05 - works
%Lambda_s = Lambda_t*(1E5);
%Mu_s = Mu_t*(3E3);

%%

Lambda_s_resample_lb = Lambda_s*0.90;
Lambda_s_resample_ub = Lambda_s*1.10;

Mu_s_resample_lb = Mu_s*0.90;
Mu_s_resample_ub = Mu_s*1.10;

cell_trans = {};
for i = 1:Ncomps
    cell_trans{i} = [0 Lambda_s(i); Mu_s(i) 0];
end

mat_trans = cell2mat(cell_trans);
mat_trans = reshape(mat_trans,2,2,5);
mat_trans = reshape(mat_trans, 1, 20);

% Faulty states: (i.e., the minimal cut sets)
failed_states = [2 1 1 1 1; 1 2 1 1 1; 1 1 2 1 1; 1 1 1 2 1; 1 1 1 1 2];

% Initial states
initial_state = [1 1 1 1 1];
parent_last_state = [0 1];

%% Branching
initial_branch_state = [0 1];
parent_last_time = -0.1;
complementary_branch_flag = 0;

%% Updates (Dec 19, 2024)
% Loca modelled as in Park et al. (2021)

% B: Parameter controlling the LOCA nonlinear dynamics
B = 12;

% Sampling parameter B from a uniform distribution
B_l = 9;
B_s = 12;

% A: Shape parameter
A = 0.3;

% Sampling parameter A from a uniform distribution
A_l = 0.2;
A_s = 0.4;

%% VARIABLES FOR INTEGRATION WITH GBM BASED PREDICTIONS

length_seq = 50; t_sam = 0.5;

%timeseries(x(t),t)
pred_mdot_c_leak = timeseries(zeros(1,50), 0.5:0.5:25);
pred_delta_tw_max = timeseries(zeros(1,50), 0.5:0.5:25);

%% POLYFIT, DELTA TP AGAINST MDOT C
load("fit_poly.mat");
poly_vector = coeffvalues(fit_poly);
poly_vector(end) = poly_vector(end) - 2.5;