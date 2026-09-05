%% Useful files to save the data
%simulation_core_dualfluid_FI3_sim_rom_extract;
%save('simcore_DF_sim_ROM_extract');

% Example paths (replace with actual paths used previously):
setenv('LD_LIBRARY_PATH', '/g100_work/PROJECTS/spack/v0.17/prod/0.17.1/install/0.17/linux-centos8-skylake_avx512/gcc-8.4.1/cuda-11.5.0-ktwkkqqhebhe64r4ial5g632vefweb4i/lib64:/cineca/prod/opt/tools/matlab/r2022b/lib:/g100_work/PROJECTS/spack/v0.17/prod/0.17.1/install/0.17/linux-centos8-cascadelake/gcc-10.2.0/python-3.8.12-doaq4i2pjuhhbrv7qlkeavdeglbce7n7/lib');

%% Parallel trees
Ncomps = 5;
%parpool;
N_par = 59; M_workers = 59;

%% Resampling: parameters
% Number of branches for each resampling
N_resampling = 2;
N_branches = 2;

%% importing maximum and minimum values for BiLSTM model ((de)normalization)

load('max_min_valuess.mat')

%%
t_sam = 0.5;
time_snapshot_init = 1160;
time_snapshots = time_snapshot_init:t_sam:time_snapshot_init + 20;

times_dropout = zeros(1,size(time_snapshots,2));
times_det = zeros(1,size(time_snapshots,2));
times_postproc = zeros(1,size(time_snapshots,2));
snapshots_data = {};

%postproc, IS and MCDET PSMs
%M_U_ya_IS_vec = zeros(1,size(time_snapshots,2));
%M_U_ya_MC_vec = zeros(1,size(time_snapshots,2));

finalOP_Raices = {};
romOPs = {};

% Taking final V2P trajectory
u_V2P_trajectory = zeros(size(time_snapshots));

% Taking sliding windows for MC Dropout
X_1_trajectory = zeros(1);
X_2_trajectory = zeros(1);

%time, u_mpc, active_mpc
control_series = zeros(1,3);

%time, tw_max_signal
tw_max_V2P_trajectory = zeros(1,2);

snapshot_id = 1;

%%
%for snapshot_id = 1:size(time_snapshots,2)
% A set of (independent) trees is...
forest_data = {};
forest_roots = {};
forest_roots_hrm = {};
time_snapshot_root = time_snapshots(snapshot_id);

%% Running the "ground truth" simulation
% At each time step, we will compute this MCDET
len_seq = 200;
% sampling time
t_sam = 0.5;
%mission time (forecast)
tsim = time_snapshot_root + len_seq*t_sam;

%% Loading LOCA initial conditions (HRM)

%fixing the seed (accident) for reproducibility
rng(1,"twister");

% loading variables of interest: coolant mass flow rate (leakage)
%Ncomps = 5;
name_root = 'SMDFR_hrm_MCDET_2022b_control';
name_Wks = 'simulation_core_dualfluid_FI3_sim_hrm_extract_MCDET';

load_system(name_root);
simRoot = Simulink.SimulationInput(name_root);

% Getting the variables from the model workspace
mdlWks = get_param(name_root,'ModelWorkspace');
unif_dist_drift = getVariable(mdlWks,'unif_dist_drift');

% Getting sampling parameters for B (LOCA nonlinear dynamics)
B_l = getVariable(mdlWks,'B_l');
B_s = getVariable(mdlWks, 'B_s');
%sampled_B = B_l + (B_s - B_l)*rand();
sampled_B = 10;

% Getting sampling parameters for B (LOCA nonlinear dynamics)
A_l = getVariable(mdlWks,'A_l');
A_s = getVariable(mdlWks, 'A_s');
sampled_A = A_l + (A_s - A_l)*rand();

% Setting a new sample for mdot_c_drift
%sampled_mdot_c_drift = unif_dist_drift(1,2) + (unif_dist_drift(2,2) - unif_dist_drift(1,2))*rand();
sampled_mdot_c_drift = 10;
simRoot = setVariable(simRoot, 'B', sampled_B, 'Workspace', name_root);
simRoot = setVariable(simRoot, 'A', sampled_A, 'Workspace', name_root);
simRoot = setVariable(simRoot, 'mdot_c_drift', sampled_mdot_c_drift, 'Workspace', name_root);
%simRoot = setVariable(simRoot, 'A', 0.18*sampled_mdot_c_drift, 'Workspace', name_root);
%disp(string(sampled_mdot_c_drift));

% Getting the variables from the model workspace
vel_ACS_ub = getVariable(mdlWks,'vel_ACS_ub');
vel_ACS_lb = getVariable(mdlWks,'vel_ACS_lb');

% Setting a new sample for mdot_c_drift
%sampled_vel_ACS = vel_ACS_lb + (vel_ACS_ub - vel_ACS_lb)*rand();
%disp(string(sampled_vel_ACS))
sampled_vel_ACS = 2;
simRoot = setVariable(simRoot,'vel_ACS', sampled_vel_ACS, 'Workspace', name_root);

% Lambda_t / Mu_t: TARGET failure and repair rates
Lambda_t = getVariable(mdlWks,'Lambda_t');
Mu_t = getVariable(mdlWks,'Mu_t');

% Including the initial parent last state to start the simulation
%parent_last_state = [0 1];

% Feedback loop effect
if snapshot_id == 1
    %simRoot = setInitialState(simRoot, finalOP_Raices(snapshot_id-1));
    simRoot = setModelParameter(simRoot, "StopTime", string(time_snapshot_root));
    inport_data = [transpose(0:t_sam:time_snapshot_root), ...
                zeros(1,length(0:t_sam:time_snapshot_root))', ...
                zeros(1,length(0:t_sam:time_snapshot_root))'];
    control_series = cat(1, control_series, inport_data(1:end-1,:));
    inport_data_ext = control_series(2:end,:);
    simRoot = simRoot.setExternalInput(inport_data_ext);
else
    % time at snapshot_id
    %simRoot = setInitialState(simRoot, finalOP_Raices(snapshot_id-1));
    simRoot = setModelParameter(simRoot, "StopTime", string(time_snapshot_root));
    inport_data = [time_snapshots(snapshot_id-1)' ...
                  u_V2P_trajectory(1,snapshot_id - 1)' ...
                  1];
    control_series = cat(1, control_series, inport_data);
    inport_data_ext = control_series(2:end,:);
    simRoot = simRoot.setExternalInput(inport_data_ext);
end

% Setting intial state of SimulationInput object
simRoot = setModelParameter(simRoot, "SaveFinalState", "on", "SaveOperatingPoint", "on", "FinalStateName", "finalOPRoot");

% Running the model
outRoot = sim(simRoot);

% Running hrm2 model for initial OP of ROM
delayed_time_rom = 10;
simRoot2 = simRoot;
simRoot2 = setModelParameter(simRoot2, "StopTime", string(time_snapshot_root-delayed_time_rom));
simRoot2 = setModelParameter(simRoot2, "SaveState", "off", "SaveFinalState", "on", "SaveOperatingPoint", "on", "FinalStateName", "finalOPRoot2");
outRoot2 = sim(simRoot2);
finalOP2_Raiz = outRoot2.finalOPRoot2;

close_system(name_root);

disp('sim root done.');

% storing information for MC dropout and control series
if snapshot_id == 1
    X_1 = outRoot.mdot_c_leak_signal.Data(:);
    X_2 = outRoot.Te_c_signal.Data(12,:,:);
    X_1_trajectory = cat(1,X_1_trajectory, X_1);
    X_2_trajectory = cat(1,X_2_trajectory, X_2(:));

    %time, u_mpc, active_mpc...
    %control_series = cat(1, control_series,...
    %   [outRoot.tout(1:end-1) zeros(length(outRoot.tout)-1,1) zeros(length(outRoot.tout)-1,1)]);

    %saving final tw_max trajectory (V2P effect)
    tw_max_V2P_trajectory = cat(1, tw_max_V2P_trajectory,...
        [outRoot.tout(1:end) outRoot.tw_max_signal.Data(1:end,:)]);
    
else
    X_1 = outRoot.mdot_c_leak_signal.Data(end);
    X_2 = outRoot.Te_c_signal.Data(12,:,end);
    X_1_trajectory = cat(1,X_1_trajectory, X_1);
    X_2_trajectory = cat(1,X_2_trajectory, X_2(:));

    %time, u_mpc, active_mpc...
    %control_series = cat(1, control_series,...
    %   [outRoot.tout(end) u_V2P_trajectory(1,snapshot_id - 1) 1]);

    %saving final tw_max trajectory (V2P effect)
    tw_max_V2P_trajectory = cat(1, tw_max_V2P_trajectory,...
        [outRoot.tout(end) outRoot.tw_max_signal.Data(end,:)]);
end

%Saving the data
finalOP_Raiz = outRoot.finalOPRoot;

finalOP_Raices{end+1} = finalOP_Raiz;

% data of the root. struct:{'id': time snapshot integer, 'data':out type, 'final_OP': OP dtype, 'init_state': list}
init_state = finalOP_Raiz.loggedStates{13}.Values;
forest_roots_hrm{1} = struct('id', time_snapshot_root, 'data', outRoot, 'final_OPRoot', finalOP_Raiz, 'init_state', init_state, 'Lambda_t', Lambda_t, 'Mu_t', Mu_t);

%% MC Dropout

% (Assuming outRoot, len_seq, X_minmax, N_par, snapshot_id, times_dropout are defined)
[Y_denorm, times_dropout] = run_mc_dropout_inference(...
    X_1_trajectory, X_2_trajectory, len_seq, X_minmax, Y_minmax, N_par, snapshot_id, times_dropout);

disp('mc dropout done.');
    
%% FUNCTIONS

function X_norm = norm_function(X_input, X_minmax)
    %copy X_input
    X_norm = zeros(size(X_input));
    %dim(X) = [n_samples,sequence_lenght, n_features]
    for i = 1:size(X_input,2)
        X_norm(:,i) = (X_input(:,i) - X_minmax(i,1))/(X_minmax(i,2) - X_minmax(i,1));
    end
end

function Y_denorm = denorm_function(Y_input, Y_minmax)
    %copy Y_input
    Y_denorm = zeros(size(Y_input));
    %dim(Y) = [n_samples,sequence_lenght, n_features]
    for i = 1:size(Y_input,3)
        Y_denorm(:,:,i) = Y_input(:,:,i)*(Y_minmax(i,2) - Y_minmax(i,1)) + Y_minmax(i,1);
    end
end

%%
function [Y_denorm, times_dropout] = run_mc_dropout_inference(...
    X_1, X_2, len_seq, X_minmax, Y_minmax, N_par, snapshot_id, times_dropout)
% RUN_MC_DROPOUT_INFERENCE  Performs MC dropout inference using mc_inference.py
% 
% Syntax:
%   [Y_denorm, times_dropout] = run_mc_dropout_inference(...
%       outRoot, len_seq, X_minmax, N_par, snapshot_id, times_dropout)
%
% Inputs:
%   outRoot       - Struct or object containing time-series data
%   len_seq       - Integer length of the sequence to extract
%   X_minmax      - Normalization min/max data for input signals
%   N_par         - Number of trajectories for MC dropout
%   snapshot_id   - Index or identifier for storing the inference time
%   times_dropout - Array for storing MC dropout timings
%
% Outputs:
%   Y_denorm      - Denormalized output from MC dropout inference
%   times_dropout - Updated array with new timing value inserted

    X_input = zeros(len_seq, 2);
    X_input(:,1) = X_1(end-len_seq+1:end);
    X_input(:,2) = X_2(end-len_seq+1:end);

    % 2) Normalize
    X_norm = norm_function(X_input, X_minmax);

    % 3) Import the Python module
    py_module = py.importlib.import_module('mc_inference');
    % If you need to ensure the module reloads each time:
    py.importlib.reload(py_module);

    % 4) Convert input to Python objects
    X_py = py.numpy.array(X_norm);
    num_trajs = py.int(N_par);  % e.g., number of trajectories = 3

    % 5) Run MC dropout inference and measure timing
    dropout_time = tic;
    result_py = py_module.mc_dropout_trajectories(py_module.model, X_py, num_trajs);
    dropout_finish = toc(dropout_time);

    % 6) Store the elapsed time
    times_dropout(snapshot_id) = dropout_finish;

    % 7) Convert the Python result to a MATLAB double
    Y_norm = double(py.numpy.array(result_py));

    % 8) Denormalize
    Y_denorm = denorm_function(Y_norm, Y_minmax);

    % 9) Clear Python references to avoid serialization issues
    clear result_py;
    clear py_module;
    clear X_py; 
    clear num_trajs;
end
