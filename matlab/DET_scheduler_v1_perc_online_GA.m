%% Useful files to save the data
%simulation_core_dualfluid_FI3_sim_rom_extract;
%save('simcore_DF_sim_ROM_extract');

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
time_snapshots = time_snapshot_init:t_sam:time_snapshot_init + 30;

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
for snapshot_id = 1:length(time_snapshots)
    %% PO terminal condition (debug)
    % if snapshot_id > 1
    %     if u_apply == 0 && tw_max_V2P_trajectory(snapshot_id-1) <= target
    %         break
    %     end
    % end

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
    name_root = 'SMDFR_hrm_MCDET_2022b_control_CMEX';
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
        inport_data = [time_snapshots(snapshot_id - 1)' ...
                      u_V2P_trajectory(1,snapshot_id - 1)' ...
                      1];
        control_series = cat(1, control_series, inport_data);
        inport_data_ext = control_series(2:end,:);
        simRoot = simRoot.setExternalInput(inport_data_ext);
    end
    
    % Setting intial state of SimulationInput object
    simRoot = setModelParameter(simRoot, "SaveFinalState", "on", "SaveOperatingPoint", "on", "FinalStateName", "finalOPRoot");

    %tic;
    % Running the model
    outRoot = sim(simRoot);
    %toc;

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
    
    %% Preprocessing signals: distribution of signals
    
    % Assuming Y_denorm is a 3D array: (n_trials, sequence_length, n_features)
    % Replace with your actual Y_denorm data
    [n_trials, len_seq, n_features] = size(Y_denorm);
    
    % Initialize the result array
    Y_denorm_95 = zeros(len_seq, n_features);
    
    % Loop over each feature to compute the 95th percentile for each timestep
    for feature_idx = 1:n_features
        % Extract data for the current feature
        feature_data = squeeze(Y_denorm(:, :, feature_idx)); % (n_trials x sequence_length)
    
        % Compute the 95th percentile for each timestep
        Y_denorm_95(:, feature_idx) = prctile(feature_data, 95, 1)'; % (sequence_length x 1)
    end
    
    %% ROM simulation (just for getting the romOP at the beginning)
    
    %if snapshot_id == 1
        %Ncomps = 5;
        name_root = 'rom_MCDET_v2_2022b_control_CMEX';
        name_Wks = 'simulation_core_dualfluid_FI3_sim_rom_extract_MCDET_v2';
        %save_system('rom_MCDET_v2_2022b_control_CMEX.slx');
        load_system('rom_MCDET_v2_2022b_control_CMEX.slx');
        simromRoot = Simulink.SimulationInput(name_root);
        simromRoot = setModelParameter(simromRoot,'FastRestart', 'off');
        simromRoot = setModelParameter(simromRoot,"SaveFinalState", "on", "SaveOperatingPoint", "on", "FinalStateName", "finalOP");

        %tsim = time_snapshot_root + length_seq*t_sam;
        %predicted window:seconds
        simromRoot = setModelParameter(simromRoot, "StopTime", string(time_snapshot_root - delayed_time_rom));
        
        % if snapshot_id == 1
        % Combine into a single input matrix (for use in Inport). 
        %1: mdot_c_leak; 2:tw_max; 3:control signal (ACS)
        inport_data = [transpose(0:t_sam:time_snapshot_root), outRoot.mdot_c_leak_signal.Data(:),...
                        zeros(1,length(0:t_sam:time_snapshot_root))', zeros(1,length(0:t_sam:time_snapshot_root))',...
                        zeros(1,length(0:t_sam:time_snapshot_root))'];
        simromRoot = simromRoot.setExternalInput(inport_data);
        % else
        %     simromRoot = setInitialState(simromRoot, romOPs(snapshot_id-1));
        %     %1: mdot_c_leak; 2:tw_max; 3:control signal (ACS)
        %     inport_data = [transpose(time_snapshot_root - t_sam:t_sam:time_snapshot_root),...
        %                     outRoot.mdot_c_leak_signal.Data(:),...
        %                     [0 0]', [0 0]',...
        %                     [0 0]'];
        %     simromRoot = simromRoot.setExternalInput(inport_data);
        % end

        % Get target distributions parameters
        mdl_Wks = get_param(name_root,'ModelWorkspace');
        %Lambda_t = getVariable(mdl_Wks,'Lambda_t');
        %Mu_t = getVariable(mdl_Wks,'Mu_t');
        % set sampled ACS velocity
        simromRoot = setVariable(simromRoot,'vel_ACS', sampled_vel_ACS, 'Workspace', name_Wks);
        
        outrom = sim(simromRoot);
        
        close_system('rom_MCDET_v2_2022b_control_CMEX.slx');
        
        romOP = outrom.finalOP;

        romOPs{end+1} = romOP;
        
        disp('sim rom done.');
        
        %save('OPs.mat')
    %end
    
    %% Luus-Jaakola optimizer and DET (Simulink-based)
    %rng(42,"twister");
    %init_loop = tic; 
    
    % MPC Parameters
    %N_pred can be different than the predicted window length of GBM
    N_pred = 25;
    t_sam = 0.5;
    x0 = [0];
    
    % Control parameters
    u_max = 56.75; u_min = 56.68;
    u_bounds = repmat([u_min, u_max], N_pred+1, 1);
    
    % Luus-Jaakola parameters
    max_passes = 30;
    M = 5;
    R = 15;
    r_init = 0.5*(u_max - u_min);
    gamma = 0.95;
    Tave = 543;
    %target = 1144;
    target = CalculateTs(100,543);
    
    % Initial guess
    u0 = (u_min+r_init)*ones(N_pred+1, 1)';
    
    % vector of discrete mpc vector
    d_mpc = ones(N_pred+1, 1)';
    
    % Additional parameters for DET
    N_par                         = 1;
    M_workers                     = 8;
    
    % Adapting parameters for DET (prediction length horizon)
    Y_denorm_95_pred              = Y_denorm_95(1:N_pred,:);
    tsim_pred                     = time_snapshot_root + N_pred*t_sam;
    
    % DET parameters struct setup (populate with actual simulation parameters)
    det_params = struct();
    
    % Explicit variable assignment from base workspace
    % Standalone struct creation from base workspace variables
    det_params.N_par              = N_par;
    det_params.M_workers          = M_workers;
    det_params.finalOP_Raiz       = finalOP_Raiz;
    det_params.finalOP2_Raiz      = finalOP2_Raiz;
    det_params.romOP              = romOP;
    det_params.tsim               = tsim_pred;
    det_params.time_snapshot_root = time_snapshot_root;
    det_params.t_sam              = t_sam;
    det_params.outRoot            = outRoot;
    det_params.outrom             = outrom;
    det_params.Y_denorm           = reshape(Y_denorm(1,1:N_pred,:),[N_pred 2]);
    det_params.Y_denorm_95        = Y_denorm_95_pred;
    det_params.sampled_vel_ACS    = sampled_vel_ACS;
    det_params.Lambda_t           = Lambda_t;
    det_params.Mu_t               = Mu_t;
    det_params.Ncomps             = Ncomps;
    det_params.time_snapshots     = time_snapshots;
    det_params.snapshot_id        = snapshot_id;
    det_params.N_resampling       = N_resampling;
    det_params.N_branches         = N_branches;
    det_params.name_root          = name_root;
    det_params.name_Wks           = name_Wks;
    det_params.delayed_time_rom   = delayed_time_rom;
    det_params.control_series     = control_series;
    det_params.snapshot_id        = snapshot_id;
    
    %% Testing DET
    [forest_data, forest_roots] = DET_2(det_params, 56.69*ones(1,N_pred+1), d_mpc);
    
    %%
    
    %Run optimization (CHECK INTERNAL PARAMETERS!!)
    %[Tp_opt, u_opt, f_opt, history, t_vector] = luus_jaakola_det(u0, u_bounds, max_passes, M, r_init, gamma, R, det_params, target, d_mpc,...
    %                                           time_snapshot_root, delayed_time_rom, t_sam);
    
    %% Genetic Algorithm

    instantiate_ga_params;

    % --- Run GA_DET ---
    [Tp_opt, u_opt, f_opt, history, t_vector] = ga_det( ...
        u0, bounds, num_generations, pop_size, ...
        crossover_frac, mutation_rate, ...
        det_params, target, d_mpc, ...
        time_snapshot_root, delayed_time_rom, t_sam);

   
    %%
    % Apply first control input
    %u_apply = u_opt(1);

    %% PID controller

    % % PID parameters setup
    % pid_params = struct();
    % pid_params.Kp = -0.5;   % Proportional gain
    % pid_params.Ki = -0.4;  % Integral gain  
    % pid_params.Kd = 0;   % Derivative gain
    % pid_params.PCT_ref = target; %setpoint, considering the opening time
    % pid_params.PCT_max = 1244;  % PCT MAX for performance tests
    % pid_params.tolerance = 0.2; % tolerance of the target value, for performance tests
    % pid_params.I_max = 100;
    % pid_params.I_min = 56.68;
    % pid_params.dI_max = 60;
    % pid_current_state = outRoot.tw_max_signal.Data(end);
    % if snapshot_id == 1
    %     pid_current_control = 0;
    % else 
    %     pid_current_control = u_V2P_trajectory(1,snapshot_id-1);
    % end
    % [u_opt, pid_info] = pid_controller(pid_current_state, pid_current_control, 0.5, pid_params);

    %%
    % Apply first control input
    u_apply = u_opt(1);

    %%
    % Display last temperature
    %disp('Last T_p');
    %disp(pid_current_state)
    % Display result
    disp('First optimal control input:');
    disp(u_apply);
    disp('m_ACS');
    disp(CurrentToMassFlow(u_apply));
    
    %% check
    %[forest_data_opt, forest_roots_opt] = DET_2(det_params, u_opt*ones(1,N_pred+1), d_mpc);
    [forest_data_opt, forest_roots_opt] = DET_2(det_params, u_opt, d_mpc);
    
    %%
    [forest_data_best, forest_roots_best] = DET_2(det_params, u_max*ones(1,N_pred+1), d_mpc);
    
    %% Call postproc routine (ISMCDET)
    % init_postproc = tic;
    % run("DET_scheduler_postprocessing_perc_online.m")
    % end_postproc = toc(init_postproc);
    % times_postproc(snapshot_id) = end_postproc;
    % %storing SMs for each snapshot to track evolution in time
    % M_U_ya_IS_vec(snapshot_id) = M_U_ya_IS;
    % M_U_ya_MC_vec(snapshot_id) = M_U_ya_MC;

    %% Keeping track of control trajectory
    u_V2P_trajectory(1,snapshot_id) = u_apply;

    %% Save global workspace
    save(['GA_test_sampled_vel_ACS_2_loop_' num2str(time_snapshot_root) '.mat'])
end

%% Performance metrics

GA_time = time_snapshots;
[control_series_performance_idx, col] = find(control_series(:,1) == time_snapshots);
[V2P_performance_idx, col] = find(tw_max_V2P_trajectory(:,1) == time_snapshots);
GA_PCT = tw_max_V2P_trajectory(V2P_performance_idx,2);
GA_control = control_series(control_series_performance_idx,2);

GA_params = struct();
GA_params.PCT_ref = target; %setpoint, considering the opening time
GA_params.PCT_max = 1244;  % PCT MAX for performance tests
GA_params.tolerance = 0.3; % tolerance of the target value, for performance tests

metrics_pid = calculate_performance_metrics(GA_time, GA_PCT, GA_control, GA_params);

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
    % py.importlib.reload(py_module);

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
