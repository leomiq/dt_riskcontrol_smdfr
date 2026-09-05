function KPI = run_phase2c_full_loop(solver, seed, saveTag)
% RUN_PHASE2C_FULL_LOOP  One full 61-query closed-loop realization for the
% Phase 2C KPI comparison (overshoot, settling time, steady-state error,
% control total variation, min DSM, peak PCT, IAE).
%
%   solver : 'LJ', 'GA' or 'PID'
%   seed   : realization seed controlling the optimizer draws at every query
%            (ignored for the deterministic PID)
%   saveTag: optional string prefix for the saved MAT (default 'phase2c_fullloop')
%
% The loop replicates the archived DET_scheduler_v1_perc_online_LJ/GA flow:
%   - accident stream fixed: rng(1,"twister") before scenario initialization,
%     exactly as archived (sampled_A, B=10, mdot_c_drift=10, vel_ACS=2);
%   - optimizer stream: rng(seed,"twister") immediately before each optimizer
%     call, so every query of a realization uses the same recorded seed;
%   - current DET_2.m semantics (Y_denorm_50 median perturbation path);
%   - no per-query workspace save; one compact MAT per realization.
%
% KPI computation follows Performance_analysis.txt:
%   calculate_performance_metrics(time_snapshots, PCT, control, params) with
%   PCT_ref = target, PCT_max = 1244 K, tolerance = 0.3 K (0.2 K for PID).
%
% Requires (current folder or path): the V5 experiment files, mc_inference.py,
% model_200x200_drop1_1.pth, and a MATLAB Python environment with torch+CUDA.

    if nargin < 3 || isempty(saveTag)
        saveTag = 'phase2c_fullloop';
    end

    %% ---- Fixed experiment parameters (archived primary series) ----
    load('max_min_valuess.mat');          % X_minmax, Y_minmax
    Ncomps = 5;
    N_par_dropout = 59;                   % MC-dropout trajectories
    N_resampling = 2;
    N_branches = 2;
    t_sam = 0.5;
    time_snapshot_init = 1160;
    time_snapshots = time_snapshot_init:t_sam:time_snapshot_init + 30;  % 61 queries
    N_pred = 25;
    % MHD pump saturation limits. The archived bounds [56.69, 56.75] sit only
    % 0.003 .. 0.063 A above the pump's zero-flow current (DeltaP/C = 56.6872 A
    % in the ECM model of CurrentToMassFlow.m), so they cap the coolant at
    % 12.86 kg/s even though the slope is ~205 kg/s per A. The archived PID
    % series already drove this same plant model to 57.26 A, so a wider upper
    % bound stays inside a range the model has been through.
    % Override with the environment variable U_MAX_OVERRIDE.
    u_min = 56.69; u_max = 56.75;
    umax_env = getenv('U_MAX_OVERRIDE');
    if ~isempty(umax_env)
        u_max = str2double(umax_env);
        fprintf('[u_max OVERRIDE] upper bound = %.5f A (default 56.75)\n', u_max);
    end
    u_bounds = repmat([u_min, u_max], N_pred+1, 1);
    % LJ parameters. RETUNED 2026-09-04 when u_max was widened: r_init scales
    % with the range (0.03 -> 0.105 A), so the original gamma = 0.95 left the
    % radius ceiling at r_init*0.95^30 ~ 0.023 A after 30 passes -- about the
    % width of the entire original range, and too coarse for the narrow interior
    % minimum the penalised cost has (~0.05 A wide, cost x179 at its edges).
    % gamma = 0.90 brings that ceiling to ~4e-3 A. R is raised from 15 to 24
    % because the same number of samples now has to cover a 3.5x wider range.
    % Override either from the environment to reproduce the original settings.
    max_passes = 30; M = 5; r_init = 0.5*(u_max - u_min);
    R = 24; gamma = 0.90;
    if ~isempty(getenv('LJ_GAMMA')), gamma = str2double(getenv('LJ_GAMMA')); end
    if ~isempty(getenv('LJ_R')),     R = round(str2double(getenv('LJ_R'))); end
    fprintf('[LJ params] R = %d, gamma = %.3f, r_init = %.5f A, max_passes = %d\n', ...
        R, gamma, r_init, max_passes);
    Tave = 543;
    target = CalculateTs(100, Tave);
    u0 = (u_min + r_init) * ones(N_pred+1, 1)';
    d_mpc = ones(N_pred+1, 1)';
    delayed_time_rom = 10;
    name_hrm = 'SMDFR_hrm_MCDET_2022b_control_CMEX';
    name_rom = 'rom_MCDET_v2_2022b_control_CMEX';
    name_wks = 'simulation_core_dualfluid_FI3_sim_rom_extract_MCDET_v2';

    %% ---- Bookkeeping (as archived) ----
    times_dropout = zeros(1, numel(time_snapshots));
    times_det = zeros(1, numel(time_snapshots));
    times_postproc = zeros(1, numel(time_snapshots));
    u_V2P_trajectory = zeros(size(time_snapshots));
    X_1_trajectory = zeros(1);
    X_2_trajectory = zeros(1);
    control_series = zeros(1,3);
    tw_max_V2P_trajectory = zeros(1,2);
    first_actions = zeros(size(time_snapshots));
    optimizer_runtimes = zeros(size(time_snapshots));

    % ---- Checkpointing ----
    % A power loss during the 2026-09-04 run killed a parallel worker at query
    % 33 of 61 and, because results were saved only at the very end, all 33
    % completed queries were lost. State is now written after every query, and
    % an existing checkpoint is resumed automatically.
    ckptFile = sprintf('%s_%s_seed%d_checkpoint.mat', saveTag, solver, seed);
    start_id = 1;
    if isfile(ckptFile)
        try
            C = load(ckptFile);
            if C.snapshot_id_done >= numel(time_snapshots)
                fprintf('Checkpoint %s is already complete (%d queries).\n', ...
                    ckptFile, C.snapshot_id_done);
            end
            u_V2P_trajectory      = C.u_V2P_trajectory;
            first_actions         = C.first_actions;
            optimizer_runtimes    = C.optimizer_runtimes;
            control_series        = C.control_series;
            tw_max_V2P_trajectory = C.tw_max_V2P_trajectory;
            X_1_trajectory        = C.X_1_trajectory;
            X_2_trajectory        = C.X_2_trajectory;
            times_dropout         = C.times_dropout;
            start_id              = C.snapshot_id_done + 1;
            fprintf('RESUMING from checkpoint %s at query %d/%d (t=%.1f s)\n', ...
                ckptFile, start_id, numel(time_snapshots), time_snapshots(min(start_id,end)));
        catch ME
            fprintf('Checkpoint unreadable (%s) -- starting from query 1\n', ME.message);
            start_id = 1;
        end
    end

    for snapshot_id = start_id:numel(time_snapshots)
        time_snapshot_root = time_snapshots(snapshot_id);
        len_seq = 200;

        %% --- Accident scenario (seed fixed, as archived) ---
        rng(1, "twister");
        load_system(name_hrm);
        simRoot = Simulink.SimulationInput(name_hrm);
        mdlWks = get_param(name_hrm, 'ModelWorkspace');
        unif_dist_drift = getVariable(mdlWks, 'unif_dist_drift');
        B_l = getVariable(mdlWks, 'B_l'); B_s = getVariable(mdlWks, 'B_s');
        A_l = getVariable(mdlWks, 'A_l'); A_s = getVariable(mdlWks, 'A_s');
        vel_ACS_ub = getVariable(mdlWks, 'vel_ACS_ub');
        vel_ACS_lb = getVariable(mdlWks, 'vel_ACS_lb');
        Lambda_t = getVariable(mdlWks, 'Lambda_t');
        Mu_t = getVariable(mdlWks, 'Mu_t');
        sampled_B = 10;
        sampled_A = A_l + (A_s - A_l)*rand();
        sampled_mdot_c_drift = 10;
        sampled_vel_ACS = 2;
        simRoot = setVariable(simRoot, 'B', sampled_B, 'Workspace', name_hrm);
        simRoot = setVariable(simRoot, 'A', sampled_A, 'Workspace', name_hrm);
        simRoot = setVariable(simRoot, 'mdot_c_drift', sampled_mdot_c_drift, 'Workspace', name_hrm);
        simRoot = setVariable(simRoot, 'vel_ACS', sampled_vel_ACS, 'Workspace', name_hrm);

        if snapshot_id == 1
            simRoot = setModelParameter(simRoot, "StopTime", string(time_snapshot_root));
            inport_data = [transpose(0:t_sam:time_snapshot_root), ...
                zeros(1, numel(0:t_sam:time_snapshot_root))', ...
                zeros(1, numel(0:t_sam:time_snapshot_root))'];
            control_series = cat(1, control_series, inport_data(1:end-1,:));
            inport_data_ext = control_series(2:end,:);
            simRoot = simRoot.setExternalInput(inport_data_ext);
        else
            simRoot = setModelParameter(simRoot, "StopTime", string(time_snapshot_root));
            inport_data = [time_snapshots(snapshot_id - 1)', ...
                u_V2P_trajectory(1, snapshot_id - 1)', 1];
            control_series = cat(1, control_series, inport_data);
            inport_data_ext = control_series(2:end,:);
            simRoot = simRoot.setExternalInput(inport_data_ext);
        end
        simRoot = setModelParameter(simRoot, "SaveFinalState", "on", ...
            "SaveOperatingPoint", "on", "FinalStateName", "finalOPRoot");
        outRoot = sim(simRoot);

        simRoot2 = simRoot;
        simRoot2 = setModelParameter(simRoot2, "StopTime", string(time_snapshot_root - delayed_time_rom));
        simRoot2 = setModelParameter(simRoot2, "SaveState", "off", ...
            "SaveFinalState", "on", "SaveOperatingPoint", "on", ...
            "FinalStateName", "finalOPRoot2");
        outRoot2 = sim(simRoot2);
        finalOP2_Raiz = outRoot2.finalOPRoot2;
        close_system(name_hrm);

        if snapshot_id == 1
            X_1 = outRoot.mdot_c_leak_signal.Data(:);
            X_2 = outRoot.Te_c_signal.Data(12,:,:);
            X_1_trajectory = cat(1, X_1_trajectory, X_1);
            X_2_trajectory = cat(1, X_2_trajectory, X_2(:));
            tw_max_V2P_trajectory = cat(1, tw_max_V2P_trajectory, ...
                [outRoot.tout(1:end) outRoot.tw_max_signal.Data(1:end,:)]);
        else
            X_1 = outRoot.mdot_c_leak_signal.Data(end);
            X_2 = outRoot.Te_c_signal.Data(12,:,end);
            X_1_trajectory = cat(1, X_1_trajectory, X_1);
            X_2_trajectory = cat(1, X_2_trajectory, X_2(:));
            tw_max_V2P_trajectory = cat(1, tw_max_V2P_trajectory, ...
                [outRoot.tout(end) outRoot.tw_max_signal.Data(end,:)]);
        end
        finalOP_Raiz = outRoot.finalOPRoot;

        %% --- MC Dropout (GPU) ---
        [Y_denorm, times_dropout] = run_mc_dropout_inference( ...
            X_1_trajectory, X_2_trajectory, len_seq, X_minmax, Y_minmax, ...
            N_par_dropout, snapshot_id, times_dropout);
        [n_trials, len_seq2, n_features] = size(Y_denorm);
        Y_denorm_95 = zeros(len_seq2, n_features);
        Y_denorm_50 = zeros(len_seq2, n_features);
        for feature_idx = 1:n_features
            feature_data = squeeze(Y_denorm(:, :, feature_idx));
            Y_denorm_95(:, feature_idx) = prctile(feature_data, 95, 1)';
            Y_denorm_50(:, feature_idx) = prctile(feature_data, 50, 1)';
        end

        %% --- ROM initial operating point ---
        load_system([name_rom '.slx']);
        simromRoot = Simulink.SimulationInput(name_rom);
        simromRoot = setModelParameter(simromRoot, 'FastRestart', 'off');
        simromRoot = setModelParameter(simromRoot, "SaveFinalState", "on", ...
            "SaveOperatingPoint", "on", "FinalStateName", "finalOP");
        simromRoot = setModelParameter(simromRoot, "StopTime", string(time_snapshot_root - delayed_time_rom));
        inport_data = [transpose(0:t_sam:time_snapshot_root), ...
            outRoot.mdot_c_leak_signal.Data(:), ...
            zeros(1, numel(0:t_sam:time_snapshot_root))', ...
            zeros(1, numel(0:t_sam:time_snapshot_root))', ...
            zeros(1, numel(0:t_sam:time_snapshot_root))'];
        simromRoot = simromRoot.setExternalInput(inport_data);
        simromRoot = setVariable(simromRoot, 'vel_ACS', sampled_vel_ACS, 'Workspace', name_wks);
        outrom = sim(simromRoot);
        romOP = outrom.finalOP;
        % NOTE: the ROM model is intentionally left loaded for DET_2's
        % persistent accelerator session (the archived flow also relies on it).

        %% --- det_params ---
        Y_denorm_95_pred = Y_denorm_95(1:N_pred, :);
        Y_denorm_50_pred = Y_denorm_50(1:N_pred, :);
        tsim_pred = time_snapshot_root + N_pred*t_sam;
        det_params = struct();
        det_params.N_par = 1;
        det_params.M_workers = 8;
        det_params.finalOP_Raiz = finalOP_Raiz;
        det_params.finalOP2_Raiz = finalOP2_Raiz;
        det_params.romOP = romOP;
        det_params.tsim = tsim_pred;
        det_params.time_snapshot_root = time_snapshot_root;
        det_params.t_sam = t_sam;
        det_params.outRoot = outRoot;
        det_params.outrom = outrom;
        det_params.Y_denorm = reshape(Y_denorm(1,1:N_pred,:), [N_pred 2]);
        det_params.Y_denorm_50 = Y_denorm_50_pred;
        det_params.Y_denorm_95 = Y_denorm_95_pred;
        det_params.sampled_vel_ACS = sampled_vel_ACS;
        det_params.Lambda_t = Lambda_t;
        det_params.Mu_t = Mu_t;
        det_params.Ncomps = Ncomps;
        det_params.time_snapshots = time_snapshots;
        det_params.snapshot_id = snapshot_id;
        det_params.N_resampling = N_resampling;
        det_params.N_branches = N_branches;
        det_params.name_root = name_hrm;
        det_params.name_Wks = name_wks;
        det_params.delayed_time_rom = delayed_time_rom;
        det_params.control_series = control_series;

        % Optional DSM safety-margin penalty (run_penalised_trajectory stages
        % these in the environment). They must ride inside det_params: the cost
        % is evaluated inside a parfor in luus_jaakola_det, and client globals
        % do not reach parallel workers.
        dsm_thr_env = getenv('DSM_THR');
        dsm_lam_env = getenv('DSM_LAMBDA');
        if ~isempty(dsm_thr_env) && ~isempty(dsm_lam_env)
            det_params.dsm_thr    = str2double(dsm_thr_env);
            det_params.lambda_dsm = str2double(dsm_lam_env);
            if snapshot_id == 1
                fprintf('[DSM penalty ACTIVE] dsm_thr = %.8f, lambda = %.0e\n', ...
                    det_params.dsm_thr, det_params.lambda_dsm);
            end
        end

        %% --- Controller / optimizer ---
        if strcmpi(solver, 'PID')
            pid_params = struct();
            pid_params.Kp = -0.5; pid_params.Ki = -0.4; pid_params.Kd = 0;
            pid_params.PCT_ref = target;
            pid_params.PCT_max = 1244;
            pid_params.tolerance = 0.2;
            pid_params.I_max = 100;
            pid_params.I_min = 56.68;
            pid_params.dI_max = 60;
            pid_current_state = outRoot.tw_max_signal.Data(end);
            if snapshot_id == 1
                pid_current_control = 0;
            else
                pid_current_control = u_V2P_trajectory(1, snapshot_id-1);
            end
            opt_tic = tic;
            [u_opt, ~] = pid_controller(pid_current_state, pid_current_control, 0.5, pid_params);
            optimizer_runtimes(snapshot_id) = toc(opt_tic);
        else
            rng(seed, "twister");          % recorded optimizer stream per realization
            opt_tic = tic;
            if strcmpi(solver, 'LJ')
                [~, u_opt, ~, ~, tvec] = luus_jaakola_det( ...
                    u0, u_bounds, max_passes, M, r_init, gamma, R, ...
                    det_params, target, d_mpc, time_snapshot_root, delayed_time_rom, t_sam);
            else % GA (equivalent budget: pop 20 x 20 generations = 420 evals)
                [~, u_opt, ~, ~, tvec] = ga_det( ...
                    u0, u_bounds, 20, 20, 0.80, 0.30, ...
                    det_params, target, d_mpc, time_snapshot_root, delayed_time_rom, t_sam);
            end
            optimizer_runtimes(snapshot_id) = sum(tvec);
        end
        u_apply = u_opt(1);
        u_V2P_trajectory(1, snapshot_id) = u_apply;
        first_actions(snapshot_id) = u_apply;
        fprintf('[%s] t=%.1f s  u=%.5f A  (optimizer %.1f s)\n', ...
            solver, time_snapshot_root, u_apply, optimizer_runtimes(snapshot_id));

        % Persist state so a crash or power loss costs one query, not the run.
        snapshot_id_done = snapshot_id; %#ok<NASGU>
        try
            save(ckptFile, 'snapshot_id_done', 'u_V2P_trajectory', 'first_actions', ...
                'optimizer_runtimes', 'control_series', 'tw_max_V2P_trajectory', ...
                'X_1_trajectory', 'X_2_trajectory', 'times_dropout', ...
                'time_snapshots', 'solver', 'seed', '-v7.3');
        catch ME
            fprintf('   !! checkpoint save failed: %s\n', ME.message);
        end
    end

    %% ---- KPI computation (Performance_analysis.txt pipeline) ----
    [ctrl_idx, ~] = find(control_series(:,1) == time_snapshots);
    [v2p_idx, ~] = find(tw_max_V2P_trajectory(:,1) == time_snapshots);
    PCT = tw_max_V2P_trajectory(v2p_idx, 2);
    control = control_series(ctrl_idx, 2);
    params = struct();
    params.PCT_ref = target;
    params.PCT_max = 1244;
    params.tolerance = 0.3;
    if strcmpi(solver, 'PID')
        params.tolerance = 0.2;
    end
    metrics = calculate_performance_metrics(time_snapshots, PCT, control, params);
    metrics.optimizer_runtime_median_s = median(optimizer_runtimes);
    metrics.optimizer_runtime_mean_s = mean(optimizer_runtimes);
    metrics.optimizer_runtime_max_s = max(optimizer_runtimes);

    %% ---- Save compact record ----
    outFile = sprintf('%s_%s_seed%d.mat', saveTag, solver, seed);
    save(outFile, 'metrics', 'solver', 'seed', 'time_snapshots', 'PCT', 'control', ...
        'u_V2P_trajectory', 'first_actions', 'optimizer_runtimes', ...
        'tw_max_V2P_trajectory', 'control_series');
    fprintf('Saved %s\n', outFile);

    KPI = metrics;
end

%% ----------------------------------------------------------------------
function X_norm = norm_function(X_input, X_minmax)
    X_norm = zeros(size(X_input));
    for i = 1:size(X_input,2)
        X_norm(:,i) = (X_input(:,i) - X_minmax(i,1)) / (X_minmax(i,2) - X_minmax(i,1));
    end
end

function Y_denorm = denorm_function(Y_input, Y_minmax)
    Y_denorm = zeros(size(Y_input));
    for i = 1:size(Y_input,3)
        Y_denorm(:,:,i) = Y_input(:,:,i)*(Y_minmax(i,2) - Y_minmax(i,1)) + Y_minmax(i,1);
    end
end

function [Y_denorm, times_dropout] = run_mc_dropout_inference( ...
    X_1, X_2, len_seq, X_minmax, Y_minmax, N_par, snapshot_id, times_dropout)
    X_input = zeros(len_seq, 2);
    X_input(:,1) = X_1(end-len_seq+1:end);
    X_input(:,2) = X_2(end-len_seq+1:end);
    X_norm = norm_function(X_input, X_minmax);
    py_module = py.importlib.import_module('mc_inference');
    X_py = py.numpy.array(X_norm);
    num_trajs = py.int(N_par);
    dropout_time = tic;
    result_py = py_module.mc_dropout_trajectories(py_module.model, X_py, num_trajs);
    times_dropout(snapshot_id) = toc(dropout_time);
    Y_norm = double(py.numpy.array(result_py));
    Y_denorm = denorm_function(Y_norm, Y_minmax);
    clear result_py; clear py_module; clear X_py; clear num_trajs;
end
