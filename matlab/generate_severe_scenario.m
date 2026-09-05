function G = generate_severe_scenario(cfg)
% GENERATE_SEVERE_SCENARIO  Build a severe LOCA case from the HRM, run the
% BiLSTM MC-dropout inference on it, and assemble a det_params ready for the
% DET / DSM-penalty study.
%
% WHY THIS IS NEEDED. Varying ROM parameters cannot make the scenario severe:
% `vel_ACS` is an HRM scenario parameter and the ROM receives the ACS action
% through its `u_ACS` inport (verified -- a 50x degradation of vel_ACS written
% straight into the ROM workspace changed the predicted trajectory by nothing at
% all, and no ROM block references vel_ACS). The accident severity lives in the
% HRM, and reaches the ROM only through the network predictions
% (`pred_mdot_c_leak`, `pred_delta_tw_max`). So a severe case must be generated
% at the source: HRM -> MC-dropout inference -> DET.
%
% PHYSICAL BOUNDS, from simulation_core_dualfluid_FI3_sim_hrm_extract_MCDET.m:
%   A (LOCA shape parameter)      : A_l = 0.2  .. A_s = 0.4
%   B (LOCA nonlinear dynamics)   : B_l = 9    .. B_s = 15
%   mdot_c_drift                  : unif_dist_drift column 2 spans 8 .. 10
%   vel_ACS                       : 0.09 .. 0.11 (nominal 0.1)
% The nominal archived case sits low in those ranges (A = 0.2834, B = 10,
% mdot_c_drift = 10, vel_ACS = 2 declared in the workspace but 0.1 in the model).
% The severe case takes A and B to their upper bounds and degrades vel_ACS to its
% lower bound -- all inside the declared distributions, so the scenario stays
% physically admissible rather than fabricated.
%
% cfg fields (all optional):
%   A            LOCA shape parameter          (default 0.4  = A_s)
%   B            LOCA dynamics parameter       (default 15   = B_s)
%   mdot_c_drift coolant drift                 (default 10)
%   vel_ACS      ACS velocity in the HRM       (default 0.09 = vel_ACS_ub)
%   t_query      query instant [s]             (default 1180)
%   delayed      ROM anticipation [s]          (default 10)
%   t_sam        sample time [s]               (default 0.5)
%   N_pred       prediction horizon steps      (default 25)
%   N_dropout    MC-dropout trajectories       (default 59)
%   tag          output file tag               (default 'severe')

    if nargin < 1 || isempty(cfg), cfg = struct(); end
    def = struct('A', 0.4, 'B', 15, 'mdot_c_drift', 10, 'vel_ACS', 0.09, ...
                 't_query', 1180, 'delayed', 10, 't_sam', 0.5, ...
                 'N_pred', 25, 'N_dropout', 59, 'tag', 'severe');
    fn = fieldnames(def);
    for i = 1:numel(fn)
        if ~isfield(cfg, fn{i}), cfg.(fn{i}) = def.(fn{i}); end
    end

    G = struct();
    G.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    G.cfg = cfg;
    fprintf('\n===== Severe-scenario generation (HRM -> BiLSTM -> DET) =====\n');
    fprintf('Started %s\n\n', G.timestamp);
    fprintf('LOCA A            = %.4f   (bounds 0.2 .. 0.4)\n', cfg.A);
    fprintf('LOCA B            = %.4f   (bounds 9 .. 15)\n', cfg.B);
    fprintf('mdot_c_drift      = %.4f\n', cfg.mdot_c_drift);
    fprintf('vel_ACS (HRM)     = %.4f   (bounds 0.09 .. 0.11)\n', cfg.vel_ACS);
    fprintf('query instant     = %.1f s, delayed %.1f s, t_sam %.2f s\n', ...
        cfg.t_query, cfg.delayed, cfg.t_sam);
    fprintf('N_pred %d, N_dropout %d\n\n', cfg.N_pred, cfg.N_dropout);

    name_hrm = 'SMDFR_hrm_MCDET_2022b_control_CMEX';
    name_rom = 'rom_MCDET_v2_2022b_control_CMEX';
    name_wks = 'simulation_core_dualfluid_FI3_sim_rom_extract_MCDET_v2';
    t_q     = cfg.t_query;
    t_sam   = cfg.t_sam;
    delayed = cfg.delayed;
    tsim    = t_q + cfg.N_pred * t_sam;

    load('max_min_valuess.mat', 'X_minmax', 'Y_minmax');

    %% ---- Stage 1: HRM accident up to the query instant ----
    fprintf('-- Stage 1: HRM accident simulation --\n');
    rng(1, 'twister');                     % same accident stream as the archive
    load_system(name_hrm);
    mdlWks = get_param(name_hrm, 'ModelWorkspace');
    Lambda_t = getVariable(mdlWks, 'Lambda_t');
    Mu_t     = getVariable(mdlWks, 'Mu_t');

    simRoot = Simulink.SimulationInput(name_hrm);
    simRoot = setVariable(simRoot, 'A', cfg.A, 'Workspace', name_hrm);
    simRoot = setVariable(simRoot, 'B', cfg.B, 'Workspace', name_hrm);
    simRoot = setVariable(simRoot, 'mdot_c_drift', cfg.mdot_c_drift, 'Workspace', name_hrm);
    simRoot = setVariable(simRoot, 'vel_ACS', cfg.vel_ACS, 'Workspace', name_hrm);
    simRoot = setModelParameter(simRoot, 'StopTime', string(t_q));

    tt = transpose(0:t_sam:t_q);
    n  = numel(tt);
    control_series = [tt, zeros(n,1), zeros(n,1)];       % no control before t_q
    simRoot = simRoot.setExternalInput(control_series);
    simRoot = setModelParameter(simRoot, 'SaveFinalState', 'on', ...
        'SaveOperatingPoint', 'on', 'FinalStateName', 'finalOPRoot');

    ttic = tic;
    outRoot = sim(simRoot);
    G.hrm_seconds = toc(ttic);
    finalOP_Raiz = outRoot.finalOPRoot;
    fprintf('   HRM to t=%.1f s done in %.1f s\n', t_q, G.hrm_seconds);

    tw = outRoot.tw_max_signal.Data(:);
    G.hrm_peak_tw = max(tw);
    G.hrm_final_tw = tw(end);
    fprintf('   HRM tw_max: final %.4f K, peak %.4f K\n', G.hrm_final_tw, G.hrm_peak_tw);

    % second HRM run to the ROM anticipation instant
    simRoot2 = simRoot;
    simRoot2 = setModelParameter(simRoot2, 'StopTime', string(t_q - delayed));
    simRoot2 = setModelParameter(simRoot2, 'SaveState', 'off', ...
        'SaveFinalState', 'on', 'SaveOperatingPoint', 'on', ...
        'FinalStateName', 'finalOPRoot2');
    outRoot2 = sim(simRoot2);
    finalOP2_Raiz = outRoot2.finalOPRoot2;
    close_system(name_hrm, 0);
    fprintf('   HRM to t=%.1f s done (ROM init state)\n\n', t_q - delayed);

    %% ---- Stage 2: BiLSTM MC-dropout inference ----
    fprintf('-- Stage 2: BiLSTM MC-dropout inference (%d trajectories) --\n', cfg.N_dropout);
    len_seq = 200;
    X_1 = outRoot.mdot_c_leak_signal.Data(:);
    X_2 = outRoot.Te_c_signal.Data(12,:,:);
    X_2 = X_2(:);
    if numel(X_1) < len_seq || numel(X_2) < len_seq
        error('generate:shortSignal', ...
            'HRM signals shorter than len_seq (%d vs %d/%d)', len_seq, numel(X_1), numel(X_2));
    end
    X_input = zeros(len_seq, 2);
    X_input(:,1) = X_1(end-len_seq+1:end);
    X_input(:,2) = X_2(end-len_seq+1:end);
    X_norm = zeros(size(X_input));
    for i = 1:2
        X_norm(:,i) = (X_input(:,i) - X_minmax(i,1)) / (X_minmax(i,2) - X_minmax(i,1));
    end

    py_module = py.importlib.import_module('mc_inference');
    ttic = tic;
    res = py_module.mc_dropout_trajectories(py_module.model, py.numpy.array(X_norm), ...
        py.int(cfg.N_dropout));
    G.inference_seconds = toc(ttic);
    Y_norm = double(py.numpy.array(res));
    Y_denorm = zeros(size(Y_norm));
    for i = 1:size(Y_norm,3)
        Y_denorm(:,:,i) = Y_norm(:,:,i)*(Y_minmax(i,2)-Y_minmax(i,1)) + Y_minmax(i,1);
    end
    fprintf('   inference done in %.2f s, output %s\n', G.inference_seconds, mat2str(size(Y_denorm)));

    [~, len2, nfeat] = size(Y_denorm);
    Y50 = zeros(len2, nfeat); Y95 = zeros(len2, nfeat);
    for f = 1:nfeat
        fd = squeeze(Y_denorm(:,:,f));
        Y50(:,f) = prctile(fd, 50, 1)';
        Y95(:,f) = prctile(fd, 95, 1)';
    end
    G.Y50_pred = Y50(1:cfg.N_pred, :);
    G.Y95_pred = Y95(1:cfg.N_pred, :);
    fprintf('   predicted delta_tw_max (median) over horizon: min %.4f, max %.4f K\n\n', ...
        min(G.Y50_pred(:,2)), max(G.Y50_pred(:,2)));

    %% ---- Stage 3: ROM initial operating point ----
    fprintf('-- Stage 3: ROM operating point --\n');
    load_system([name_rom '.slx']);
    simrom = Simulink.SimulationInput(name_rom);
    simrom = setModelParameter(simrom, 'FastRestart', 'off');
    simrom = setModelParameter(simrom, 'SaveFinalState', 'on', ...
        'SaveOperatingPoint', 'on', 'FinalStateName', 'finalOP');
    simrom = setModelParameter(simrom, 'StopTime', string(t_q - delayed));
    inport = [tt, outRoot.mdot_c_leak_signal.Data(:), zeros(n,1), zeros(n,1), zeros(n,1)];
    simrom = simrom.setExternalInput(inport);
    simrom = setVariable(simrom, 'vel_ACS', cfg.vel_ACS, 'Workspace', name_rom);
    ttic = tic;
    outrom = sim(simrom);
    G.rom_seconds = toc(ttic);
    romOP = outrom.finalOP;
    fprintf('   ROM to t=%.1f s done in %.1f s\n\n', t_q - delayed, G.rom_seconds);

    %% ---- Stage 4: assemble det_params ----
    fprintf('-- Stage 4: assembling det_params --\n');
    dp = struct();
    dp.N_par = 1;  dp.M_workers = 8;
    dp.finalOP_Raiz = finalOP_Raiz;
    dp.finalOP2_Raiz = finalOP2_Raiz;
    dp.romOP = romOP;
    dp.tsim = tsim;
    dp.time_snapshot_root = t_q;
    dp.t_sam = t_sam;
    dp.delayed_time_rom = delayed;
    dp.outRoot = outRoot;
    dp.outrom = outrom;
    dp.Y_denorm = reshape(Y_denorm(1,1:cfg.N_pred,:), [cfg.N_pred 2]);
    dp.Y_denorm_50 = G.Y50_pred;
    dp.Y_denorm_95 = G.Y95_pred;
    dp.sampled_vel_ACS = cfg.vel_ACS;
    dp.Lambda_t = Lambda_t;  dp.Mu_t = Mu_t;
    dp.Ncomps = 5;
    dp.time_snapshots = t_q;
    dp.snapshot_id = 1;
    dp.N_resampling = 2;  dp.N_branches = 2;
    dp.name_root = name_hrm;  dp.name_Wks = name_wks;
    dp.control_series = control_series;
    G.det_params = dp;

    outfile = sprintf('severe_det_params_%s_A%.2f_B%.0f_vel%.3f.mat', ...
        cfg.tag, cfg.A, cfg.B, cfg.vel_ACS);
    det_params = dp;  %#ok<NASGU>
    save(outfile, 'det_params', 'G', '-v7.3');
    fprintf('   saved %s\n', outfile);
    G.outfile = outfile;

    %% ---- Stage 5: quick severity check ----
    fprintf('\n-- Stage 5: severity check (DET at u_min and u_max) --\n');
    target = CalculateTs(100, 543);
    U_Y = 1244; span = U_Y - target;
    d_mpc = ones(1, cfg.N_pred + 1);
    for uv = [56.69 56.75]
        u_trial = uv * ones(1, cfg.N_pred + 1);
        try
            [fd, fr] = DET_2(dp, u_trial, d_mpc);
            opts = struct('dsm_thr', 1.0, 'lambda_pct', 0, 'U_y', U_Y);
            [Tp, ~, ~, tm] = fitness_cost_dsm(fd, fr, cfg.N_pred, target, ...
                t_q, delayed, t_sam, u_trial, opts);
            dsm = (tm.peak_Tp - target)/span;
            fprintf('   u=%.5f -> peak Tp %.4f K, DSM %.6f, margin %.4f K, J_error %.4e\n', ...
                uv, tm.peak_Tp, dsm, U_Y - tm.peak_Tp, tm.J_error);
            if uv == 56.75
                G.dsm_at_umax = dsm; G.peak_Tp_umax = tm.peak_Tp;
            else
                G.dsm_at_umin = dsm; G.peak_Tp_umin = tm.peak_Tp;
            end
        catch ME
            fprintf('   u=%.5f FAILED: %s\n', uv, ME.message);
        end
    end

    fprintf('\n===== VERDICT =====\n');
    if isfield(G, 'dsm_at_umax')
        fprintf('DSM at full coolant (u_max): %.6f  (%.2f%% of the margin)\n', ...
            G.dsm_at_umax, 100*G.dsm_at_umax);
        if G.dsm_at_umax > 0.05
            G.conclusion = sprintf(['SEVERE CASE ACHIEVED: even at maximum coolant the PCT ' ...
                'consumes %.2f%% of the safety margin, so tracking alone cannot resolve the ' ...
                'excursion and the DSM penalty is no longer redundant. Run the joint ' ...
                '(dsm_thr, lambda) sweep on %s.'], 100*G.dsm_at_umax, G.outfile);
        else
            G.conclusion = sprintf(['Not yet severe: DSM at u_max is %.6f (below 5%%). ' ...
                'A and B are already at their upper bounds, so consider a later query ' ...
                'instant (closer to the thermal peak) or a longer horizon.'], G.dsm_at_umax);
        end
        fprintf('\n%s\n', G.conclusion);
    end
    fprintf('=============================================\n\n');
end
