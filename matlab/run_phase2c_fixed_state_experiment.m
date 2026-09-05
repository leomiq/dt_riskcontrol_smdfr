%% run_phase2c_fixed_state_experiment.m
% Fixed-state repeated-seed experiment for Chapter 7 MPC/optimisation
% (closes P-M4, PAT-C-065 and PAT-C-066 once executed and reported).
%
% REQUIRED FILES (copy into this folder or add to the MATLAB path):
%   LJ_test_sampled_vel_ACS_2_loop_1180.mat   archived LJ workspace, "LJ test 11_10"
%   luus_jaakola_det.m, ga_det.m, DET_2.m, fitness_cost.m,
%   CurrentToMassFlow.m, rom_MCDET_v2_2022b_control_CMEX.slx (+ .slxc),
%   all *.mexw64 from the archived CINECA folder, max_min_valuess.mat
%
% DESIGN
%   - One fixed query state: t_q = 1180 s. The archived det_params, target,
%     initial guess, bounds and stopping rule are reused as-is, so the
%     prediction model, bounds, initial state and stopping rule are FIXED.
%   - LJ: 20 independent recorded seeds with the archived configuration
%     (max_passes=30, M=5, R=15, r_init=0.03 A, gamma=0.95, early stop after
%     5 consecutive passes below 1) and therefore the archived evaluation budget.
%   - GA: 20 independent recorded seeds with the archived operators
%     (roulette-wheel selection, crossover fraction 0.80, mutation 0.30,
%     one elite) but pop_size = 20 so that the evaluation budget is
%     20 + 20*20 = 420, i.e. inside the LJ range 376-451 recorded for the
%     archived primary series (equivalent objective-evaluation budget).
%   - Every seed is recorded; raw per-run arrays are saved to MAT.
%
% PER-RUN METRICS (fixed query state)
%   f_opt       optimizer objective value
%   u_apply     first applied ACS current of the optimized sequence
%   minDSM_K    minimum DSM of the optimal predicted trajectory
%               (1244 K safety limit minus the predicted peak Tp)
%   runtime_s   optimizer-loop time = sum of nonzero t_vector entries
%   passes / generations and evaluation counts
%
% NOTES
%   - The archived 11_10 series predates the Y_denorm_50 field of DET_2.m;
%     its DET perturbation signal was stored in det_params.Y_denorm (see
%     DET_2.asv, lines 124/408). The harness maps that signal into the field
%     read by the current DET_2.m so the archived inputs are reproduced.
%   - The MAT files do not encode hardware. Record CPU, RAM, MATLAB release
%     and worker count in the run log for the provenance statement.
%   - Closed-loop trajectory metrics (overshoot, settling time, steady-state
%     error, control total variation) require full 61-query closed-loop
%     realizations; a stage-2 skeleton is provided at the end of this file.

function run_phase2c_fixed_state_experiment()
    tStart = tic;
    outName = 'phase2c_fixed_state_results.mat';

    %% ---- Load the archived fixed query state (t_q = 1180 s) ----
    wsFile = 'LJ_test_sampled_vel_ACS_2_loop_1180.mat';
    if ~isfile(wsFile)
        error('Archived workspace not found in current folder: %s', wsFile);
    end
    S = load(wsFile);
    det_params = refresh_operating_points(S);
    target     = S.target;
    d_mpc      = S.d_mpc;
    u0         = S.u0;
    u_bounds   = S.u_bounds;
    t0         = S.time_snapshot_root;
    t_delayed  = S.delayed_time_rom;
    ts         = S.t_sam;
    popEquiv   = 20;                 % 20 + 20*20 = 420 evals (LJ archived 376-451)
    gens       = S.num_generations;  % 20
    nSeeds     = 20;
    seeds      = 1:nSeeds;
    config     = struct('query_time_s', t0, 'seeds', seeds, 'nSeeds', nSeeds, ...
        'LJ', struct('max_passes', S.max_passes, 'M', S.M, 'R', S.R, ...
                     'r_init', S.r_init, 'gamma', S.gamma), ...
        'GA', struct('generations', gens, 'pop_size', popEquiv, ...
                     'crossover_frac', S.crossover_frac, 'mutation_rate', S.mutation_rate), ...
        'u_min', S.u_min, 'u_max', S.u_max, 'target', target, ...
        'Y_denorm_50_mapping', 'archived single-trial det_params.Y_denorm (11_10 series)', ...
        'hardware_notes', 'record CPU, RAM, MATLAB release and worker count in the run log', ...
        'run_date', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

    % Map the archived single-trial perturbation signal into the field read by
    % the current DET_2.m (the 11_10 series stored it as det_params.Y_denorm).
    if ~isfield(det_params, 'Y_denorm_50')
        det_params.Y_denorm_50 = det_params.Y_denorm;
    end

    %% ---- Parallel pool (reused by parfor in both optimizers) ----
    nWorkers = min(feature('numcores'), 8);
    if isempty(gcp('nocreate'))
        parpool('local', nWorkers);
    end


    %% ---- Stage 1: LJ repeated runs (archived configuration) ----
    fprintf('=== LJ: %d fixed-state runs (recorded seeds 1..%d) ===\n', nSeeds, nSeeds);
    LJ = repmat(struct(), 1, nSeeds);
    for k = 1:nSeeds
        rng(seeds(k), 'twister');
        [Tp, u_opt, f_opt, hist, tvec] = luus_jaakola_det( ...
            u0, u_bounds, S.max_passes, S.M, S.r_init, S.gamma, S.R, ...
            det_params, target, d_mpc, t0, t_delayed, ts);
        dsm = 1244 - max(Tp(:,1));
        LJ(k).seed        = seeds(k);
        LJ(k).f_opt       = f_opt;
        LJ(k).u_apply     = u_opt(1);
        LJ(k).minDSM_K    = dsm;
        LJ(k).runtime_s   = sum(tvec);
        LJ(k).passes      = nnz(tvec);
        LJ(k).evaluations = 1 + nnz(tvec) * S.M * S.R;
        LJ(k).Tp_opt      = Tp;
        LJ(k).u_opt       = u_opt;
        LJ(k).history     = hist;
        LJ(k).t_vector    = tvec;
        fprintf('LJ seed %2d: f=%.6f  u1=%.6f  minDSM=%.3f K  t=%.2f s  passes=%d\n', ...
            seeds(k), f_opt, u_opt(1), dsm, sum(tvec), nnz(tvec));
    end

    % Intermediate save so the LJ stage survives any later failure
    save('phase2c_fixed_state_results_LJ_only.mat', 'LJ', 'seeds', 'config', '-v7.3');

    %% ---- Stage 2: GA repeated runs (equivalent evaluation budget) ----
    fprintf('=== GA: %d fixed-state runs, pop_size=%d, %d generations ===\n', ...
        nSeeds, popEquiv, gens);
    GA = repmat(struct(), 1, nSeeds);
    for k = 1:nSeeds
        rng(seeds(k), 'twister');
        [Tp, u_opt, f_opt, hist, tvec] = ga_det( ...
            u0, u_bounds, gens, popEquiv, S.crossover_frac, S.mutation_rate, ...
            det_params, target, d_mpc, t0, t_delayed, ts);
        dsm = 1244 - max(Tp(:,1));
        GA(k).seed        = seeds(k);
        GA(k).f_opt       = f_opt;
        GA(k).u_apply     = u_opt(1);
        GA(k).minDSM_K    = dsm;
        GA(k).runtime_s   = sum(tvec);
        GA(k).generations = nnz(tvec);
        GA(k).evaluations = popEquiv + nnz(tvec) * popEquiv;
        GA(k).Tp_opt      = Tp;
        GA(k).u_opt       = u_opt;
        GA(k).history     = hist;
        GA(k).t_vector    = tvec;
        fprintf('GA  seed %2d: f=%.6f  u1=%.6f  minDSM=%.3f K  t=%.2f s  gen=%d\n', ...
            seeds(k), f_opt, u_opt(1), dsm, sum(tvec), nnz(tvec));
    end

    %% ---- Summary statistics (mean, SD, median, IQR, range) ----
    stat = @(v) struct('mean', mean(v), 'sd', std(v), 'median', median(v), ...
                       'iqr', iqr(v), 'min', min(v), 'max', max(v));
    LJstats = struct('f_opt', stat([LJ.f_opt]), 'u_apply', stat([LJ.u_apply]), ...
                     'minDSM_K', stat([LJ.minDSM_K]), 'runtime_s', stat([LJ.runtime_s]), ...
                     'passes', stat([LJ.passes]));
    GAstats = struct('f_opt', stat([GA.f_opt]), 'u_apply', stat([GA.u_apply]), ...
                     'minDSM_K', stat([GA.minDSM_K]), 'runtime_s', stat([GA.runtime_s]), ...
                     'generations', stat([GA.generations]));

    fprintf('\n=== LJ summary (n=%d) ===\n', nSeeds);
    fprintf('f_opt : mean %.4f  sd %.4f  med %.4f  IQR %.4f  range [%.4f, %.4f]\n', ...
        LJstats.f_opt.mean, LJstats.f_opt.sd, LJstats.f_opt.median, ...
        LJstats.f_opt.iqr, LJstats.f_opt.min, LJstats.f_opt.max);
    fprintf('u_apply: mean %.6f  sd %.6f  med %.6f  IQR %.6f  range [%.6f, %.6f] A\n', ...
        LJstats.u_apply.mean, LJstats.u_apply.sd, LJstats.u_apply.median, ...
        LJstats.u_apply.iqr, LJstats.u_apply.min, LJstats.u_apply.max);
    fprintf('minDSM: mean %.3f  sd %.3f  med %.3f  IQR %.3f  range [%.3f, %.3f] K\n', ...
        LJstats.minDSM_K.mean, LJstats.minDSM_K.sd, LJstats.minDSM_K.median, ...
        LJstats.minDSM_K.iqr, LJstats.minDSM_K.min, LJstats.minDSM_K.max);
    fprintf('runtime: mean %.2f  sd %.2f  med %.2f  IQR %.2f  range [%.2f, %.2f] s\n', ...
        LJstats.runtime_s.mean, LJstats.runtime_s.sd, LJstats.runtime_s.median, ...
        LJstats.runtime_s.iqr, LJstats.runtime_s.min, LJstats.runtime_s.max);

    fprintf('\n=== GA summary (n=%d) ===\n', nSeeds);
    fprintf('f_opt : mean %.4f  sd %.4f  med %.4f  IQR %.4f  range [%.4f, %.4f]\n', ...
        GAstats.f_opt.mean, GAstats.f_opt.sd, GAstats.f_opt.median, ...
        GAstats.f_opt.iqr, GAstats.f_opt.min, GAstats.f_opt.max);
    fprintf('u_apply: mean %.6f  sd %.6f  med %.6f  IQR %.6f  range [%.6f, %.6f] A\n', ...
        GAstats.u_apply.mean, GAstats.u_apply.sd, GAstats.u_apply.median, ...
        GAstats.u_apply.iqr, GAstats.u_apply.min, GAstats.u_apply.max);
    fprintf('minDSM: mean %.3f  sd %.3f  med %.3f  IQR %.3f  range [%.3f, %.3f] K\n', ...
        GAstats.minDSM_K.mean, GAstats.minDSM_K.sd, GAstats.minDSM_K.median, ...
        GAstats.minDSM_K.iqr, GAstats.minDSM_K.min, GAstats.minDSM_K.max);
    fprintf('runtime: mean %.2f  sd %.2f  med %.2f  IQR %.2f  range [%.2f, %.2f] s\n', ...
        GAstats.runtime_s.mean, GAstats.runtime_s.sd, GAstats.runtime_s.median, ...
        GAstats.runtime_s.iqr, GAstats.runtime_s.min, GAstats.runtime_s.max);

    %% ---- Save raw arrays, seeds and summary ----
    % config is defined before the LJ stage (above)
    save(outName, 'LJ', 'GA', 'LJstats', 'GAstats', 'config', 'seeds');
    fprintf('\nSaved %s  (elapsed %.1f min)\n', outName, toc(tStart)/60);
end

% --------------------------------------------------------------------------
% STAGE 2 SKELETON (optional, much longer): full closed-loop realizations.
% Trajectory-level metrics (overshoot, settling time, steady-state error,
% control total variation) require the 61-query closed loop with the same
% optimizer invoked at every query time. That loop is the body of
% DET_scheduler_v1_perc_online_LJ.m / DET_scheduler_v1_perc_online_GA.m;
% wrap it as a function, set rng(seed,'twister') before each realization,
% keep the accident seed fixed (rng(1,"twister") before scenario init), and
% run calculate_performance_metrics on the resulting tw_max_V2P_trajectory
% and control_series, as the archived Performance_analysis.txt does. At ~50 s
% per LJ query, 20 realizations of 61 queries take of the order of 17 hours.
% --------------------------------------------------------------------------
