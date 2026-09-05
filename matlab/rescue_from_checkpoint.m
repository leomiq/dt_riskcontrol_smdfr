function R = rescue_from_checkpoint(ckptFile, outFile)
% RESCUE_FROM_CHECKPOINT  Rebuild the final result MAT from a checkpoint whose
% run completed all queries but died in the KPI step.
%
% The widened-bound penalised realization finished all 61 queries and then threw
% MATLAB:badsubscript inside calculate_performance_metrics (settling-time index
% out of range when the trajectory is still outside the tolerance band at the
% last sample -- now fixed there). The trajectory data survived in the
% checkpoint, so there is no need to repeat ~5 h of simulation.

    if nargin < 1 || isempty(ckptFile)
        ckptFile = 'widepen_umax56p90_LJ_seed1_checkpoint.mat';
    end
    if nargin < 2 || isempty(outFile)
        outFile = 'widepen_umax56p90_LJ_seed1.mat';
    end

    fprintf('\n===== Rescue from checkpoint =====\n');
    fprintf('checkpoint: %s\n', ckptFile);

    C = load(ckptFile);
    time_snapshots        = C.time_snapshots;
    u_V2P_trajectory      = C.u_V2P_trajectory;
    first_actions         = C.first_actions;
    optimizer_runtimes    = C.optimizer_runtimes;
    control_series        = C.control_series;
    tw_max_V2P_trajectory = C.tw_max_V2P_trajectory;
    solver                = C.solver;
    seed                  = C.seed;
    n_done                = C.snapshot_id_done;

    fprintf('queries completed: %d of %d\n', n_done, numel(time_snapshots));
    if n_done < numel(time_snapshots)
        fprintf('!! incomplete -- KPIs would cover only the finished part\n');
    end

    % Same extraction as run_phase2c_full_loop's KPI block
    [ctrl_idx, ~] = find(control_series(:,1) == time_snapshots);
    [v2p_idx, ~]  = find(tw_max_V2P_trajectory(:,1) == time_snapshots);
    PCT     = tw_max_V2P_trajectory(v2p_idx, 2);
    control = control_series(ctrl_idx, 2);
    target  = CalculateTs(100, 543);

    fprintf('PCT samples: %d, control samples: %d\n', numel(PCT), numel(control));
    fprintf('PCT: first %.4f, last %.4f, min %.4f, max %.4f K\n', ...
        PCT(1), PCT(end), min(PCT), max(PCT));

    params = struct('PCT_ref', target, 'PCT_max', 1244, 'tolerance', 0.3);
    metrics = calculate_performance_metrics(time_snapshots, PCT, control, params);
    metrics.optimizer_runtime_median_s = median(optimizer_runtimes);
    metrics.optimizer_runtime_mean_s   = mean(optimizer_runtimes);
    metrics.optimizer_runtime_max_s    = max(optimizer_runtimes);

    fprintf('\n-- KPIs --\n');
    fn = fieldnames(metrics);
    for i = 1:numel(fn)
        v = metrics.(fn{i});
        if isscalar(v) && isnumeric(v)
            fprintf('  %-32s = %.6f\n', fn{i}, v);
        elseif islogical(v)
            fprintf('  %-32s = %d\n', fn{i}, v);
        end
    end

    save(outFile, 'metrics', 'solver', 'seed', 'time_snapshots', 'PCT', 'control', ...
        'u_V2P_trajectory', 'first_actions', 'optimizer_runtimes', ...
        'tw_max_V2P_trajectory', 'control_series');
    fprintf('\nsaved %s\n', outFile);
    R = struct('metrics', metrics, 'PCT', PCT, 'control', control, ...
               'time_snapshots', time_snapshots, 'outFile', outFile);
    fprintf('==================================\n\n');
end
