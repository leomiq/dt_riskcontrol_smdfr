function C = test_control_authority(matfile)
% TEST_CONTROL_AUTHORITY  How much PCT can the actuator actually move?
%
% This is the precondition for the whole J_psm / DSM demonstration. If the full
% admissible current range cannot shift the predicted PCT by more than the
% distance to a defensible DSM threshold, then no (DSM_thr, lambda) pair can
% make the safety penalty steer the decision -- the optimisation problem would
% be ill-posed for that purpose, regardless of how the penalty is written.
%
% Evaluates the DET at the fixed query state for the extremes and the midpoint
% of u in [u_min, u_max], with NO optimizer (3 evaluations instead of 420), and
% reports the resulting PCT / DSM spread plus the mass flow each current gives.

    if nargin < 1 || isempty(matfile)
        matfile = 'LJ_test_sampled_vel_ACS_2_loop_1180.mat';
    end

    C = struct();
    C.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fprintf('\n===== Control-authority test (fixed state) =====\n');
    fprintf('Started %s\n\n', C.timestamp);

    S  = load(matfile);
    fprintf('Regenerating operating points...\n');
    dp = refresh_operating_points(S);
    if ~isfield(dp, 'Y_denorm_50'), dp.Y_denorm_50 = dp.Y_denorm; end

    t_root  = double(dp.time_snapshot_root);
    t_sam   = double(dp.t_sam);
    delayed = double(dp.delayed_time_rom);
    N_pred  = size(dp.Y_denorm_50, 1);
    target  = CalculateTs(100, 543);
    U_Y     = 1244;

    u_min = 56.69; u_max = 56.75;
    C.u_min = u_min; C.u_max = u_max; C.u_width = u_max - u_min;
    C.target = target; C.U_Y = U_Y; C.dsm_span = U_Y - target;
    fprintf('u range      : [%.5f, %.5f] A   (width %.5f A)\n', u_min, u_max, C.u_width);
    fprintf('target       : %.4f K\n', target);
    fprintf('DSM span     : %.4f K  (U_y - target)\n\n', C.dsm_span);

    % Mass flow the analytical model gives for each extreme -- the physical
    % coolant injection the current commands.
    try
        m_lo = CurrentToMassFlow(u_min * ones(1, N_pred+1));
        m_hi = CurrentToMassFlow(u_max * ones(1, N_pred+1));
        C.mdot_min = m_lo(1); C.mdot_max = m_hi(1);
        C.mdot_span = C.mdot_max - C.mdot_min;
        fprintf('coolant mass flow: u_min -> %.6f, u_max -> %.6f  (span %.6f)\n\n', ...
            C.mdot_min, C.mdot_max, C.mdot_span);
    catch ME
        fprintf('(CurrentToMassFlow unavailable: %s)\n\n', ME.message);
    end

    time_delayed = t_root - delayed : t_sam : t_root - t_sam;
    d_mpc = ones(1, N_pred + 1);

    u_cases = [u_min, 0.5*(u_min+u_max), u_max];
    labels  = {'u_min', 'u_mid', 'u_max'};

    C.cases = struct('label', {}, 'u', {}, 'peak_Tp', {}, 'min_Tp', {}, ...
                     'mean_Tp', {}, 'peak_DSM', {}, 'J', {}, 'n_samples', {}, ...
                     'n_nodes', {}, 'sim_s', {}, 'note', {});

    for k = 1:numel(u_cases)
        uv = u_cases(k);
        fprintf('---- %s: u = %.5f A ----\n', labels{k}, uv);
        rec = struct('label', labels{k}, 'u', uv, 'peak_Tp', NaN, 'min_Tp', NaN, ...
                     'mean_Tp', NaN, 'peak_DSM', NaN, 'J', NaN, 'n_samples', NaN, ...
                     'n_nodes', NaN, 'sim_s', NaN, 'note', '');
        u_trial = uv * ones(1, N_pred + 1);
        tt = tic;
        try
            [fd, fr] = DET_2(dp, u_trial, d_mpc);
            rec.n_nodes = numel(fd{1}.tree);
            [Tp, J, ~] = fitness_cost(fd, fr, N_pred, target, t_root, delayed, t_sam, u_trial);
            rec.sim_s = toc(tt);
            rec.n_samples = numel(Tp);
            rec.peak_Tp = max(Tp);
            rec.min_Tp = min(Tp);
            rec.mean_Tp = mean(Tp);
            rec.peak_DSM = (rec.peak_Tp - target) / C.dsm_span;
            rec.J = J;
            fprintf('   nodes=%d  samples=%d  peak Tp=%.4f K  DSM=%.6f  J=%.6e  (%.1f s)\n', ...
                rec.n_nodes, rec.n_samples, rec.peak_Tp, rec.peak_DSM, J, rec.sim_s);
        catch ME
            rec.sim_s = toc(tt);
            rec.note = sprintf('[%s] %s', ME.identifier, ME.message);
            fprintf('   FAILED: %s\n', rec.note);
        end
        C.cases(end+1) = rec;
        fprintf('\n');
    end

    %% ---- Authority verdict ----
    fprintf('===== CONTROL AUTHORITY =====\n');
    peaks = [C.cases.peak_Tp];
    dsms  = [C.cases.peak_DSM];
    Js    = [C.cases.J];
    ok = ~isnan(peaks);
    if sum(ok) < 2
        fprintf('Not enough successful evaluations to judge.\n');
        C.conclusion = 'inconclusive';
        finish(C); return;
    end

    C.pct_authority_K = max(peaks(ok)) - min(peaks(ok));
    C.dsm_authority   = max(dsms(ok)) - min(dsms(ok));
    C.J_spread        = max(Js(ok)) - min(Js(ok));
    fprintf('PCT swing over the FULL current range : %.6f K\n', C.pct_authority_K);
    fprintf('DSM swing over the FULL current range : %.6f\n', C.dsm_authority);
    fprintf('J spread                              : %.6e\n', C.J_spread);
    for k = 1:numel(C.cases)
        if ~isnan(C.cases(k).peak_Tp)
            fprintf('   %-6s u=%.5f -> peak Tp %.4f K (DSM %.6f), J %.6e\n', ...
                C.cases(k).label, C.cases(k).u, C.cases(k).peak_Tp, ...
                C.cases(k).peak_DSM, C.cases(k).J);
        end
    end

    fprintf('\n-- Can a DSM threshold be steered by the actuator? --\n');
    % A threshold is "reachable" if the actuator can move the peak DSM across it.
    dsm_lo = min(dsms(ok)); dsm_hi = max(dsms(ok));
    fprintf('actuator can place peak DSM anywhere in [%.6f, %.6f]\n', dsm_lo, dsm_hi);
    for thr = [0.20 0.10 0.05 0.02 0.01 0.005 0.002]
        if thr > dsm_hi
            verdict = 'never violated -> penalty always 0, lambda irrelevant';
        elseif thr < dsm_lo
            verdict = 'always violated -> penalty is a constant offset, no steering';
        else
            verdict = '>>> STEERABLE: actuator can cross this threshold <<<';
        end
        fprintf('   DSM_thr %.3f : %s\n', thr, verdict);
    end

    C.steerable = any(arrayfun(@(x) x >= dsm_lo && x <= dsm_hi, ...
        [0.20 0.10 0.05 0.02 0.01 0.005 0.002]));
    if C.steerable
        C.conclusion = ['At least one DSM threshold lies inside the actuator''s ' ...
            'reachable DSM band, so a (DSM_thr, lambda) pair can steer the decision.'];
    else
        C.conclusion = sprintf(['No tested DSM threshold lies inside the reachable ' ...
            'band [%.6f, %.6f]. The actuator moves the peak DSM by only %.6f, so the ' ...
            'safety penalty cannot reorient the decision at this operating point.'], ...
            dsm_lo, dsm_hi, C.dsm_authority);
    end
    fprintf('\n%s\n', C.conclusion);
    finish(C);
end

function finish(C)
    outfile = sprintf('test_control_authority_%s.mat', datestr(now,'yyyymmdd_HHMMSS'));
    save(outfile, 'C');
    fprintf('\nSaved %s\n=================================\n\n', outfile);
end
