function V = verify_lj_suboptimality(tq_list, n_grid)
% VERIFY_LJ_SUBOPTIMALITY  Is the late-horizon PCT drift caused by LJ returning
% sub-optimal solutions, rather than by the cost function?
%
% EVIDENCE THAT MOTIVATES THIS (diagnose_drift_cause.py, on the widened
% penalised realization):
%   * the drift happens entirely AFTER the penalty switches off at t=1169.5
%     (J_psm = 0 for the remaining 41 queries), so lambda cannot be driving it;
%   * the drift ACCELERATES (-104 -> -152 mK/s), so a longer horizon would make
%     it worse, not better;
%   * in 13 of 40 queries where the PCT is already BELOW setpoint with no
%     penalty, LJ commanded MORE coolant than the archived u_max -- inconsistent
%     with minimising J_error;
%   * LJ finished in about half the baseline's runtime (185 s vs 346 s), i.e. it
%     is stopping early on a search space that got 3.5x wider (0.06 -> 0.21 A)
%     without any retuning: r_init = 0.5*range grew 0.03 -> 0.105 A, and with
%     gamma = 0.95 over 30 passes the final radius is still ~0.023 A, about the
%     width of the whole original range.
%
% METHOD. For a few query instants in the drift region, evaluate the true cost on
% a fine grid of constant-u candidates (the DET is deterministic given the state,
% so this is a clean reference), then compare:
%     J(u_LJ)      what LJ actually paid
%     min_u J(u)   what was available on the grid
% A large gap, together with u_LJ far from argmin, confirms sub-optimality.
%
% This does NOT re-run the closed loop: each query is re-evaluated from the
% archived operating point at that instant, so the states are the archived ones
% rather than the penalised run's own. It therefore answers "can LJ find the
% optimum of a query like this one", which is exactly the optimizer question.
%
%   tq_list  query instants to probe (default [1176 1180 1184 1188])
%   n_grid   grid points across [u_min, u_max] (default 22)

    if nargin < 1 || isempty(tq_list), tq_list = [1176 1180 1184 1188]; end
    if nargin < 2 || isempty(n_grid),  n_grid = 22; end

    V = struct();
    V.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fprintf('\n===== LJ sub-optimality check (drift region) =====\n');
    fprintf('Started %s\n', V.timestamp);
    fprintf('query instants: %s\n', mat2str(tq_list));
    fprintf('grid: %d points over [56.69, 56.90]\n\n', n_grid);

    try, suppress_simulink_prompts(); catch, end

    % Penalty settings of the widened run
    dsm_thr = 0.00237737;
    lambda  = 1e6;
    U_Y = 1244;
    u_min = 56.69; u_max = 56.90;
    u_grid = linspace(u_min, u_max, n_grid);
    target = CalculateTs(100, 543);
    span = U_Y - target;

    % What LJ applied at those instants in the widened run
    runfile = 'widepen_umax56p90_LJ_seed1.mat';
    if ~isfile(runfile)
        error('verify:noRun', 'realization %s not found in %s', runfile, pwd);
    end
    W = load(runfile);
    t_all = W.time_snapshots(:)';
    u_all = W.u_V2P_trajectory(:)';

    V.tq = tq_list;
    V.u_lj = nan(size(tq_list));
    V.u_best = nan(size(tq_list));
    V.J_lj = nan(size(tq_list));
    V.J_best = nan(size(tq_list));
    V.gap_abs = nan(size(tq_list));
    V.gap_rel = nan(size(tq_list));
    V.grids = cell(size(tq_list));

    for q = 1:numel(tq_list)
        tq = tq_list(q);
        fprintf('---- t_q = %.1f s ----\n', tq);

        k = find(abs(t_all - tq) < 1e-9, 1);
        if isempty(k)
            fprintf('   not in the realization; skipping\n\n'); continue;
        end
        u_lj = u_all(k);
        V.u_lj(q) = u_lj;
        fprintf('   LJ applied u = %.5f A\n', u_lj);

        wsfile = sprintf('LJ_test_sampled_vel_ACS_2_loop_%g.mat', tq);
        if ~isfile(wsfile)
            fprintf('   archived workspace %s absent; skipping\n\n', wsfile); continue;
        end
        S = load(wsfile);
        dp = refresh_operating_points(S);
        if ~isfield(dp, 'Y_denorm_50'), dp.Y_denorm_50 = dp.Y_denorm; end
        dp.dsm_thr = dsm_thr; dp.lambda_dsm = lambda;

        t_sam = double(dp.t_sam);
        delayed = double(dp.delayed_time_rom);
        N_pred = size(dp.Y_denorm_50, 1);
        d_mpc = ones(1, N_pred + 1);

        Jg = nan(1, n_grid); Peak = nan(1, n_grid); Dsm = nan(1, n_grid);
        Jerr = nan(1, n_grid); Jpsm = nan(1, n_grid);
        for i = 1:n_grid
            u_trial = u_grid(i) * ones(1, N_pred + 1);
            try
                [fd, fr] = DET_2(dp, u_trial, d_mpc);
                opts = struct('dsm_thr', dsm_thr, 'lambda_pct', lambda, 'U_y', U_Y);
                [Tp, Jv, ~, tm] = fitness_cost_dsm(fd, fr, N_pred, target, ...
                    tq, delayed, t_sam, u_trial, opts);
                Jg(i) = Jv; Peak(i) = tm.peak_Tp;
                Dsm(i) = (tm.peak_Tp - target)/span;
                Jerr(i) = tm.J_error; Jpsm(i) = tm.J_psm;
            catch ME
                fprintf('   u=%.5f failed: %s\n', u_grid(i), ME.message);
            end
        end

        % cost at LJ's own action, evaluated the same way
        J_at_lj = NaN;
        try
            u_trial = u_lj * ones(1, N_pred + 1);
            [fd, fr] = DET_2(dp, u_trial, d_mpc);
            opts = struct('dsm_thr', dsm_thr, 'lambda_pct', lambda, 'U_y', U_Y);
            [~, J_at_lj] = fitness_cost_dsm(fd, fr, N_pred, target, ...
                tq, delayed, t_sam, u_trial, opts);
        catch ME
            fprintf('   evaluating LJ''s action failed: %s\n', ME.message);
        end

        ok = ~isnan(Jg);
        if sum(ok) < 3
            fprintf('   too few grid evaluations\n\n'); continue;
        end
        [Jb, ib] = min(Jg);
        V.u_best(q) = u_grid(ib);
        V.J_best(q) = Jb;
        V.J_lj(q) = J_at_lj;
        V.gap_abs(q) = J_at_lj - Jb;
        V.gap_rel(q) = (J_at_lj - Jb) / max(abs(Jb), eps);
        V.grids{q} = struct('u', u_grid, 'J', Jg, 'J_error', Jerr, ...
                            'J_psm', Jpsm, 'peak_Tp', Peak, 'dsm', Dsm);

        fprintf('\n   %8s %12s %12s %12s %11s\n', 'u [A]', 'J', 'J_error', 'J_psm', 'peak DSM');
        for i = 1:n_grid
            if isnan(Jg(i)), continue; end
            mk = '';
            if i == ib, mk = '  <-- grid optimum'; end
            if abs(u_grid(i) - u_lj) < (u_grid(2)-u_grid(1))/2
                mk = [mk '  <-- nearest to LJ'];
            end
            fprintf('   %8.5f %12.4e %12.4e %12.4e %+11.6f%s\n', ...
                u_grid(i), Jg(i), Jerr(i), Jpsm(i), Dsm(i), mk);
        end
        fprintf('\n   LJ:   u = %.5f  J = %.4e\n', u_lj, J_at_lj);
        fprintf('   best: u = %.5f  J = %.4e\n', u_grid(ib), Jb);
        fprintf('   gap:  %+.4e absolute, %+.2f%% relative\n', ...
            V.gap_abs(q), 100*V.gap_rel(q));
        fprintf('   distance in u: %+.5f A\n\n', u_lj - u_grid(ib));
    end

    %% ---- verdict ----
    fprintf('===== VERDICT =====\n');
    ok = ~isnan(V.gap_rel);
    if ~any(ok)
        V.conclusion = 'no query could be evaluated';
        fprintf('%s\n', V.conclusion); finish(V); return;
    end
    fprintf('%8s %11s %11s %13s %11s\n', 't_q', 'u_LJ', 'u_best', 'rel gap %', 'du');
    fprintf('%s\n', repmat('-', 1, 58));
    for q = find(ok)
        fprintf('%8.1f %11.5f %11.5f %13.2f %+11.5f\n', ...
            V.tq(q), V.u_lj(q), V.u_best(q), 100*V.gap_rel(q), V.u_lj(q)-V.u_best(q));
    end
    V.mean_gap = mean(V.gap_rel(ok));
    V.max_gap = max(V.gap_rel(ok));
    fprintf('\nmean relative gap: %+.2f%%   max: %+.2f%%\n', ...
        100*V.mean_gap, 100*V.max_gap);

    if V.max_gap > 0.10
        V.conclusion = sprintf(['LJ SUB-OPTIMALITY CONFIRMED: the applied action costs up ' ...
            'to %.1f%% more than a plain %d-point grid search at the same state. The late ' ...
            'drift is an optimizer artefact, so retune LJ for the widened range (smaller ' ...
            'gamma, relative early-stop threshold) rather than lowering lambda or the ' ...
            'pump limit.'], 100*V.max_gap, n_grid);
    elseif V.max_gap > 0.02
        V.conclusion = sprintf(['Mild sub-optimality: worst gap %.1f%%. LJ is close to the ' ...
            'grid optimum, so retuning would help but is unlikely to remove the drift on ' ...
            'its own.'], 100*V.max_gap);
    else
        V.conclusion = sprintf(['LJ is NOT the problem: worst gap %.2f%%, i.e. it finds ' ...
            'essentially the grid optimum. The drift must come from the cost function or ' ...
            'the plant dynamics -- revisit lambda/threshold or the horizon instead.'], ...
            100*V.max_gap);
    end
    fprintf('\n%s\n', V.conclusion);
    finish(V);
end

function finish(V)
    outfile = sprintf('verify_lj_suboptimality_%s.mat', datestr(now,'yyyymmdd_HHMMSS'));
    save(outfile, 'V');
    fprintf('\nSaved %s\n==================================\n\n', outfile);
end
