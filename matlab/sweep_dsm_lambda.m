function W = sweep_dsm_lambda(matfile, n_u)
% SWEEP_DSM_LAMBDA  Joint (dsm_thr, lambda_pct) sweep showing whether the
% deterministic safety-margin penalty reorients the MPC decision toward more
% ACS coolant injection.
%
% Purpose (Leo): demonstrate that the margin activates in critical situations
% and smooths the control action -- i.e. that lowering the DSM threshold while
% raising lambda shifts the optimal current toward more coolant, keeping the PCT
% inside the imposed margin.
%
% Method, and why it is cheap: the DET trajectory for a given u_trial does not
% depend on (dsm_thr, lambda) at all -- those only enter the cost algebra. So we
% simulate the DET ONCE per candidate current (n_u evaluations, ~0.5-3 s each)
% and then evaluate every (dsm_thr, lambda) combination analytically on the
% cached trajectories. A 6x6 sweep over 11 currents costs 11 simulations, not
% 396.
%
% Reports, for each (dsm_thr, lambda):
%   - the argmin current and its coolant mass flow
%   - the shift in argmin relative to the unpenalised baseline
%   - whether the penalty is active (violating samples > 0) and its share of J

    if nargin < 1 || isempty(matfile), matfile = 'LJ_test_sampled_vel_ACS_2_loop_1180.mat'; end
    if nargin < 2 || isempty(n_u),     n_u = 11; end

    W = struct();
    W.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fprintf('\n========= Joint (DSM threshold, lambda) sweep =========\n');
    fprintf('Started %s\n\n', W.timestamp);

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
    span    = U_Y - target;

    u_min = 56.69; u_max = 56.75;
    u_grid = linspace(u_min, u_max, n_u);
    d_mpc  = ones(1, N_pred + 1);

    W.matfile = matfile; W.t_root = t_root; W.target = target;
    W.U_Y = U_Y; W.dsm_span = span; W.u_grid = u_grid;
    fprintf('query state t_root = %.1f s\n', t_root);
    fprintf('u grid      : %d points in [%.5f, %.5f] A\n', n_u, u_min, u_max);
    fprintf('target      : %.4f K,  DSM span %.4f K\n\n', target, span);

    %% ---- Stage 1: simulate the DET once per candidate current ----
    fprintf('-- Stage 1: caching DET trajectories (%d evaluations) --\n', n_u);
    traj = cell(1, n_u);
    W.mdot = nan(1, n_u);
    W.peak_Tp = nan(1, n_u);
    W.peak_DSM = nan(1, n_u);
    W.J_error = nan(1, n_u);
    for k = 1:n_u
        uv = u_grid(k);
        u_trial = uv * ones(1, N_pred + 1);
        tt = tic;
        try
            [fd, fr] = DET_2(dp, u_trial, d_mpc);
            opts = struct('dsm_thr', 1.0, 'lambda_pct', 0, 'verbose', false, 'U_y', U_Y);
            [Tp, ~, ~, tm] = fitness_cost_dsm(fd, fr, N_pred, target, ...
                t_root, delayed, t_sam, u_trial, opts);
            traj{k} = Tp;
            W.J_error(k) = tm.J_error;
            W.peak_Tp(k) = tm.peak_Tp;
            W.peak_DSM(k) = tm.peak_DSM;
            try
                m = CurrentToMassFlow(u_trial);
                W.mdot(k) = m(1);
            catch
            end
            fprintf('   u=%.5f -> peak Tp %.4f K, DSM %.6f, J_error %.6e, mdot %.4f (%.1f s)\n', ...
                uv, tm.peak_Tp, tm.peak_DSM, tm.J_error, W.mdot(k), toc(tt));
        catch ME
            fprintf('   u=%.5f FAILED: %s\n', uv, ME.message);
        end
    end

    valid = ~cellfun(@isempty, traj);
    if sum(valid) < 3
        fprintf('\nToo few successful evaluations (%d). Aborting.\n', sum(valid));
        W.conclusion = 'inconclusive'; finish(W); return;
    end

    dsm_lo = min(W.peak_DSM(valid)); dsm_hi = max(W.peak_DSM(valid));
    W.reachable_dsm = [dsm_lo dsm_hi];
    fprintf('\nreachable peak-DSM band: [%.6f, %.6f]  (width %.6f)\n', ...
        dsm_lo, dsm_hi, dsm_hi - dsm_lo);
    fprintf('PCT authority over the grid: %.6f K\n', ...
        max(W.peak_Tp(valid)) - min(W.peak_Tp(valid)));

    %% ---- Baseline: unpenalised optimum ----
    [~, i0] = min(W.J_error(valid));
    vi = find(valid);
    W.u_baseline = u_grid(vi(i0));
    W.dsm_baseline = W.peak_DSM(vi(i0));
    fprintf('\nunpenalised optimum: u = %.5f A (peak DSM %.6f, J_error %.6e)\n', ...
        W.u_baseline, W.dsm_baseline, W.J_error(vi(i0)));

    %% ---- Stage 2: analytic sweep over (dsm_thr, lambda) ----
    % Thresholds spanning the reachable band, so some bind and some do not.
    thr_grid = unique(round([dsm_lo*0.5, dsm_lo, ...
        dsm_lo + 0.25*(dsm_hi-dsm_lo), dsm_lo + 0.5*(dsm_hi-dsm_lo), ...
        dsm_lo + 0.75*(dsm_hi-dsm_lo), dsm_hi, dsm_hi*1.5], 8));
    lam_grid = [0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6];
    W.thr_grid = thr_grid; W.lam_grid = lam_grid;

    fprintf('\n-- Stage 2: analytic (dsm_thr, lambda) sweep --\n');
    fprintf('thresholds: %s\n', mat2str(thr_grid, 6));
    fprintf('lambdas   : %s\n\n', mat2str(lam_grid, 4));

    nT = numel(thr_grid); nL = numel(lam_grid);
    W.u_opt      = nan(nT, nL);
    W.dsm_opt    = nan(nT, nL);
    W.mdot_opt   = nan(nT, nL);
    W.J_opt      = nan(nT, nL);
    W.psm_share  = nan(nT, nL);
    W.n_active   = nan(nT, nL);
    W.shift_A    = nan(nT, nL);

    for a = 1:nT
        for b = 1:nL
            thr = thr_grid(a); lam = lam_grid(b);
            Jtot = nan(1, n_u); Jpsm = nan(1, n_u); nviol = zeros(1, n_u);
            for k = 1:n_u
                if isempty(traj{k}), continue; end
                Tp = traj{k};
                DSM = (Tp - target) ./ span;
                exc = max(0, DSM - thr);
                Jpsm(k) = lam * sum(exc.^2);
                Jtot(k) = W.J_error(k) + Jpsm(k);
                nviol(k) = sum(DSM > thr);
            end
            [Jm, im] = min(Jtot);
            W.u_opt(a,b)     = u_grid(im);
            W.dsm_opt(a,b)   = W.peak_DSM(im);
            W.mdot_opt(a,b)  = W.mdot(im);
            W.J_opt(a,b)     = Jm;
            W.psm_share(a,b) = Jpsm(im) / max(Jm, eps);
            W.n_active(a,b)  = nviol(im);
            W.shift_A(a,b)   = u_grid(im) - W.u_baseline;
        end
    end

    %% ---- Report ----
    fprintf('OPTIMAL CURRENT u*(dsm_thr, lambda)   [A]\n');
    print_matrix(W.u_opt, thr_grid, lam_grid, '%10.5f');
    fprintf('\nSHIFT vs unpenalised optimum   [A]  (positive = MORE coolant)\n');
    print_matrix(W.shift_A, thr_grid, lam_grid, '%+10.5f');
    fprintf('\nCOOLANT MASS FLOW at u*\n');
    print_matrix(W.mdot_opt, thr_grid, lam_grid, '%10.4f');
    fprintf('\nPEAK DSM at u*\n');
    print_matrix(W.dsm_opt, thr_grid, lam_grid, '%10.6f');
    fprintf('\nJ_psm SHARE of total J at u*\n');
    print_matrix(W.psm_share, thr_grid, lam_grid, '%10.4f');

    %% ---- Verdict ----
    fprintf('\n========= VERDICT =========\n');
    moved = abs(W.shift_A) > 1e-9;
    W.n_combos_moving = sum(moved(:));
    fprintf('(dsm_thr, lambda) combinations that move the optimum: %d of %d\n', ...
        W.n_combos_moving, numel(moved));

    if W.n_combos_moving > 0
        [ia, ib] = find(moved);
        shifts = arrayfun(@(p) W.shift_A(ia(p), ib(p)), 1:numel(ia));
        [~, best] = max(abs(shifts));
        W.best_thr = thr_grid(ia(best));
        W.best_lam = lam_grid(ib(best));
        W.best_shift = shifts(best);
        W.best_u = W.u_opt(ia(best), ib(best));
        W.best_mdot = W.mdot_opt(ia(best), ib(best));
        W.best_dsm = W.dsm_opt(ia(best), ib(best));
        fprintf('\nlargest reorientation:\n');
        fprintf('   dsm_thr = %.6f, lambda = %.0e\n', W.best_thr, W.best_lam);
        fprintf('   u* moves %+.5f A : %.5f -> %.5f A\n', ...
            W.best_shift, W.u_baseline, W.best_u);
        fprintf('   coolant mass flow at u*  : %.4f\n', W.best_mdot);
        fprintf('   peak DSM %.6f -> %.6f\n', W.dsm_baseline, W.best_dsm);
        if W.best_shift > 0
            W.conclusion = sprintf(['The DSM penalty REORIENTS the decision toward more ' ...
                'coolant: with dsm_thr = %.6f and lambda = %.0e the optimum moves ' ...
                '%+.5f A (mass flow %.4f), lowering the peak DSM from %.6f to %.6f. ' ...
                'This demonstrates the safety margin activating and shaping the ' ...
                'control action.'], W.best_thr, W.best_lam, W.best_shift, ...
                W.best_mdot, W.dsm_baseline, W.best_dsm);
        else
            W.conclusion = sprintf(['The penalty moves the optimum by %+.5f A, i.e. toward ' ...
                'LESS coolant. Inspect: with a non-monotonic J_error the DSM term may be ' ...
                'trading tracking error against margin in the opposite direction.'], W.best_shift);
        end
    else
        W.conclusion = ['No (dsm_thr, lambda) combination moves the optimum. Either every ' ...
            'threshold sits outside the reachable DSM band, or J_error dominates at every ' ...
            'tested lambda.'];
    end
    fprintf('\n%s\n', W.conclusion);
    finish(W);
end

function print_matrix(M, thr, lam, fmt)
    fprintf('%12s', 'dsm_thr \\ l');
    for b = 1:numel(lam), fprintf('%10.0e', lam(b)); end
    fprintf('\n');
    for a = 1:size(M,1)
        fprintf('%12.6f', thr(a));
        for b = 1:size(M,2)
            fprintf(fmt, M(a,b));
        end
        fprintf('\n');
    end
end

function finish(W)
    outfile = sprintf('sweep_dsm_lambda_%s.mat', datestr(now,'yyyymmdd_HHMMSS'));
    save(outfile, 'W');
    fprintf('\nSaved %s\n===========================\n\n', outfile);
end
