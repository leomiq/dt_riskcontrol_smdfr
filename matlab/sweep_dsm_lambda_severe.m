function W = sweep_dsm_lambda_severe(paramfile, n_u)
% SWEEP_DSM_LAMBDA_SEVERE  Joint (dsm_thr, lambda) sweep on a PREGENERATED
% severe scenario, i.e. one whose det_params already came out of
% generate_severe_scenario (HRM -> BiLSTM -> ROM).
%
% Difference from sweep_dsm_lambda.m: that one calls refresh_operating_points on
% an archived workspace; this one loads det_params directly from a
% severe_det_params_*.mat file.
%
% Difference in the metric that matters: in these severe cases the GLOBAL peak
% always lands at t = t_q - delayed (the ROM initialisation transient in the
% warm-up window), which is identical for every u by construction and therefore
% uncontrollable. The quantity the controller can actually shape is the peak
% over the FUTURE horizon (t > t_q - t_sam), so thresholds are placed inside the
% reachable FUTURE-DSM band and the penalty is evaluated on the reconstructed
% (future-only) trajectory that fitness_cost_dsm returns.

    if nargin < 1 || isempty(paramfile)
        paramfile = 'severe_det_params_scanA283_A0.28_B15_vel0.090.mat';
    end
    if nargin < 2 || isempty(n_u), n_u = 11; end

    W = struct();
    W.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    W.paramfile = paramfile;
    fprintf('\n===== Joint (DSM, lambda) sweep on a severe scenario =====\n');
    fprintf('Started %s\nfile: %s\n\n', W.timestamp, paramfile);

    L = load(paramfile, 'det_params', 'G');
    dp = L.det_params;
    if isfield(L, 'G') && isfield(L.G, 'cfg')
        W.cfg = L.G.cfg;
        fprintf('scenario: A = %.4f, B = %.1f, vel_ACS = %.4f, t_query = %.1f\n', ...
            W.cfg.A, W.cfg.B, W.cfg.vel_ACS, W.cfg.t_query);
    end

    t_q     = double(dp.time_snapshot_root);
    t_sam   = double(dp.t_sam);
    delayed = double(dp.delayed_time_rom);
    N_pred  = size(dp.Y_denorm_50, 1);
    target  = CalculateTs(100, 543);
    U_Y = 1244; span = U_Y - target;

    W.dtw_max = max(dp.Y_denorm_50(:,2));
    fprintf('predicted delta_tw_max over horizon: %.4f K', W.dtw_max);
    if W.dtw_max > 10
        fprintf('   <<< NOTE: nominal is ~2 K, so the BiLSTM is extrapolating;\n');
        fprintf('        results must be reported as a hypothetical case.\n');
    else
        fprintf('   (close to nominal ~2 K: prediction inside training range)\n');
    end
    fprintf('target %.4f K, U_y %.1f K, DSM span %.4f K\n\n', target, U_Y, span);

    u_min = 56.69; u_max = 56.75;
    u_grid = linspace(u_min, u_max, n_u);
    d_mpc  = ones(1, N_pred + 1);
    W.u_grid = u_grid;

    %% ---- Stage 1: cache one DET trajectory per candidate current ----
    fprintf('-- Stage 1: caching DET trajectories (%d evaluations) --\n', n_u);
    traj = cell(1, n_u);
    W.J_error = nan(1, n_u);
    W.peak_Tp = nan(1, n_u);
    W.peak_DSM = nan(1, n_u);
    W.mdot = nan(1, n_u);
    for k = 1:n_u
        uv = u_grid(k);
        u_trial = uv * ones(1, N_pred + 1);
        tt = tic;
        try
            [fd, fr] = DET_2(dp, u_trial, d_mpc);
            opts = struct('dsm_thr', 1.0, 'lambda_pct', 0, 'U_y', U_Y);
            [Tp, ~, ~, tm] = fitness_cost_dsm(fd, fr, N_pred, target, ...
                t_q, delayed, t_sam, u_trial, opts);
            traj{k} = Tp;                       % future-horizon only
            W.J_error(k)  = tm.J_error;
            W.peak_Tp(k)  = tm.peak_Tp;
            W.peak_DSM(k) = tm.peak_DSM;
            try, m = CurrentToMassFlow(u_trial); W.mdot(k) = m(1); catch, end
            fprintf('   u=%.5f -> future peak Tp %.4f K, DSM %.6f, J_error %.6e, mdot %.4f (%.1f s)\n', ...
                uv, tm.peak_Tp, tm.peak_DSM, tm.J_error, W.mdot(k), toc(tt));
        catch ME
            fprintf('   u=%.5f FAILED: %s\n', uv, ME.message);
        end
    end

    valid = ~cellfun(@isempty, traj);
    if sum(valid) < 3
        fprintf('\nToo few successful evaluations. Aborting.\n');
        W.conclusion = 'inconclusive'; finish(W); return;
    end

    dsm_lo = min(W.peak_DSM(valid)); dsm_hi = max(W.peak_DSM(valid));
    W.reachable_dsm = [dsm_lo dsm_hi];
    W.pct_authority = max(W.peak_Tp(valid)) - min(W.peak_Tp(valid));
    fprintf('\nreachable FUTURE peak-DSM band: [%.6f, %.6f] (width %.6f)\n', ...
        dsm_lo, dsm_hi, dsm_hi - dsm_lo);
    fprintf('PCT authority on the future horizon: %.6f K\n', W.pct_authority);

    vi = find(valid);
    [~, i0] = min(W.J_error(valid));
    W.u_baseline = u_grid(vi(i0));
    W.dsm_baseline = W.peak_DSM(vi(i0));
    W.saturated_baseline = (vi(i0) == n_u) || (vi(i0) == 1);
    fprintf('unpenalised optimum: u = %.5f A (DSM %.6f, J_error %.6e)%s\n\n', ...
        W.u_baseline, W.dsm_baseline, W.J_error(vi(i0)), ...
        tern(W.saturated_baseline, '  <<< at a bound', ''));

    %% ---- Stage 2: analytic sweep, thresholds inside the reachable band ----
    w = dsm_hi - dsm_lo;
    thr_grid = unique(round([dsm_lo - 0.25*w, dsm_lo, dsm_lo + 0.25*w, ...
        dsm_lo + 0.5*w, dsm_lo + 0.75*w, dsm_hi, dsm_hi + 0.25*w], 8));
    lam_grid = [0 1e0 1e1 1e2 1e3 1e4 1e5];
    W.thr_grid = thr_grid; W.lam_grid = lam_grid;
    fprintf('-- Stage 2: analytic (dsm_thr, lambda) sweep --\n');
    fprintf('thresholds (inside/around the reachable band): %s\n', mat2str(thr_grid, 6));
    fprintf('lambdas: %s\n\n', mat2str(lam_grid, 4));

    nT = numel(thr_grid); nL = numel(lam_grid);
    W.u_opt = nan(nT,nL); W.shift_A = nan(nT,nL); W.mdot_opt = nan(nT,nL);
    W.dsm_opt = nan(nT,nL); W.psm_share = nan(nT,nL); W.n_active = nan(nT,nL);
    for a = 1:nT
        for b = 1:nL
            thr = thr_grid(a); lam = lam_grid(b);
            Jt = nan(1,n_u); Jp = nan(1,n_u); nv = zeros(1,n_u);
            for k = 1:n_u
                if isempty(traj{k}), continue; end
                dsm = (traj{k} - target)./span;
                exc = max(0, dsm - thr);
                Jp(k) = lam*sum(exc.^2);
                Jt(k) = W.J_error(k) + Jp(k);
                nv(k) = sum(dsm > thr);
            end
            [Jm, im] = min(Jt);
            W.u_opt(a,b)     = u_grid(im);
            W.shift_A(a,b)   = u_grid(im) - W.u_baseline;
            W.mdot_opt(a,b)  = W.mdot(im);
            W.dsm_opt(a,b)   = W.peak_DSM(im);
            W.psm_share(a,b) = Jp(im)/max(Jm, eps);
            W.n_active(a,b)  = nv(im);
        end
    end

    fprintf('OPTIMAL CURRENT u*  [A]\n');       pm(W.u_opt, thr_grid, lam_grid, '%10.5f');
    fprintf('\nSHIFT vs unpenalised  [A]  (positive = MORE coolant)\n');
    pm(W.shift_A, thr_grid, lam_grid, '%+10.5f');
    fprintf('\nCOOLANT MASS FLOW at u*\n');     pm(W.mdot_opt, thr_grid, lam_grid, '%10.4f');
    fprintf('\nFUTURE peak DSM at u*\n');       pm(W.dsm_opt, thr_grid, lam_grid, '%10.6f');
    fprintf('\nJ_psm SHARE of J at u*\n');      pm(W.psm_share, thr_grid, lam_grid, '%10.4f');
    fprintf('\nVIOLATING SAMPLES at u*\n');     pm(W.n_active, thr_grid, lam_grid, '%10.0f');

    %% ---- Verdict ----
    fprintf('\n===== VERDICT =====\n');
    moved = abs(W.shift_A) > 1e-9;
    W.n_moving = sum(moved(:));
    fprintf('(dsm_thr, lambda) combinations moving the optimum: %d of %d\n', ...
        W.n_moving, numel(moved));
    if W.n_moving > 0
        [ia, ib] = find(moved);
        sh = arrayfun(@(p) W.shift_A(ia(p), ib(p)), 1:numel(ia));
        [~, bi] = max(sh);                       % largest move toward MORE coolant
        W.best = struct('thr', thr_grid(ia(bi)), 'lam', lam_grid(ib(bi)), ...
            'shift', sh(bi), 'u', W.u_opt(ia(bi),ib(bi)), ...
            'mdot', W.mdot_opt(ia(bi),ib(bi)), 'dsm', W.dsm_opt(ia(bi),ib(bi)), ...
            'share', W.psm_share(ia(bi),ib(bi)));
        fprintf('\nlargest move toward more coolant:\n');
        fprintf('   dsm_thr %.6f, lambda %.0e\n', W.best.thr, W.best.lam);
        fprintf('   u* %+.5f A : %.5f -> %.5f A\n', W.best.shift, W.u_baseline, W.best.u);
        fprintf('   coolant mass flow %.4f  (baseline %.4f)\n', W.best.mdot, ...
            W.mdot(u_grid == W.u_baseline));
        fprintf('   future peak DSM %.6f -> %.6f\n', W.dsm_baseline, W.best.dsm);
        fprintf('   J_psm share of cost: %.1f%%\n', 100*W.best.share);
        W.conclusion = sprintf(['DSM penalty reorients the decision: dsm_thr %.6f with ' ...
            'lambda %.0e moves u* by %+.5f A (mass flow %.4f -> %.4f), taking the future ' ...
            'peak DSM from %.6f to %.6f with J_psm at %.1f%% of the cost.'], ...
            W.best.thr, W.best.lam, W.best.shift, W.mdot(u_grid == W.u_baseline), ...
            W.best.mdot, W.dsm_baseline, W.best.dsm, 100*W.best.share);
    else
        W.conclusion = sprintf(['No combination moves the optimum. Unpenalised optimum is ' ...
            'u = %.5f%s; if it already sits at a bound, tracking and safety agree and the ' ...
            'penalty is redundant here.'], W.u_baseline, ...
            tern(W.saturated_baseline, ' (at a bound)', ''));
    end
    fprintf('\n%s\n', W.conclusion);
    finish(W);
end

function pm(M, thr, lam, fmt)
    fprintf('%12s', 'thr \ lambda');
    for b = 1:numel(lam), fprintf('%10.0e', lam(b)); end
    fprintf('\n');
    for a = 1:size(M,1)
        fprintf('%12.6f', thr(a));
        for b = 1:size(M,2), fprintf(fmt, M(a,b)); end
        fprintf('\n');
    end
end

function s = tern(c,a,b)
    if c, s = a; else, s = b; end
end

function finish(W)
    outfile = sprintf('sweep_dsm_lambda_severe_%s.mat', datestr(now,'yyyymmdd_HHMMSS'));
    save(outfile, 'W');
    fprintf('\nSaved %s\n=========================================\n\n', outfile);
end
