function V = verify_lj_retuned(tq_list, n_grid)
% VERIFY_LJ_RETUNED  Does the retuned LJ now find the optimum on the queries
% where the original settings failed?
%
% Runs the ACTUAL optimizer (not a grid) at each instant with the retuned
% settings, then compares its answer against a fine grid reference. The original
% settings produced relative gaps of 161% to 18932% here
% (verify_lj_suboptimality.m).
%
% Changes under test, all in response to the widened bound:
%   1. luus_jaakola_det.m pass reset: r = min(|u_best - old_u_best|,
%      r_init*gamma^(pass-1)) instead of r = |u_best - old_u_best|. The original
%      reset discarded the within-pass contraction, so the radius never became
%      fine enough to resolve a ~0.05 A wide minimum.
%   2. early stop: relative (improvement < 0.1% of f_best) instead of absolute
%      (improvement < 1), which was meaningless at late-horizon costs of ~3000.
%   3. gamma 0.95 -> 0.90 and R 15 -> 24 in run_phase2c_full_loop.m.
%
%   tq_list  instants to probe (default [1176 1180 1184 1188])
%   n_grid   grid points for the reference (default 22)

    if nargin < 1 || isempty(tq_list), tq_list = [1176 1180 1184 1188]; end
    if nargin < 2 || isempty(n_grid),  n_grid = 22; end

    V = struct();
    V.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fprintf('\n===== Retuned LJ verification =====\n');
    fprintf('Started %s\n', V.timestamp);

    try, suppress_simulink_prompts(); catch, end

    dsm_thr = 0.00237737; lambda = 1e6; U_Y = 1244;
    u_min = 56.69; u_max = 56.90;
    target = CalculateTs(100, 543);
    span = U_Y - target;
    u_grid = linspace(u_min, u_max, n_grid);

    % retuned LJ settings
    max_passes = 30; M = 5; R = 24; gamma = 0.90;
    r_init = 0.5*(u_max - u_min);
    fprintf('LJ: R = %d, gamma = %.2f, r_init = %.5f, max_passes = %d\n', ...
        R, gamma, r_init, max_passes);
    fprintf('reference grid: %d points over [%.2f, %.2f]\n\n', n_grid, u_min, u_max);

    % what the ORIGINAL settings had produced, for the comparison table
    orig = containers.Map({1176, 1180, 1184, 1188}, ...
                          {56.78838, 56.69000, 56.84061, 56.69000});

    V.tq = tq_list;
    V.u_lj = nan(size(tq_list)); V.J_lj = nan(size(tq_list));
    V.u_best = nan(size(tq_list)); V.J_best = nan(size(tq_list));
    V.gap_rel = nan(size(tq_list));
    V.u_orig = nan(size(tq_list)); V.J_orig = nan(size(tq_list));
    V.passes = nan(size(tq_list));

    for q = 1:numel(tq_list)
        tq = tq_list(q);
        fprintf('---- t_q = %.1f s ----\n', tq);
        wsfile = sprintf('LJ_test_sampled_vel_ACS_2_loop_%g.mat', tq);
        if ~isfile(wsfile)
            fprintf('   %s absent; skipping\n\n', wsfile); continue;
        end
        S = load(wsfile);
        dp = refresh_operating_points(S);
        if ~isfield(dp, 'Y_denorm_50'), dp.Y_denorm_50 = dp.Y_denorm; end
        dp.dsm_thr = dsm_thr; dp.lambda_dsm = lambda;

        t_sam = double(dp.t_sam);
        delayed = double(dp.delayed_time_rom);
        N_pred = size(dp.Y_denorm_50, 1);
        d_mpc = ones(1, N_pred + 1);
        u_bounds = repmat([u_min, u_max], N_pred+1, 1);
        u0 = (u_min + r_init) * ones(1, N_pred+1);

        % --- reference grid ---
        Jg = nan(1, n_grid);
        for i = 1:n_grid
            ut = u_grid(i) * ones(1, N_pred + 1);
            try
                [fd, fr] = DET_2(dp, ut, d_mpc);
                opts = struct('dsm_thr', dsm_thr, 'lambda_pct', lambda, 'U_y', U_Y);
                [~, Jg(i)] = fitness_cost_dsm(fd, fr, N_pred, target, ...
                    tq, delayed, t_sam, ut, opts);
            catch
            end
        end
        [Jb, ib] = min(Jg);
        V.u_best(q) = u_grid(ib); V.J_best(q) = Jb;
        fprintf('   grid optimum: u = %.5f, J = %.4e\n', u_grid(ib), Jb);

        % --- retuned LJ, the real optimizer ---
        rng(1, 'twister');
        tt = tic;
        try
            [~, u_opt, f_opt, hist, tvec] = luus_jaakola_det( ...
                u0, u_bounds, max_passes, M, r_init, gamma, R, ...
                dp, target, d_mpc, tq, delayed, t_sam);
            V.u_lj(q) = u_opt(1);
            V.J_lj(q) = f_opt;
            V.passes(q) = sum(tvec > 0);
            fprintf('   retuned LJ : u = %.5f, J = %.4e (%d passes, %.1f s)\n', ...
                u_opt(1), f_opt, V.passes(q), toc(tt));
        catch ME
            fprintf('   retuned LJ FAILED: %s\n', ME.message);
            continue;
        end

        % --- what the original settings had chosen, costed the same way ---
        if isKey(orig, tq)
            V.u_orig(q) = orig(tq);
            try
                ut = V.u_orig(q) * ones(1, N_pred + 1);
                [fd, fr] = DET_2(dp, ut, d_mpc);
                opts = struct('dsm_thr', dsm_thr, 'lambda_pct', lambda, 'U_y', U_Y);
                [~, V.J_orig(q)] = fitness_cost_dsm(fd, fr, N_pred, target, ...
                    tq, delayed, t_sam, ut, opts);
                fprintf('   original LJ: u = %.5f, J = %.4e\n', V.u_orig(q), V.J_orig(q));
            catch
            end
        end

        V.gap_rel(q) = (V.J_lj(q) - Jb) / max(abs(Jb), eps);
        fprintf('   gap vs grid: %+.2f%%   (original was %+.2f%%)\n\n', ...
            100*V.gap_rel(q), 100*(V.J_orig(q)-Jb)/max(abs(Jb),eps));
    end

    %% ---- verdict ----
    fprintf('===== VERDICT =====\n');
    ok = ~isnan(V.gap_rel);
    fprintf('%8s %11s %11s %11s %13s %13s\n', ...
        't_q', 'u_grid', 'u_retuned', 'u_original', 'gap retuned', 'gap original');
    fprintf('%s\n', repmat('-', 1, 72));
    for q = find(ok)
        go = 100*(V.J_orig(q)-V.J_best(q))/max(abs(V.J_best(q)),eps);
        fprintf('%8.1f %11.5f %11.5f %11.5f %12.2f%% %12.2f%%\n', ...
            V.tq(q), V.u_best(q), V.u_lj(q), V.u_orig(q), 100*V.gap_rel(q), go);
    end
    if any(ok)
        V.mean_gap = mean(V.gap_rel(ok));
        V.max_gap = max(V.gap_rel(ok));
        fprintf('\nretuned: mean gap %+.2f%%, max %+.2f%%\n', ...
            100*V.mean_gap, 100*V.max_gap);
        if V.max_gap < 0.05
            V.conclusion = sprintf(['RETUNING WORKS: worst gap now %.2f%% (was up to ' ...
                '18932%%). LJ resolves the interior minimum, so the closed-loop drift ' ...
                'should disappear without touching lambda or the pump bound.'], ...
                100*V.max_gap);
        elseif V.max_gap < 0.5
            V.conclusion = sprintf(['Substantially better: worst gap %.1f%% (was up to ' ...
                '18932%%), but not yet tight. Consider gamma 0.85 or more passes before ' ...
                'committing to the full realization.'], 100*V.max_gap);
        else
            V.conclusion = sprintf(['Still sub-optimal: worst gap %.1f%%. The radius ' ...
                'schedule is not the only issue -- inspect the within-pass contraction ' ...
                'and R.'], 100*V.max_gap);
        end
        fprintf('\n%s\n', V.conclusion);
    end
    outfile = sprintf('verify_lj_retuned_%s.mat', datestr(now,'yyyymmdd_HHMMSS'));
    save(outfile, 'V');
    fprintf('\nSaved %s\n===================================\n\n', outfile);
end
