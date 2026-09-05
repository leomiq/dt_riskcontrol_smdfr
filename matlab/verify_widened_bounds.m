function V = verify_widened_bounds(u_max_new, dsm_thr, lambda_pct, n_u)
% VERIFY_WIDENED_BOUNDS  Does the DSM penalty reorient the decision once the
% MHD pump is allowed to deliver more coolant?
%
% MOTIVATION. In the 61-query realization the penalty dominated the cost for 36
% queries (J_psm/J_error up to 197) yet the PCT trajectory barely moved. The
% reason is not a weak penalty: the control was SATURATED at u_max most of that
% window, so the safety term had already obtained everything it could ask for.
% Saturation was the binding constraint.
%
% The pump model (CurrentToMassFlow.m, ECM of Nowak et al. 2022) shows why the
% window is so tight:
%     zero-flow current  DeltaP/C = 56.6872 A
%     bounds used        [56.69, 56.75] -> only 0.0028 .. 0.0628 A above onset
%     slope              dmdot/dI ~ 204.7 kg/s per A (essentially linear)
% So [56.69, 56.75] gives just 12.86 kg/s of authority, and only because the
% bounds sit right at the pump's onset -- not because the pump is weak.
%
%     u_max 56.75 ->  12.86 kg/s   (current)
%     u_max 56.90 ->  43.55 kg/s   (3.4x)
%     u_max 57.00 ->  63.99 kg/s   (5.0x)
% The archived PID series already drove this same plant model up to 57.26 A, so
% widening toward ~57.0 A uses a region the model has already been through.
%
% This checks, at the t=1180 fixed state and with NO optimizer, whether the
% penalty now moves the argmin -- before committing ~8 h to two realizations.
%
%   u_max_new  widened upper bound (default 56.90)
%   dsm_thr    DSM threshold (default 0.00237737)
%   lambda_pct penalty weight (default 1e6)
%   n_u        grid points (default 13)

    if nargin < 1 || isempty(u_max_new),  u_max_new = 56.90;      end
    if nargin < 2 || isempty(dsm_thr),    dsm_thr = 0.00237737;   end
    if nargin < 3 || isempty(lambda_pct), lambda_pct = 1e6;       end
    if nargin < 4 || isempty(n_u),        n_u = 13;               end

    V = struct();
    V.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fprintf('\n===== Widened-bounds verification (fixed state) =====\n');
    fprintf('Started %s\n', V.timestamp);
    fprintf('u range: [56.69, %.2f]  (was [56.69, 56.75])\n', u_max_new);
    fprintf('dsm_thr = %.8f, lambda = %.0e, grid %d points\n\n', ...
        dsm_thr, lambda_pct, n_u);

    try, suppress_simulink_prompts(); catch, end

    S = load('LJ_test_sampled_vel_ACS_2_loop_1180.mat');
    dp = refresh_operating_points(S);
    if ~isfield(dp, 'Y_denorm_50'), dp.Y_denorm_50 = dp.Y_denorm; end

    t_q     = double(dp.time_snapshot_root);
    t_sam   = double(dp.t_sam);
    delayed = double(dp.delayed_time_rom);
    N_pred  = size(dp.Y_denorm_50, 1);
    target  = CalculateTs(100, 543);
    U_Y = 1244; span = U_Y - target;
    d_mpc = ones(1, N_pred + 1);

    u_min = 56.69;
    u_grid = linspace(u_min, u_max_new, n_u);
    V.u_grid = u_grid;
    V.mdot = arrayfun(@(u) firstflow(u, N_pred), u_grid);

    fprintf('-- Stage 1: DET trajectories over the widened grid --\n');
    traj = cell(1, n_u);
    V.J_error = nan(1, n_u); V.peak_Tp = nan(1, n_u); V.peak_DSM = nan(1, n_u);
    for k = 1:n_u
        uv = u_grid(k);
        u_trial = uv * ones(1, N_pred + 1);
        tt = tic;
        try
            [fd, fr] = DET_2(dp, u_trial, d_mpc);
            opts = struct('dsm_thr', 1.0, 'lambda_pct', 0, 'U_y', U_Y);
            [Tp, ~, ~, tm] = fitness_cost_dsm(fd, fr, N_pred, target, ...
                t_q, delayed, t_sam, u_trial, opts);
            traj{k} = Tp;
            V.J_error(k) = tm.J_error;
            V.peak_Tp(k) = tm.peak_Tp;
            V.peak_DSM(k) = tm.peak_DSM;
            fprintf('   u=%.5f  mdot=%8.3f  peak Tp=%.4f  DSM=%+.6f  J_error=%.6e  (%.1f s)\n', ...
                uv, V.mdot(k), tm.peak_Tp, tm.peak_DSM, tm.J_error, toc(tt));
        catch ME
            fprintf('   u=%.5f FAILED: %s\n', uv, ME.message);
        end
    end

    ok = ~cellfun(@isempty, traj);
    if sum(ok) < 3
        fprintf('\nToo few evaluations. Aborting.\n');
        V.conclusion = 'inconclusive'; finish(V); return;
    end
    vi = find(ok);

    fprintf('\n-- Authority with the widened bound --\n');
    V.pct_authority = max(V.peak_Tp(ok)) - min(V.peak_Tp(ok));
    V.dsm_band = [min(V.peak_DSM(ok)) max(V.peak_DSM(ok))];
    V.mdot_span = [min(V.mdot(ok)) max(V.mdot(ok))];
    fprintf('PCT authority    : %.6f K   (was 0.9749 K with u_max 56.75)\n', V.pct_authority);
    fprintf('reachable DSM    : [%+.6f, %+.6f]\n', V.dsm_band(1), V.dsm_band(2));
    fprintf('mass-flow range  : [%.3f, %.3f] kg/s\n', V.mdot_span(1), V.mdot_span(2));

    [~, i0] = min(V.J_error(ok));
    V.u_baseline = u_grid(vi(i0));
    V.dsm_baseline = V.peak_DSM(vi(i0));
    V.mdot_baseline = V.mdot(vi(i0));
    V.baseline_at_bound = (vi(i0) == 1) || (vi(i0) == n_u);
    fprintf('\nunpenalised optimum: u = %.5f A (mdot %.3f, DSM %+.6f)%s\n', ...
        V.u_baseline, V.mdot_baseline, V.dsm_baseline, ...
        tern(V.baseline_at_bound, '   <<< AT A BOUND', '   <<< INTERIOR'));

    %% ---- penalised optimum ----
    fprintf('\n-- Stage 2: with the DSM penalty --\n');
    Jt = nan(1, n_u); Jp = nan(1, n_u); nv = zeros(1, n_u);
    for k = 1:n_u
        if isempty(traj{k}), continue; end
        DSM = (traj{k} - target)./span;
        exc = max(0, DSM - dsm_thr);
        Jp(k) = lambda_pct * sum(exc.^2);
        Jt(k) = V.J_error(k) + Jp(k);
        nv(k) = sum(DSM > dsm_thr);
    end
    [Jm, im] = min(Jt);
    V.u_pen = u_grid(im);
    V.mdot_pen = V.mdot(im);
    V.dsm_pen = V.peak_DSM(im);
    V.psm_share = Jp(im)/max(Jm, eps);
    V.shift = V.u_pen - V.u_baseline;
    V.extra_mdot = V.mdot_pen - V.mdot_baseline;
    for k = 1:n_u
        if isnan(Jt(k)), continue; end
        mark = '';
        if k == im, mark = '  <<< penalised optimum'; end
        if u_grid(k) == V.u_baseline, mark = [mark '  (unpenalised optimum)']; end
        fprintf('   u=%.5f  J_error=%.4e  J_psm=%.4e  J=%.4e  viol=%3d%s\n', ...
            u_grid(k), V.J_error(k), Jp(k), Jt(k), nv(k), mark);
    end

    %% ---- Verdict ----
    fprintf('\n===== VERDICT =====\n');
    fprintf('unpenalised u* = %.5f A (mdot %.3f kg/s)\n', V.u_baseline, V.mdot_baseline);
    fprintf('penalised   u* = %.5f A (mdot %.3f kg/s)\n', V.u_pen, V.mdot_pen);
    fprintf('shift          = %+.5f A   extra coolant = %+.3f kg/s\n', V.shift, V.extra_mdot);
    fprintf('peak DSM       : %+.6f -> %+.6f\n', V.dsm_baseline, V.dsm_pen);
    fprintf('J_psm share    : %.1f%%\n', 100*V.psm_share);
    fprintf('penalised optimum at a bound: %d\n', (im == 1) || (im == n_u));

    if abs(V.shift) < 1e-9
        V.conclusion = ['NO reorientation even with the widened bound. The penalised and ' ...
            'unpenalised optima coincide, so the tight bound was not the limiting factor ' ...
            'at this operating point -- do not spend 8 h on the two realizations yet.'];
    elseif V.shift > 0
        V.conclusion = sprintf(['REORIENTATION CONFIRMED with the widened bound: the penalty ' ...
            'moves u* by %+.5f A, injecting %+.3f kg/s more coolant (%.2fx) and taking the ' ...
            'peak DSM from %+.6f to %+.6f. Widening the pump limit gives the safety term ' ...
            'room to act. Worth running the two full realizations.'], ...
            V.shift, V.extra_mdot, V.mdot_pen/max(V.mdot_baseline,eps), ...
            V.dsm_baseline, V.dsm_pen);
    else
        V.conclusion = sprintf(['The penalty moves u* by %+.5f A, i.e. toward LESS coolant. ' ...
            'Check the sign convention before running the realizations.'], V.shift);
    end
    fprintf('\n%s\n', V.conclusion);
    finish(V);
end

function m = firstflow(u, N_pred)
    m = NaN;
    try
        seq = CurrentToMassFlow(u * ones(1, N_pred + 1));
        m = seq(1);
    catch
    end
end

function s = tern(c, a, b)
    if c, s = a; else, s = b; end
end

function finish(V)
    outfile = sprintf('verify_widened_bounds_%s.mat', datestr(now,'yyyymmdd_HHMMSS'));
    save(outfile, 'V');
    fprintf('\nSaved %s\n=====================================\n\n', outfile);
end
