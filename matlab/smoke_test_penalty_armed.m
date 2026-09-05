function S = smoke_test_penalty_armed(dsm_thr, lambda_pct)
% SMOKE_TEST_PENALTY_ARMED  Confirm the DSM penalty actually fires inside
% fitness_cost before committing ~6 h to a full closed-loop realization.
%
% Checks, at the t=1180 fixed state and with NO optimizer:
%   1. with the global unset, J_psm is 0 (historical behaviour preserved)
%   2. with the global set, J_psm is non-zero and J differs
%   3. the argmin over a small current grid moves the way the fixed-state sweep
%      predicted (toward more coolant)
%
% Takes about a minute. If this passes, run_penalised_trajectory is safe to launch.

    if nargin < 1 || isempty(dsm_thr),    dsm_thr = 0.00237737; end
    if nargin < 2 || isempty(lambda_pct), lambda_pct = 1e4;     end

    S = struct();
    fprintf('\n===== Smoke test: is the DSM penalty armed? =====\n');
    fprintf('dsm_thr = %.8f, lambda_pct = %.0e\n\n', dsm_thr, lambda_pct);

    try, suppress_simulink_prompts(); catch, end

    W = load('LJ_test_sampled_vel_ACS_2_loop_1180.mat');
    dp = refresh_operating_points(W);
    if ~isfield(dp, 'Y_denorm_50'), dp.Y_denorm_50 = dp.Y_denorm; end
    t_q = double(dp.time_snapshot_root);
    t_sam = double(dp.t_sam);
    delayed = double(dp.delayed_time_rom);
    N_pred = size(dp.Y_denorm_50,1);
    target = CalculateTs(100,543);
    d_mpc = ones(1, N_pred+1);

    u_grid = [56.69 56.714 56.726 56.744 56.75];
    S.u_grid = u_grid;

    %% ---- Pass 1: penalty DISARMED ----
    % The settings ride inside det_params (not a global) because the real run
    % evaluates the cost inside a parfor, where client globals are invisible.
    fprintf('-- Pass 1: det_params without penalty fields (historical) --\n');
    dp_off = dp;
    if isfield(dp_off,'dsm_thr'),    dp_off = rmfield(dp_off,'dsm_thr');    end
    if isfield(dp_off,'lambda_dsm'), dp_off = rmfield(dp_off,'lambda_dsm'); end
    J_off = nan(1, numel(u_grid));
    for k = 1:numel(u_grid)
        u_trial = u_grid(k)*ones(1, N_pred+1);
        [fd, fr] = DET_2(dp_off, u_trial, d_mpc);
        [~, J_off(k)] = evalc_cost(fd, fr, N_pred, target, t_q, delayed, t_sam, u_trial);
        fprintf('   u=%.5f  J=%.6e\n', u_grid(k), J_off(k));
    end
    [~, i_off] = min(J_off);
    S.J_off = J_off; S.u_off = u_grid(i_off);
    fprintf('   argmin: u = %.5f\n\n', S.u_off);

    %% ---- Pass 2: penalty ARMED ----
    fprintf('-- Pass 2: det_params.dsm_thr / .lambda_dsm set (penalty active) --\n');
    dp_on = dp;
    dp_on.dsm_thr    = dsm_thr;
    dp_on.lambda_dsm = lambda_pct;
    J_on = nan(1, numel(u_grid));
    for k = 1:numel(u_grid)
        u_trial = u_grid(k)*ones(1, N_pred+1);
        [fd, fr] = DET_2(dp_on, u_trial, d_mpc);
        [~, J_on(k)] = evalc_cost(fd, fr, N_pred, target, t_q, delayed, t_sam, u_trial);
        fprintf('   u=%.5f  J=%.6e  (delta vs disarmed %+.6e)\n', ...
            u_grid(k), J_on(k), J_on(k)-J_off(k));
    end
    [~, i_on] = min(J_on);
    S.J_on = J_on; S.u_on = u_grid(i_on);
    fprintf('   argmin: u = %.5f\n\n', S.u_on);

    %% ---- Verdict ----
    fprintf('===== VERDICT =====\n');
    S.penalty_changes_J = any(abs(J_on - J_off) > 1e-12);
    S.argmin_moved = (S.u_on ~= S.u_off);
    S.shift = S.u_on - S.u_off;
    fprintf('penalty changes J : %d\n', S.penalty_changes_J);
    fprintf('argmin moved      : %d  (%.5f -> %.5f, %+.5f A)\n', ...
        S.argmin_moved, S.u_off, S.u_on, S.shift);

    if S.penalty_changes_J && S.shift > 0
        S.conclusion = sprintf(['ARMED AND EFFECTIVE: the penalty raises J and moves the ' ...
            'optimum %+.5f A toward more coolant. Safe to launch the full realization.'], S.shift);
    elseif S.penalty_changes_J
        S.conclusion = sprintf(['ARMED but the optimum moved %+.5f A. Check the sign before ' ...
            'committing 6 h.'], S.shift);
    else
        S.conclusion = ['NOT ARMED: J is unchanged. The global is not reaching fitness_cost ' ...
            '-- do not launch the realization.'];
    end
    fprintf('\n%s\n=================================\n\n', S.conclusion);
end

function [Tp, J] = evalc_cost(fd, fr, nu, target, t_q, delayed, t_sam, u_trial)
    % fitness_cost prints its own J breakdown line; evalc keeps it out of the
    % smoke-test log. The captured text is discarded.
    evalc('[Tp, J] = fitness_cost(fd, fr, nu, target, t_q, delayed, t_sam, u_trial);');
end
