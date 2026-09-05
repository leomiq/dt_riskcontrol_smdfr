function A = audit_psm_future_only(matfile)
% AUDIT_PSM_FUTURE_ONLY  Verify that the DSM/J_psm penalty applies ONLY to the
% FUTURE propagated instants of each DET branch, never to past instants.
%
% Leo flagged this as a bug that existed previously and is believed fixed. It
% must be confirmed before any lambda sweep, because a penalty that also sees
% past instants would make every sweep result invalid: the optimizer would be
% charged for a trajectory it cannot influence.
%
% What "past" means here: the ROM is initialized delayed_time_rom seconds before
% the query instant (the warm-up window time_delayed = [t_root-delayed ..
% t_root-t_sam]) to avoid an initialization transient. Those samples precede the
% decision instant t_root, so u_trial cannot affect them and they must NOT enter
% J_error, J_psm or the DSM feasibility test.
%
% This runs DET_2 once (no optimizer) and checks the reconstructed time vector
% against t_root.

    if nargin < 1 || isempty(matfile)
        matfile = 'LJ_test_sampled_vel_ACS_2_loop_1180.mat';
    end

    A = struct();
    A.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fprintf('\n===== J_psm past/future audit =====\n');
    fprintf('Started %s\n\n', A.timestamp);

    S  = load(matfile);
    fprintf('Regenerating operating points...\n');
    dp = refresh_operating_points(S);
    if ~isfield(dp, 'Y_denorm_50'), dp.Y_denorm_50 = dp.Y_denorm; end

    t_root  = double(dp.time_snapshot_root);
    t_sam   = double(dp.t_sam);
    delayed = double(dp.delayed_time_rom);
    tsim    = double(dp.tsim);
    N_pred  = size(dp.Y_denorm_50, 1);
    target  = CalculateTs(100, 543);

    time_delayed = t_root - delayed : t_sam : t_root - t_sam;
    A.t_root = t_root;
    A.warmup_window = [time_delayed(1) time_delayed(end)];
    A.n_warmup = numel(time_delayed);
    fprintf('decision instant t_root   = %.1f s\n', t_root);
    fprintf('warm-up (PAST) window     = [%.1f .. %.1f]  (%d samples)\n', ...
        time_delayed(1), time_delayed(end), numel(time_delayed));
    fprintf('prediction (FUTURE) window= [%.1f .. %.1f]\n\n', t_root, tsim);

    u_trial = 56.72 * ones(1, N_pred + 1);
    d_mpc   = ones(1, N_pred + 1);

    fprintf('-- Running DET_2 --\n');
    [forest_data, forest_roots] = DET_2(dp, u_trial, d_mpc);
    fprintf('done\n\n');

    % Reproduce the reconstruction exactly as fitness_cost does, but keep the
    % time vector so past/future membership can be checked.
    single_tree = forest_data{1};
    root_init_state = single_tree(1).tree_data.init_state;
    Tp_rec = []; t_rec = []; owner = [];
    for i = 1:length(single_tree)
        nd = single_tree(i).tree_data;
        node_id = single_tree(i).tree.id;
        Tp_signal = nd.data.pred_tw_max_signal.Data(:);
        t_signal  = nd.data.pred_tw_max_signal.Time(:);
        if isempty(t_signal), continue; end
        if node_id == 0
            if isequal(root_init_state, [1,1])
                keep = t_signal > time_delayed(end);
            else
                keep = true(size(t_signal));
            end
        else
            keep = t_signal > time_delayed(end);
            if ~any(keep), keep = t_signal > t_signal(1); end
        end
        Tp_rec = cat(1, Tp_rec, Tp_signal(keep));
        t_rec  = cat(1, t_rec,  t_signal(keep));
        owner  = cat(1, owner, repmat(node_id, sum(keep), 1));
    end

    A.n_samples = numel(t_rec);
    A.t_min = min(t_rec);
    A.t_max = max(t_rec);
    fprintf('-- Reconstructed vector --\n');
    fprintf('samples   : %d\n', A.n_samples);
    fprintf('time range: [%.2f .. %.2f]\n', A.t_min, A.t_max);

    %% ---- THE CHECK ----
    fprintf('\n-- Past-instant contamination check --\n');
    past_mask = t_rec < t_root;
    A.n_past = sum(past_mask);
    A.n_future = sum(~past_mask);
    fprintf('samples STRICTLY BEFORE t_root (%.1f) : %d\n', t_root, A.n_past);
    fprintf('samples at or after t_root            : %d\n', A.n_future);

    warm_mask = t_rec >= time_delayed(1) & t_rec <= time_delayed(end);
    A.n_in_warmup = sum(warm_mask);
    fprintf('samples inside the warm-up window     : %d\n', A.n_in_warmup);

    if A.n_past == 0
        fprintf('\n>> CLEAN: no past instants enter the cost. The penalty applies to\n');
        fprintf('   the propagated FUTURE horizon only. Lambda sweeps are valid.\n');
        A.clean = true;
    else
        fprintf('\n>> CONTAMINATED: %d past samples reach the cost.\n', A.n_past);
        fprintf('   Times: %s\n', mat2str(t_rec(past_mask)', 8));
        fprintf('   These precede the decision instant, so u_trial cannot affect\n');
        fprintf('   them; J_error/J_psm/DSM would charge the optimizer for a\n');
        fprintf('   trajectory it does not control. FIX BEFORE ANY LAMBDA SWEEP.\n');
        A.clean = false;
        A.past_times = t_rec(past_mask)';
    end

    %% ---- Where does the penalty stand today? ----
    fprintf('\n-- Penalty activation status (current settings) --\n');
    U_y = 1244;
    A.peak_Tp = max(Tp_rec);
    A.margin_to_limit = U_y - A.peak_Tp;
    nominal = target;
    DSM = (Tp_rec - nominal) ./ (U_y - nominal);
    A.max_DSM = max(DSM);
    A.min_DSM = min(DSM);
    fprintf('target (setpoint)   : %.4f K\n', nominal);
    fprintf('peak Tp             : %.4f K\n', A.peak_Tp);
    fprintf('absolute limit U_y  : %.1f K\n', U_y);
    fprintf('margin to limit     : %.4f K\n', A.margin_to_limit);
    fprintf('DSM = (Tp-target)/(U_y-target)\n');
    fprintf('  max DSM           : %.6f   (penalty needs DSM > 1 to fire)\n', A.max_DSM);
    fprintf('  min DSM           : %.6f\n', A.min_DSM);

    % What threshold WOULD activate the penalty?
    fprintf('\n-- What it takes to activate J_psm --\n');
    fprintf('J_psm = lambda_pct * sum(max(0, Tp - U_y)^2), so it is identically\n');
    fprintf('zero while peak Tp (%.4f) stays below U_y (%.1f).\n', A.peak_Tp, U_y);
    A.U_y_activation = A.peak_Tp;
    fprintf('>> U_y would have to drop below %.4f K for the penalty to engage\n', A.peak_Tp);
    fprintf('   at this trajectory. Raising lambda_pct alone changes NOTHING\n');
    fprintf('   while the violation term is zero (lambda * 0 = 0).\n');

    % Effective threshold expressed as a DSM cap
    for cap = [0.99 0.5 0.2 0.1 0.05 0.02 0.01]
        U_eff = nominal + cap*(U_y - nominal);
        n_viol = sum(Tp_rec > U_eff);
        pen = sum(max(0, Tp_rec - U_eff).^2);
        fprintf('   DSM cap %.2f -> U_eff %9.4f K : %3d samples violate, penalty term %.6e\n', ...
            cap, U_eff, n_viol, pen);
    end

    outfile = sprintf('audit_psm_future_only_%s.mat', datestr(now,'yyyymmdd_HHMMSS'));
    save(outfile, 'A', 't_rec', 'Tp_rec', 'owner');
    fprintf('\nSaved %s\n===================================\n\n', outfile);
end
