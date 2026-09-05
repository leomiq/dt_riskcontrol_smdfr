function A = check_severe_authority(paramfile)
% CHECK_SEVERE_AUTHORITY  In the severe scenario the peak PCT came out IDENTICAL
% at u_min and u_max (1195.7108 K to the 15th digit), i.e. the actuator appears
% to have lost all authority. Before running any (dsm_thr, lambda) sweep this
% must be explained: a sweep over a trajectory the control cannot move would be
% meaningless.
%
% Hypotheses:
%   H1  The predicted delta_tw_max (74.4 K, versus ~2 K nominal) dominates the
%       ROM output so completely that the coolant term is numerically
%       irrelevant -- the network extrapolated far outside its training range.
%   H2  The peak occurs in the WARM-UP window (before t_q), which is identical
%       for every u by construction, while the controllable tail differs.
%   H3  The ROM saturates (a limiter/clamp) at this temperature.
%
% H2 is the one that would still allow a meaningful study: the peak would be
% uncontrollable but the post-decision horizon would not.

    if nargin < 1 || isempty(paramfile)
        paramfile = 'severe_det_params_sevA40B15_A0.40_B15_vel0.090.mat';
    end

    A = struct();
    A.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fprintf('\n===== Severe-case control authority =====\n');
    fprintf('Started %s\nfile: %s\n\n', A.timestamp, paramfile);

    L = load(paramfile, 'det_params');
    dp = L.det_params;
    t_q     = double(dp.time_snapshot_root);
    t_sam   = double(dp.t_sam);
    delayed = double(dp.delayed_time_rom);
    N_pred  = size(dp.Y_denorm_50, 1);
    target  = CalculateTs(100, 543);
    U_Y = 1244; span = U_Y - target;
    time_delayed = t_q - delayed : t_sam : t_q - t_sam;

    fprintf('t_q %.1f, warm-up window [%.1f .. %.1f], horizon to %.1f\n', ...
        t_q, time_delayed(1), time_delayed(end), double(dp.tsim));
    fprintf('predicted delta_tw_max over horizon: min %.4f max %.4f K\n', ...
        min(dp.Y_denorm_50(:,2)), max(dp.Y_denorm_50(:,2)));
    fprintf('target %.4f K, U_y %.1f K\n\n', target, U_Y);

    d_mpc = ones(1, N_pred + 1);
    u_list = [56.69 56.71 56.73 56.75];
    A.u_list = u_list;
    A.full_peak = nan(1, numel(u_list));
    A.full_argmax_t = nan(1, numel(u_list));
    A.future_peak = nan(1, numel(u_list));
    A.future_mean = nan(1, numel(u_list));
    A.future_last = nan(1, numel(u_list));
    A.n_future = nan(1, numel(u_list));
    A.traj = cell(1, numel(u_list));
    A.times = cell(1, numel(u_list));

    for k = 1:numel(u_list)
        uv = u_list(k);
        u_trial = uv * ones(1, N_pred + 1);
        try
            [fd, ~] = DET_2(dp, u_trial, d_mpc);
            nd = fd{1}(1).tree_data;
            Tp = nd.data.pred_tw_max_signal.Data(:);
            ts = nd.data.pred_tw_max_signal.Time(:);
            A.traj{k} = Tp; A.times{k} = ts;

            [pk, im] = max(Tp);
            A.full_peak(k) = pk;
            A.full_argmax_t(k) = ts(im);

            fut = ts > time_delayed(end);
            A.n_future(k) = sum(fut);
            if any(fut)
                A.future_peak(k) = max(Tp(fut));
                A.future_mean(k) = mean(Tp(fut));
                A.future_last(k) = Tp(find(fut, 1, 'last'));
            end
            fprintf('u=%.5f : full peak %.6f K at t=%.2f | future peak %.6f, mean %.6f, last %.6f (n=%d)\n', ...
                uv, pk, ts(im), A.future_peak(k), A.future_mean(k), A.future_last(k), A.n_future(k));
        catch ME
            fprintf('u=%.5f FAILED: %s\n', uv, ME.message);
        end
    end

    %% ---- Diagnosis ----
    fprintf('\n-- Where is the peak? --\n');
    ok = ~isnan(A.full_peak);
    if any(ok)
        tpk = A.full_argmax_t(ok);
        in_warm = all(tpk <= time_delayed(end));
        fprintf('peak times: %s\n', mat2str(tpk, 8));
        if in_warm
            fprintf('>> H2 HOLDS: the peak sits INSIDE the warm-up window (<= %.1f),\n', time_delayed(end));
            fprintf('   which is identical for every u by construction. The peak is\n');
            fprintf('   therefore uncontrollable, but the FUTURE horizon may still be.\n');
            A.peak_in_warmup = true;
        else
            fprintf('>> peak is in the controllable horizon (> %.1f)\n', time_delayed(end));
            A.peak_in_warmup = false;
        end
    end

    fprintf('\n-- Authority on the FUTURE horizon only --\n');
    okf = ~isnan(A.future_peak);
    if sum(okf) >= 2
        A.future_peak_spread = max(A.future_peak(okf)) - min(A.future_peak(okf));
        A.future_mean_spread = max(A.future_mean(okf)) - min(A.future_mean(okf));
        A.future_last_spread = max(A.future_last(okf)) - min(A.future_last(okf));
        fprintf('future peak spread : %.8f K\n', A.future_peak_spread);
        fprintf('future mean spread : %.8f K\n', A.future_mean_spread);
        fprintf('future last spread : %.8f K\n', A.future_last_spread);
        A.future_dsm_lo = (min(A.future_peak(okf)) - target)/span;
        A.future_dsm_hi = (max(A.future_peak(okf)) - target)/span;
        fprintf('future peak DSM band: [%.6f, %.6f]\n', A.future_dsm_lo, A.future_dsm_hi);
    end

    %% ---- Verdict ----
    fprintf('\n===== VERDICT =====\n');
    if isfield(A,'future_peak_spread') && A.future_peak_spread > 1e-6
        A.conclusion = sprintf(['Authority survives on the future horizon: the peak DSM ' ...
            'there spans [%.6f, %.6f] (%.6f K of PCT). A (dsm_thr, lambda) sweep is ' ...
            'meaningful if the threshold is placed inside that band.'], ...
            A.future_dsm_lo, A.future_dsm_hi, A.future_peak_spread);
    else
        A.conclusion = ['NO authority anywhere: the predicted trajectory is identical for ' ...
            'every u. With delta_tw_max at ~74 K the network prediction dominates the ROM ' ...
            'output and the coolant term is numerically irrelevant -- the BiLSTM is ' ...
            'extrapolating far outside its training range, so this severe case is not ' ...
            'usable for a control study. Reduce the severity until authority reappears.'];
    end
    fprintf('%s\n', A.conclusion);

    outfile = sprintf('check_severe_authority_%s.mat', datestr(now,'yyyymmdd_HHMMSS'));
    save(outfile, 'A');
    fprintf('\nSaved %s\n=========================================\n\n', outfile);
end
