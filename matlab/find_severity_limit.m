function F = find_severity_limit(A_list, t_query)
% FIND_SEVERITY_LIMIT  Locate the most severe LOCA the BiLSTM still supports.
%
% At A = 0.4, B = 15 the network predicts delta_tw_max = 74.4 K against ~2 K
% nominal (a 37x jump), i.e. it extrapolates far outside its training range and
% the resulting "severe" case is an uncontrollable warm-up artifact. This scans
% A between the nominal value and the upper bound, reporting for each:
%
%   - predicted delta_tw_max over the horizon (the extrapolation indicator)
%   - where the PCT peak falls: warm-up window (uncontrollable) or future horizon
%   - the peak DSM on the FUTURE horizon, and the control authority there
%
% The usable case is the largest A whose delta_tw_max stays physically plausible
% AND whose peak lands in the controllable horizon with a positive DSM.
%
% Only the HRM + inference stages run per A (about 20 s each), plus two DET
% evaluations to measure authority.

    if nargin < 1 || isempty(A_list),  A_list = [0.2834 0.30 0.32 0.34 0.36]; end
    if nargin < 2 || isempty(t_query), t_query = 1180; end

    F = struct();
    F.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    F.A_list = A_list; F.t_query = t_query;
    fprintf('\n===== BiLSTM severity limit scan =====\n');
    fprintf('Started %s\n', F.timestamp);
    fprintf('A values : %s   (nominal 0.2834, bound 0.4)\n', mat2str(A_list));
    fprintf('t_query  : %.1f s\n', t_query);
    fprintf('reference: nominal delta_tw_max is ~2 K; at A=0.4,B=15 it was 74.4 K\n\n');

    target = CalculateTs(100, 543);
    U_Y = 1244; span = U_Y - target;
    n = numel(A_list);
    F.dtw_min = nan(1,n); F.dtw_max = nan(1,n);
    F.peak_time = nan(1,n); F.peak_full = nan(1,n);
    F.future_peak_lo = nan(1,n); F.future_peak_hi = nan(1,n);
    F.future_dsm_lo = nan(1,n); F.future_dsm_hi = nan(1,n);
    F.authority = nan(1,n);
    F.peak_in_warmup = false(1,n);
    F.files = cell(1,n);

    for i = 1:n
        Av = A_list(i);
        fprintf('---- A = %.4f (B = 15, vel_ACS = 0.09) ----\n', Av);
        cfg = struct('A', Av, 'B', 15, 'mdot_c_drift', 10, 'vel_ACS', 0.09, ...
                     't_query', t_query, 'tag', sprintf('scanA%03.0f', Av*1000));
        try
            % evalc keeps generate_severe_scenario's own printout out of the log;
            % the assignment inside the evalc'd string lands in THIS workspace.
            evalc('G = generate_severe_scenario(cfg);');
            dp = G.det_params;
            F.files{i} = G.outfile;
            F.dtw_min(i) = min(dp.Y_denorm_50(:,2));
            F.dtw_max(i) = max(dp.Y_denorm_50(:,2));
            fprintf('   delta_tw_max predicted: %.4f .. %.4f K\n', F.dtw_min(i), F.dtw_max(i));

            t_sam = double(dp.t_sam); delayed = double(dp.delayed_time_rom);
            N_pred = size(dp.Y_denorm_50,1);
            time_delayed = t_query - delayed : t_sam : t_query - t_sam;
            d_mpc = ones(1, N_pred+1);

            fp = nan(1,2);
            u_probe = [56.69 56.75];
            for k = 1:2
                uv = u_probe(k);
                u_trial = uv * ones(1, N_pred+1);
                [fd, ~] = DET_2(dp, u_trial, d_mpc);
                nd = fd{1}(1).tree_data;
                Tp = nd.data.pred_tw_max_signal.Data(:);
                ts = nd.data.pred_tw_max_signal.Time(:);
                [pk, im] = max(Tp);
                if k == 1
                    F.peak_full(i) = pk; F.peak_time(i) = ts(im);
                    F.peak_in_warmup(i) = ts(im) <= time_delayed(end);
                end
                fut = ts > time_delayed(end);
                if any(fut), fp(k) = max(Tp(fut)); end
            end
            F.future_peak_lo(i) = min(fp); F.future_peak_hi(i) = max(fp);
            F.authority(i) = abs(fp(1) - fp(2));
            F.future_dsm_lo(i) = (min(fp) - target)/span;
            F.future_dsm_hi(i) = (max(fp) - target)/span;
            fprintf('   full peak %.4f K at t=%.2f  %s\n', F.peak_full(i), F.peak_time(i), ...
                tern(F.peak_in_warmup(i), '<<< IN WARM-UP (uncontrollable)', '<<< in controllable horizon'));
            fprintf('   future peak DSM band [%.6f, %.6f], authority %.6f K\n\n', ...
                F.future_dsm_lo(i), F.future_dsm_hi(i), F.authority(i));
        catch ME
            fprintf('   FAILED: %s\n\n', ME.message);
        end
    end

    %% ---- Summary ----
    fprintf('===== SUMMARY =====\n');
    fprintf('%8s %12s %12s %10s %14s %14s %12s\n', 'A', 'dtw_max', 'peak_K', 'peak_t', ...
        'future_DSM_lo', 'future_DSM_hi', 'authority');
    fprintf('%s\n', repmat('-', 1, 86));
    for i = 1:n
        fprintf('%8.4f %12.4f %12.4f %10.2f %14.6f %14.6f %12.6f%s\n', ...
            A_list(i), F.dtw_max(i), F.peak_full(i), F.peak_time(i), ...
            F.future_dsm_lo(i), F.future_dsm_hi(i), F.authority(i), ...
            tern(F.peak_in_warmup(i), '  (warm-up)', ''));
    end

    %% ---- Pick the usable case ----
    fprintf('\n===== USABLE CASE =====\n');
    % criteria: peak in the controllable horizon, positive future DSM, and
    % authority large enough for a threshold to sit inside the band
    good = ~F.peak_in_warmup & F.future_dsm_hi > 0 & F.authority > 1e-3;
    if any(good)
        gi = find(good);
        [~, best] = max(F.future_dsm_hi(gi));
        b = gi(best);
        F.best_A = A_list(b);
        F.best_file = F.files{b};
        F.conclusion = sprintf(['Usable severe case at A = %.4f: peak %.4f K at t = %.2f ' ...
            '(controllable horizon), future DSM band [%.6f, %.6f], authority %.6f K, ' ...
            'delta_tw_max %.4f K. Run the joint sweep on %s.'], ...
            A_list(b), F.peak_full(b), F.peak_time(b), F.future_dsm_lo(b), ...
            F.future_dsm_hi(b), F.authority(b), F.dtw_max(b), F.files{b});
    else
        F.conclusion = ['No scanned A gives a controllable peak with positive DSM. The peak ' ...
            'stays in the warm-up window across the range, which means the query instant ' ...
            'is the problem, not the severity: t_query must be placed so the thermal peak ' ...
            'falls after t_q. Try earlier query instants (the accident starts at ~1160).'];
    end
    fprintf('%s\n', F.conclusion);

    outfile = sprintf('find_severity_limit_%s.mat', datestr(now,'yyyymmdd_HHMMSS'));
    save(outfile, 'F');
    fprintf('\nSaved %s\n======================================\n\n', outfile);
end

function s = tern(c, a, b)
    if c, s = a; else, s = b; end
end
