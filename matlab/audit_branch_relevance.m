function A = audit_branch_relevance(matfile)
% AUDIT_BRANCH_RELEVANCE  Decide whether the commented-out sim at DET_2.m:450
% is DELIBERATE PRUNING of decision-irrelevant branches, or an accidental defect.
%
% Leo's argument: branches where the ACS fails are not worth visiting, because
% they cannot change the optimal control decision. If the skipped child nodes
% are exactly those, the commented sim is a design choice. If productive
% (decision-relevant) nodes are also being skipped, it is a defect.
%
% This walks the child-branch expansion at the suture instant t = 1170 s WITHOUT
% the optimizer (a handful of ROM sims at ~10 s each), classifying every node the
% queue would generate:
%
%   ACS-FAILED / TERMINAL : parent_last_state(1)==1 && new_queue_state(1)==0
%                           or state_var_ps == 2        (DET_2.m:474-479)
%   PRODUCTIVE            : anything else -> the control action still matters
%
% Key structural point being tested: the termination tests at DET_2.m:474-479
% run AFTER line 462 and depend on finalOP, which requires the sim. So a
% commented-out sim cannot prune selectively — it kills the whole child branch.
% This script measures whether that distinction has practical consequences.

    if nargin < 1 || isempty(matfile)
        matfile = 'LJ_test_sampled_vel_ACS_2_loop_1180.mat';
    end

    A = struct();
    A.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fprintf('\n=========== DET branch-relevance audit (t=1170 s) ===========\n');
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
    N_branches = double(dp.N_branches);
    name_rom = 'rom_MCDET_v2_2022b_control_CMEX';
    name_wks = dp.name_Wks;

    A.t_branch = t_root - delayed;
    A.N_branches = N_branches;
    fprintf('Suture/branch instant = %.1f s, N_branches = %d, tsim = %.1f s\n\n', ...
        A.t_branch, N_branches, tsim);

    if ~bdIsLoaded(name_rom), load_system([name_rom '.slx']); end

    time_delayed = t_root - delayed : t_sam : t_root - t_sam;
    u_fixed = 56.72;
    u_trial = u_fixed * ones(1, N_pred + 1);
    d_mpc   = ones(1, N_pred + 1);

    % Shared inport payload (identical for every branch at this instant).
    outRoot = dp.outRoot; outrom = dp.outrom; Y50 = dp.Y_denorm_50;
    pred_mdot = cat(2, outRoot.mdot_c_leak_signal.Data(end,1), Y50(:,1)');
    time_ext  = cat(2, time_delayed, t_root : t_sam : tsim);
    pred_mdot_ext = cat(2, ...
        outRoot.mdot_c_leak_signal.Data(end-numel(time_delayed):end-1,1)', pred_mdot);
    pred_dtw_delayed = outRoot.tw_max_signal.Data(end-numel(time_delayed):end,1) - ...
                       outrom.tw_max_signal.Data(end-numel(time_delayed):end,1);
    pred_dtw_ext = cat(2, pred_dtw_delayed', Y50(:,2)');
    cs = dp.control_series;
    [ci, ~] = find(cs(:,1) == time_delayed);
    u_ext = cat(2, cs(ci,2)', u_trial);
    d_ext = cat(2, cs(ci,3)', d_mpc);
    inport_data = [time_ext', pred_mdot_ext', pred_dtw_ext', u_ext', d_ext'];

    %% ---- Expand the first generation of child nodes ----
    % Parent = the root at the suture. Its two children differ by
    % compl_branch_event_value (0 = follow the transition, 1 = complement).
    parent_finalOP = dp.romOP;
    parent_last_state = [1 1];

    A.nodes = struct('id', {}, 'compl', {}, 'init_state', {}, ...
                     'sim_ok', {}, 'sim_s', {}, 'snapshotTime', {}, ...
                     'new_queue_state', {}, 'state_var_ps', {}, ...
                     'class', {}, 'peak_Tp', {}, 'segment_len', {}, 'note', {});

    for k = 1:N_branches
        compl = double(mod(k,2) == 0);   % DET_2.m:490-497 -> i even => 1
        fprintf('---- child node %d (compl_branch_event = %d) ----\n', k, compl);

        n = struct('id', k, 'compl', compl, 'sim_ok', false, 'sim_s', NaN, ...
                   'snapshotTime', NaN, 'new_queue_state', [], 'state_var_ps', [], ...
                   'class', 'UNKNOWN', 'peak_Tp', NaN, 'segment_len', NaN, 'note', '');

        simIn2 = Simulink.SimulationInput(name_rom);
        simIn2 = setVariable(simIn2, 'parent_last_state', parent_last_state, 'Workspace', name_rom);

        % Branch-state logic verbatim from DET_2.m:315-338.
        if compl == 0
            init_state = parent_last_state;
        else
            if isequal(parent_last_state, [1 0])
                init_state = [1 1];
            elseif isequal(parent_last_state, [1 1])
                init_state = [1 0];
            else
                init_state = parent_last_state;   % terminal: [0 1] / [0 0]
                compl = 0;
            end
        end
        n.init_state = init_state;
        fprintf('   initial_branch_state = %s\n', mat2str(init_state));

        simIn2 = setVariable(simIn2, 'initial_branch_state', init_state, 'Workspace', name_rom);
        simIn2 = setVariable(simIn2, 'complementary_branch_flag', compl, 'Workspace', name_rom);
        simIn2 = setVariable(simIn2, 'parent_last_time', safe_snaptime(parent_finalOP), ...
                             'Workspace', name_rom);
        simIn2 = setModelParameter(simIn2, 'StopTime', string(tsim));

        op = parent_finalOP;
        try, op.loggedStates{13}.Values.Data = 0; catch, end
        try
            simIn2 = setInitialState(simIn2, op);
        catch ME
            n.note = ['setInitialState failed: ' ME.message];
            fprintf('   !! %s\n', n.note);
            A.nodes(end+1) = n; continue;
        end

        simIn2 = simIn2.setExternalInput(inport_data);
        simIn2 = setVariable(simIn2, 'vel_ACS', double(dp.sampled_vel_ACS), 'Workspace', name_wks);

        tt = tic;
        try
            out2 = sim(simIn2);
            n.sim_ok = true;
            n.sim_s = toc(tt);
            fOP = out2.finalOP;
            n.snapshotTime = safe_snaptime(fOP);
            try, n.new_queue_state = fOP.loggedStates{12}.Values.Data; catch, end
            try, n.state_var_ps    = fOP.loggedStates{11}.Values.Data; catch, end
            fprintf('   sim OK (%.1f s), snapshotTime = %s\n', n.sim_s, num2str(n.snapshotTime));
            fprintf('   new_queue_state = %s, state_var_ps = %s\n', ...
                mat2str(n.new_queue_state), mat2str(n.state_var_ps));

            % Classification per the termination tests at DET_2.m:474-479.
            nq = n.new_queue_state;
            if ~isempty(nq) && parent_last_state(1) == 1 && nq(1) == 0
                n.class = 'ACS-FAILED/TERMINAL (demand dropped)';
            elseif ~isempty(n.state_var_ps) && all(n.state_var_ps == 2)
                n.class = 'ACS-FAILED/TERMINAL (state_var_ps == 2)';
            else
                n.class = 'PRODUCTIVE (decision-relevant)';
            end
            fprintf('   CLASS: %s\n', n.class);

            % Contribution to the reconstructed trajectory.
            try
                tsig = out2.pred_tw_max_signal.Time(:);
                psig = out2.pred_tw_max_signal.Data(:);
                idx  = find(tsig == time_delayed(end));
                if ~isempty(idx) && idx(end) < numel(tsig)
                    n.segment_len = numel(tsig) - idx(end);
                    n.peak_Tp = max(psig(idx(end)+1:end));
                else
                    n.segment_len = 0;
                end
                fprintf('   segment = %d samples, peak Tp = %.4f K\n', ...
                    n.segment_len, n.peak_Tp);
            catch ME
                n.note = ['signal read: ' ME.message];
            end
        catch ME
            n.sim_s = toc(tt);
            n.note = sprintf('[%s] %s', ME.identifier, ME.message);
            fprintf('   sim FAILED after %.1f s: %s\n', n.sim_s, n.note);
        end
        A.nodes(end+1) = n;
        fprintf('\n');
    end

    %% ---- Verdict ----
    fprintf('=========== VERDICT ===========\n');
    cls = {A.nodes.class};
    prod_mask = contains(cls, 'PRODUCTIVE');
    A.n_total = numel(A.nodes);
    A.n_productive = sum(prod_mask);
    A.n_terminal = sum(contains(cls, 'ACS-FAILED'));
    fprintf('child nodes expanded : %d\n', A.n_total);
    fprintf('  PRODUCTIVE         : %d\n', A.n_productive);
    fprintf('  ACS-FAILED/TERMINAL: %d\n', A.n_terminal);

    peaks = [A.nodes.peak_Tp];
    A.peak_spread = NaN;
    if sum(~isnan(peaks)) > 1
        A.peak_spread = max(peaks(~isnan(peaks))) - min(peaks(~isnan(peaks)));
        fprintf('peak Tp spread across children: %.4f K\n', A.peak_spread);
    end

    if A.n_productive > 0
        A.conclusion = sprintf(['DEFECT, not pruning: %d of %d child nodes are ' ...
            'decision-relevant (PRODUCTIVE) yet the commented sim at DET_2.m:450 ' ...
            'skips them all. Deliberate ACS pruning would keep these.'], ...
            A.n_productive, A.n_total);
    elseif A.n_total > 0
        A.conclusion = ['CONSISTENT WITH PRUNING: every expanded child is ' ...
            'ACS-failed/terminal, so skipping them cannot change the optimal action.'];
    else
        A.conclusion = 'INCONCLUSIVE: no child nodes could be expanded.';
    end
    fprintf('\n%s\n', A.conclusion);

    outfile = sprintf('audit_branch_relevance_%s.mat', datestr(now,'yyyymmdd_HHMMSS'));
    save(outfile, 'A');
    fprintf('\nSaved %s\n', outfile);
    fprintf('===============================\n\n');
end

function t = safe_snaptime(op)
    t = NaN;
    try, t = op.snapshotTime; catch, end
end
