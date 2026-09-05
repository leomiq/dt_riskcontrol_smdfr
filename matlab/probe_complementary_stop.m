function P = probe_complementary_stop(matfile)
% PROBE_COMPLEMENTARY_STOP  Identify WHAT stops the DET complementary branch
% one t_sam after the suture instant (t = 1170.5 instead of tsim = 1192.5).
%
% Node 2 of the branch-relevance audit (compl_branch_event = 1,
% initial_branch_state = [1 0]) halts immediately and contributes 0 samples.
% This probe reads SimulationMetadata.ExecutionInfo for BOTH children so the
% difference localizes the trigger, then inspects the model for the blocks that
% consume complementary_branch_flag / initial_branch_state.
%
% Cost: two ROM sims (~26 s + ~11 s) plus static model queries.

    if nargin < 1 || isempty(matfile)
        matfile = 'LJ_test_sampled_vel_ACS_2_loop_1180.mat';
    end

    P = struct();
    P.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fprintf('\n======== Complementary-branch stop probe (t=1170 s) ========\n');
    fprintf('Started %s\n\n', P.timestamp);

    S  = load(matfile);
    fprintf('Regenerating operating points...\n');
    dp = refresh_operating_points(S);
    if ~isfield(dp, 'Y_denorm_50'), dp.Y_denorm_50 = dp.Y_denorm; end

    t_root  = double(dp.time_snapshot_root);
    t_sam   = double(dp.t_sam);
    delayed = double(dp.delayed_time_rom);
    tsim    = double(dp.tsim);
    N_pred  = size(dp.Y_denorm_50, 1);
    name_rom = 'rom_MCDET_v2_2022b_control_CMEX';
    name_wks = dp.name_Wks;

    if ~bdIsLoaded(name_rom), load_system([name_rom '.slx']); end

    %% ---- Static model inspection: who consumes the branch variables? ----
    fprintf('\n-- Blocks referencing the branch control variables --\n');
    vars = {'complementary_branch_flag', 'initial_branch_state', 'parent_last_state', 'parent_last_time'};
    P.var_refs = struct();
    for v = 1:numel(vars)
        vn = vars{v};
        hits = {};
        try
            f = Simulink.findBlocks(name_rom, ...
                'MatchFilter', @(b) contains_var(b, vn));
            for i = 1:numel(f)
                hits{end+1} = getfullname(f(i));
            end
        catch
        end
        if isempty(hits)
            try
                allb = find_system(name_rom, 'LookUnderMasks','all', ...
                    'FollowLinks','on', 'Type','block');
                for i = 1:numel(allb)
                    txt = '';
                    for prm = {'Value','Gain','Expression','InitialCondition','Criteria','Threshold','sfunction_parameters','Parameters'}
                        try
                            q = get_param(allb{i}, prm{1});
                            if ischar(q) || isstring(q), txt = [txt ' ' char(q)]; end
                        catch
                        end
                    end
                    if contains(txt, vn), hits{end+1} = allb{i}; end
                end
            catch ME
                fprintf('  (block scan failed: %s)\n', ME.message);
            end
        end
        P.var_refs.(vn) = hits;
        fprintf('  %-28s -> %d block(s)\n', vn, numel(hits));
        for i = 1:min(numel(hits), 12)
            fprintf('        %s\n', hits{i});
        end
    end

    %% ---- Stop blocks and assertions in the model ----
    fprintf('\n-- Stop / assertion / S-function blocks --\n');
    P.stop_blocks = {};
    P.sfun_blocks = {};
    try
        sb = find_system(name_rom, 'LookUnderMasks','all', 'FollowLinks','on', ...
            'BlockType','Stop');
        P.stop_blocks = sb;
        fprintf('  Stop blocks: %d\n', numel(sb));
        for i = 1:numel(sb), fprintf('     %s\n', sb{i}); end
    catch ME
        fprintf('  Stop scan failed: %s\n', ME.message);
    end
    try
        sf = find_system(name_rom, 'LookUnderMasks','all', 'FollowLinks','on', ...
            'BlockType','S-Function');
        P.sfun_blocks = sf;
        fprintf('  S-Function blocks: %d\n', numel(sf));
        for i = 1:numel(sf)
            fn = '';
            try, fn = get_param(sf{i}, 'FunctionName'); catch, end
            fprintf('     %-58s fcn=%s\n', sf{i}, fn);
        end
    catch ME
        fprintf('  S-Function scan failed: %s\n', ME.message);
    end

    %% ---- Shared inport payload ----
    time_delayed = t_root - delayed : t_sam : t_root - t_sam;
    u_trial = 56.72 * ones(1, N_pred + 1);
    d_mpc   = ones(1, N_pred + 1);
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

    %% ---- Run both children, capturing ExecutionInfo ----
    cases = struct( ...
        'label', {'node1_compl0', 'node2_compl1'}, ...
        'compl', {0, 1}, ...
        'init_state', {[1 1], [1 0]});

    P.cases = struct('label', {}, 'compl', {}, 'init_state', {}, ...
                     'stopEvent', {}, 'stopEventSource', {}, 'stopEventDescription', {}, ...
                     'snapshotTime', {}, 'sim_s', {}, 'tEnd', {}, ...
                     'diagnostics', {}, 'note', {});

    for c = 1:numel(cases)
        cc = cases(c);
        fprintf('\n---- %s (compl=%d, initial_branch_state=%s) ----\n', ...
            cc.label, cc.compl, mat2str(cc.init_state));
        rec = struct('label', cc.label, 'compl', cc.compl, 'init_state', cc.init_state, ...
                     'stopEvent', '', 'stopEventSource', '', 'stopEventDescription', '', ...
                     'snapshotTime', NaN, 'sim_s', NaN, 'tEnd', NaN, ...
                     'diagnostics', {{}}, 'note', '');

        simIn2 = Simulink.SimulationInput(name_rom);
        simIn2 = setVariable(simIn2, 'parent_last_state', [1 1], 'Workspace', name_rom);
        simIn2 = setVariable(simIn2, 'initial_branch_state', cc.init_state, 'Workspace', name_rom);
        simIn2 = setVariable(simIn2, 'complementary_branch_flag', cc.compl, 'Workspace', name_rom);
        simIn2 = setVariable(simIn2, 'parent_last_time', safe_snaptime(dp.romOP), 'Workspace', name_rom);
        simIn2 = setModelParameter(simIn2, 'StopTime', string(tsim));

        op = dp.romOP;
        try, op.loggedStates{13}.Values.Data = 0; catch, end
        try
            simIn2 = setInitialState(simIn2, op);
        catch ME
            rec.note = ['setInitialState: ' ME.message];
            fprintf('   !! %s\n', rec.note);
            P.cases(end+1) = rec; continue;
        end
        simIn2 = simIn2.setExternalInput(inport_data);
        simIn2 = setVariable(simIn2, 'vel_ACS', double(dp.sampled_vel_ACS), 'Workspace', name_wks);
        simIn2 = setModelParameter(simIn2, 'CaptureErrors', 'on');

        tt = tic;
        try
            out2 = sim(simIn2);
            rec.sim_s = toc(tt);
            rec.snapshotTime = safe_snaptime_out(out2);
            fprintf('   sim returned in %.1f s, snapshotTime = %s\n', ...
                rec.sim_s, num2str(rec.snapshotTime));

            % THE KEY READ: what stopped it?
            try
                ei = out2.SimulationMetadata.ExecutionInfo;
                if isfield(ei, 'StopEvent') || isprop(ei, 'StopEvent')
                    rec.stopEvent = char(string(ei.StopEvent));
                end
                try, rec.stopEventSource = char(string(ei.StopEventSource)); catch, end
                try, rec.stopEventDescription = char(string(ei.StopEventDescription)); catch, end
                fprintf('   >> StopEvent            = %s\n', rec.stopEvent);
                fprintf('   >> StopEventSource      = %s\n', rec.stopEventSource);
                fprintf('   >> StopEventDescription = %s\n', rec.stopEventDescription);

                % Diagnostics (warnings/errors raised during the run)
                try
                    dg = ei.ErrorDiagnostic;
                    if ~isempty(dg)
                        rec.diagnostics{end+1} = char(string(dg.Message));
                        fprintf('   >> ErrorDiagnostic: %s\n', char(string(dg.Message)));
                    end
                catch
                end
                try
                    wr = ei.WarningDiagnostics;
                    for w = 1:numel(wr)
                        m = char(string(wr(w).Diagnostic.Message));
                        rec.diagnostics{end+1} = m;
                        fprintf('   >> Warning %d: %s\n', w, m);
                    end
                catch
                end
            catch ME
                fprintf('   (ExecutionInfo unreadable: %s)\n', ME.message);
            end

            % Where did time actually end?
            try
                rec.tEnd = out2.tout(end);
                fprintf('   tout(end) = %.4f  (StopTime was %.1f)\n', rec.tEnd, tsim);
            catch
                try
                    rec.tEnd = out2.pred_tw_max_signal.Time(end);
                    fprintf('   pred_tw_max Time(end) = %.4f  (StopTime was %.1f)\n', rec.tEnd, tsim);
                catch
                end
            end
        catch ME
            rec.sim_s = toc(tt);
            rec.note = sprintf('[%s] %s', ME.identifier, ME.message);
            fprintf('   sim FAILED after %.1f s: %s\n', rec.sim_s, rec.note);
        end
        P.cases(end+1) = rec;
    end

    %% ---- Verdict ----
    fprintf('\n======== VERDICT ========\n');
    for c = 1:numel(P.cases)
        r = P.cases(c);
        fprintf('%s: tEnd=%s  StopEvent=%s  Source=%s\n', ...
            r.label, num2str(r.tEnd), r.stopEvent, r.stopEventSource);
    end
    P.conclusion = '';
    if numel(P.cases) == 2
        a = P.cases(1); b = P.cases(2);
        if ~isnan(b.tEnd) && ~isnan(a.tEnd) && b.tEnd < a.tEnd
            P.conclusion = sprintf(['Complementary branch ends at %.4f vs %.4f for the ' ...
                'transition branch. StopEvent="%s", Source="%s" names the trigger.'], ...
                b.tEnd, a.tEnd, b.stopEvent, b.stopEventSource);
        else
            P.conclusion = 'Both branches reached the same end time in this probe.';
        end
    end
    fprintf('\n%s\n', P.conclusion);

    outfile = sprintf('probe_complementary_stop_%s.mat', datestr(now,'yyyymmdd_HHMMSS'));
    save(outfile, 'P');
    fprintf('\nSaved %s\n=========================\n\n', outfile);
end

function tf = contains_var(b, vn)
    tf = false;
    try
        for prm = {'Value','Gain','Expression','InitialCondition','Criteria','Threshold'}
            q = get_param(b, prm{1});
            if (ischar(q) || isstring(q)) && contains(char(q), vn), tf = true; return; end
        end
    catch
    end
end

function t = safe_snaptime(op)
    t = NaN;
    try, t = op.snapshotTime; catch, end
end

function t = safe_snaptime_out(out)
    t = NaN;
    try, t = out.finalOP.snapshotTime; catch, end
end
