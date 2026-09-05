function E = probe_root_expansion(matfile)
% PROBE_ROOT_EXPANSION  Why does the DET tree stay at a single node even after
% Fix 1 restored the child simulation?
%
% Hypothesis: the root branch runs to StopTime (the probe measured
% snapshotTime = 1192.5 == tsim), so the guard at DET_2.m:243
%
%     if time_snapshot < tsim        % <- false when the root reaches tsim
%         ... expand N_resampling children ...
%     end
%
% never fires and no child nodes are ever queued. In that reading the DET
% degenerates to one node for a structural reason, independent of the
% commented-out sim that Fix 1 addressed.
%
% This reproduces the ROOT branch exactly as DET_2.m:85-233 does and reports
% every quantity the guard depends on, without running the optimizer.

    if nargin < 1 || isempty(matfile)
        matfile = 'LJ_test_sampled_vel_ACS_2_loop_1180.mat';
    end

    E = struct();
    E.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fprintf('\n======== Root-expansion guard probe ========\n');
    fprintf('Started %s\n\n', E.timestamp);

    S  = load(matfile);
    fprintf('Regenerating operating points...\n');
    dp = refresh_operating_points(S);
    if ~isfield(dp, 'Y_denorm_50'), dp.Y_denorm_50 = dp.Y_denorm; end

    t_root  = double(dp.time_snapshot_root);
    t_sam   = double(dp.t_sam);
    delayed = double(dp.delayed_time_rom);
    tsim    = double(dp.tsim);
    N_pred  = size(dp.Y_denorm_50, 1);
    N_resampling = double(dp.N_resampling);
    N_branches   = double(dp.N_branches);
    name_rom = 'rom_MCDET_v2_2022b_control_CMEX';
    name_wks = dp.name_Wks;

    E.tsim = tsim; E.N_resampling = N_resampling; E.N_branches = N_branches;
    fprintf('tsim = %.2f, N_resampling = %d, N_branches = %d\n', tsim, N_resampling, N_branches);
    fprintf('t_root = %.1f, suture = %.1f, t_sam = %.2f\n\n', t_root, t_root - delayed, t_sam);

    if ~bdIsLoaded(name_rom), load_system([name_rom '.slx']); end

    %% ---- Reproduce the ROOT branch (DET_2.m:85-220) ----
    Raiz_to_Rom = [1,2,4,5,6,7,3,8,9,10,22,23,24];
    romOP_temp = dp.romOP;
    finalOP2_Raiz = dp.finalOP2_Raiz;

    fprintf('-- Feeding the ROM operating point from the HRM (DET_2.m:90-101) --\n');
    E.transfer_ok = true;
    E.transfer_errors = {};
    for index = [1:2, 6:13]
        try
            j = Raiz_to_Rom(index);
            vec = finalOP2_Raiz.loggedStates{j}.Values.Data;
            romOP_temp.loggedStates{index}.Values.Data = vec;
        catch ME
            E.transfer_ok = false;
            E.transfer_errors{end+1} = sprintf('index %d -> %d: %s', index, Raiz_to_Rom(index), ME.message);
        end
    end
    for index = [3,4,5]
        try
            j = Raiz_to_Rom(index);
            vec = finalOP2_Raiz.loggedStates{j}.Values.Data;
            romOP_temp.loggedStates{index}.Values.Data = mean(vec(1:12));
        catch ME
            E.transfer_ok = false;
            E.transfer_errors{end+1} = sprintf('index %d -> %d: %s', index, Raiz_to_Rom(index), ME.message);
        end
    end
    fprintf('   HRM->ROM state transfer clean: %d\n', E.transfer_ok);
    for k = 1:numel(E.transfer_errors)
        fprintf('   !! %s\n', E.transfer_errors{k});
    end
    fprintf('   (this is the Raiz_to_Rom mapping; loggedStates has %d entries)\n\n', ...
        numel(dp.romOP.loggedStates));

    %% ---- Inport payload, as the root branch builds it ----
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

    %% ---- Run the root exactly as DET_2 does ----
    fprintf('-- Running the ROOT branch --\n');
    simIn = Simulink.SimulationInput(name_rom);
    simIn = setInitialState(simIn, romOP_temp);
    simIn = setModelParameter(simIn, 'StopTime', string(tsim));
    simIn = setVariable(simIn, 'initial_branch_state', ...
        finalOP2_Raiz.loggedStates{23}.Values.Data, 'Workspace', name_rom);
    simIn = simIn.setExternalInput(inport_data);
    simIn = setVariable(simIn, 'vel_ACS', double(dp.sampled_vel_ACS), 'Workspace', name_wks);

    tt = tic;
    try
        out = sim(simIn);
        E.root_sim_ok = true;
        E.root_sim_s = toc(tt);
        fprintf('   root sim OK in %.1f s\n', E.root_sim_s);
    catch ME
        E.root_sim_ok = false;
        E.root_error = sprintf('[%s] %s', ME.identifier, ME.message);
        fprintf('   root sim FAILED: %s\n', E.root_error);
        finish(E); return;
    end

    try
        ei = out.SimulationMetadata.ExecutionInfo;
        E.root_stopEvent = char(string(ei.StopEvent));
        try, E.root_stopDesc = char(string(ei.StopEventDescription)); catch, E.root_stopDesc = ''; end
        fprintf('   StopEvent = %s  (%s)\n', E.root_stopEvent, E.root_stopDesc);
    catch
    end

    %% ---- THE GUARD (DET_2.m:236-251) ----
    fprintf('\n-- Evaluating the expansion guard at DET_2.m:243 --\n');
    finalOP = out.finalOP;
    E.time_snapshot = safe(@() finalOP.snapshotTime);
    fprintf('   time_snapshot = finalOP.snapshotTime = %s\n', num2str(E.time_snapshot));
    fprintf('   tsim                                 = %.4f\n', tsim);
    E.guard_passes = ~isnan(E.time_snapshot) && E.time_snapshot < tsim;
    fprintf('   >> (time_snapshot < tsim) = %d\n', E.guard_passes);

    E.parent_last_state = safe(@() finalOP2_Raiz.loggedStates{23}.Values.Data);
    E.new_queue_state   = safe(@() finalOP.loggedStates{12}.Values.Data);
    E.state_var_ps      = safe(@() finalOP.loggedStates{11}.Values.Data);
    fprintf('   parent_last_state (finalOP2_Raiz{23}) = %s\n', mat2str(E.parent_last_state));
    fprintf('   new_queue_state   (finalOP{12})       = %s\n', mat2str(E.new_queue_state));
    fprintf('   state_var_ps      (finalOP{11})       = %s\n', mat2str(E.state_var_ps));

    if ~E.guard_passes
        fprintf('\n   >> GUARD BLOCKS EXPANSION: the root reached tsim, so DET_2 never\n');
        fprintf('      queues any child. The tree stays at 1 node regardless of Fix 1.\n');
        E.conclusion = ['Root reaches tsim, so the guard `time_snapshot < tsim` is false ' ...
            'and no children are ever queued. Single-node tree is caused by the guard, ' ...
            'not by the commented-out child sim.'];
    else
        % Which terminal test would fire?
        pls = E.parent_last_state; nqs = E.new_queue_state; svp = E.state_var_ps;
        if ~isempty(pls) && ~isempty(nqs) && pls(1) == 1 && nqs(1) == 0
            E.conclusion = 'Guard passes but terminal test 1 fires (demand dropped): no children.';
        elseif ~isempty(svp) && all(svp == 2)
            E.conclusion = 'Guard passes but terminal test 2 fires (state_var_ps == 2): no children.';
        else
            E.conclusion = sprintf('Guard passes and no terminal test fires: %d children should be queued.', N_resampling);
        end
        fprintf('\n   >> %s\n', E.conclusion);
    end

    fprintf('\n======== VERDICT ========\n%s\n', E.conclusion);
    finish(E);
end

function v = safe(f)
    v = NaN;
    try, v = f(); catch, end
end

function finish(E)
    outfile = sprintf('probe_root_expansion_%s.mat', datestr(now,'yyyymmdd_HHMMSS'));
    save(outfile, 'E');
    fprintf('\nSaved %s\n=========================\n\n', outfile);
end
