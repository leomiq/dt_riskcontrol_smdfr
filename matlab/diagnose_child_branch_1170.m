function R = diagnose_child_branch_1170(matfile)
% DIAGNOSE_CHILD_BRANCH_1170  Single-shot diagnosis of the DET child-node
% propagation failure at the HRM->ROM suture instant t = 1170 s.
%
%   R = diagnose_child_branch_1170()            % uses the archived 1180 workspace
%   R = diagnose_child_branch_1170('other.mat')
%
% WHY THIS SCRIPT: a full closed-loop realization costs ~6 h (~350 s per query
% x 61 queries). This runs ONE child-branch simulation with a fixed u_trial and
% NO optimizer in the loop (1 objective evaluation instead of 420), so the
% failure at the known-bad instant is reproduced in minutes.
%
% HYPOTHESES UNDER TEST (in order of cost to check):
%   H1  out2 = sim(simIn2) is commented out at DET_2.m:450 while out2.finalOP is
%       consumed at 462 -> the child's tree_data entry can never be written.
%   H2  parent_finalOP.snapshotTime == time_delayed(end) == 1170, so
%       find(t_signal == time_delayed(end)) in fitness_cost.m (lines 39-41,
%       58-59) returns an empty/off-by-one index and the child segment is
%       discarded during reconstruction.
%   H3  FastRestart='on' + SimulationMode='accelerator' + SaveOperatingPoint
%       (DET_2.m:8-13), combined with mutating
%       parent_finalOP.loggedStates{13}.Values.Data = 0 (line 350), invalidates
%       the operating-point checksum and setInitialState is rejected.
%
% The script does NOT modify DET_2.m. It reproduces the child-branch setup
% verbatim, runs it inside try/catch, and reports which hypothesis fires.

    % The archived workspace is the only source that carries det_params with the
    % OperatingPoint objects. phase2c_fixed_state_results.mat holds only summary
    % structs (seeds/config/LJ/GA/LJstats/GAstats), no det_params.
    if nargin < 1 || isempty(matfile)
        matfile = 'LJ_test_sampled_vel_ACS_2_loop_1180.mat';
    end

    R = struct();
    R.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    R.matfile   = matfile;
    fprintf('\n================ DET child-branch diagnosis (t=1170 s) ================\n');
    fprintf('Started %s\n', R.timestamp);

    %% ---- 0. Load the cached det_params (operating points already regenerated) ----
    S = load(matfile);

    % The archived Oct-2025 OperatingPoint objects refuse to load into R2024a
    % ("model was changed after the operating point was saved"), so rebuild them
    % exactly as the experiment harness does. This is the ~2 min cost of the run.
    fprintf('\nRegenerating operating points (archived OPs do not load in R2024a)...\n');
    dp = refresh_operating_points(S);
    if ~isfield(dp, 'Y_denorm_50')
        dp.Y_denorm_50 = dp.Y_denorm;   % 11_10 series stored it as Y_denorm
    end
    R.det_params_source = matfile;

    t_root   = double(dp.time_snapshot_root);
    t_sam    = double(dp.t_sam);
    delayed  = double(dp.delayed_time_rom);
    t_branch = t_root - delayed;              % == 1170 for t_root = 1180

    fprintf('\n-- Design constants --\n');
    fprintf('  time_snapshot_root  = %.1f s\n', t_root);
    fprintf('  delayed_time_rom    = %.1f s\n', delayed);
    fprintf('  t_sam               = %.2f s\n', t_sam);
    fprintf('  SUTURE / BRANCH t   = %.1f s   <-- instant under investigation\n', t_branch);
    R.t_root = t_root; R.delayed = delayed; R.t_sam = t_sam; R.t_branch = t_branch;

    %% ---- 1. H2: snapshotTime vs time_delayed(end) collision ----
    fprintf('\n-- H2: snapshotTime / time_delayed collision --\n');
    time_delayed = t_root - delayed : t_sam : t_root - t_sam;
    R.time_delayed_first = time_delayed(1);
    R.time_delayed_last  = time_delayed(end);
    fprintf('  time_delayed = [%.1f ... %.1f], numel = %d\n', ...
        time_delayed(1), time_delayed(end), numel(time_delayed));

    op2 = dp.finalOP2_Raiz;
    romOP = dp.romOP;
    R.finalOP2_snapshotTime = safe_snaptime(op2);
    R.romOP_snapshotTime    = safe_snaptime(romOP);
    fprintf('  finalOP2_Raiz.snapshotTime = %s\n', num2str(R.finalOP2_snapshotTime));
    fprintf('  romOP.snapshotTime         = %s\n', num2str(R.romOP_snapshotTime));

    % The first child inherits parent_last_time = parent_finalOP.snapshotTime.
    % For the first child the parent is the root, whose OP is romOP at t_branch.
    R.H2_collision = isequal_tol(R.romOP_snapshotTime, time_delayed(end)) || ...
                     isequal_tol(R.romOP_snapshotTime, t_branch);
    if isequal_tol(R.romOP_snapshotTime, t_branch)
        fprintf('  >> romOP.snapshotTime == branch instant %.1f s\n', t_branch);
    end
    if isequal_tol(R.romOP_snapshotTime, time_delayed(end))
        fprintf('  >> H2 FIRES: snapshotTime == time_delayed(end) == %.1f\n', time_delayed(end));
    else
        fprintf('  >> snapshotTime (%s) vs time_delayed(end) (%.1f): offset %s s\n', ...
            num2str(R.romOP_snapshotTime), time_delayed(end), ...
            num2str(R.romOP_snapshotTime - time_delayed(end)));
    end

    %% ---- 2. H2b: does the find() in fitness_cost actually resolve? ----
    fprintf('\n-- H2b: reconstruction index on the ROOT signal --\n');
    try
        t_sig = dp.outrom.tw_max_signal.Time(:);
        idx = find(t_sig == time_delayed(end));
        R.root_find_idx = idx;
        R.root_t_sig_range = [min(t_sig) max(t_sig)];
        fprintf('  outrom tw_max Time range = [%.2f, %.2f], n = %d\n', ...
            min(t_sig), max(t_sig), numel(t_sig));
        fprintf('  find(t_signal == %.1f) -> %s\n', time_delayed(end), mat2str(idx(:)'));
        if isempty(idx)
            fprintf('  >> H2b FIRES: empty index. Tp_signal(idx+1:end) would ERROR/return [].\n');
            R.H2b_empty_index = true;
        else
            R.H2b_empty_index = false;
            if idx(end) >= numel(t_sig)
                fprintf('  >> H2b FIRES: idx at end of signal; Tp_signal(idx+1:end) is EMPTY.\n');
                R.H2b_empty_slice = true;
            else
                R.H2b_empty_slice = false;
                fprintf('  slice length after idx+1 = %d samples\n', numel(t_sig) - idx(end));
            end
        end
    catch ME
        R.H2b_error = ME.message;
        fprintf('  could not evaluate: %s\n', ME.message);
    end

    %% ---- 3. loggedStates map (indices DET_2 accesses positionally) ----
    fprintf('\n-- loggedStates map (DET_2 uses {11},{12},{13},{23} positionally) --\n');
    R.loggedStates = dump_logged_states(romOP, 'romOP');
    R.loggedStates_op2 = dump_logged_states(op2, 'finalOP2_Raiz');

    %% ---- 4. H3: FastRestart / accelerator / OperatingPoint conflict ----
    fprintf('\n-- H3: model config vs external OperatingPoint --\n');
    name_rom = 'rom_MCDET_v2_2022b_control_CMEX';
    if ~bdIsLoaded(name_rom)
        load_system([name_rom '.slx']);
    end
    cfgnames = {'SimulationMode','FastRestart','SaveFinalState','SaveOperatingPoint','FinalStateName'};
    for i = 1:numel(cfgnames)
        try
            v = get_param(name_rom, cfgnames{i});
            R.cfg.(cfgnames{i}) = v;
            fprintf('  %-20s = %s\n', cfgnames{i}, string(v));
        catch
            fprintf('  %-20s = <unavailable>\n', cfgnames{i});
        end
    end

    %% ---- 5. THE experiment: one child-branch simulation ----
    fprintf('\n-- Child-branch simulation (1 evaluation, no optimizer) --\n');
    N_pred = size(dp.Y_denorm_50, 1);
    u_fixed = 56.72;                                   % mid-range, arbitrary but fixed
    u_trial = u_fixed * ones(1, N_pred + 1);
    d_mpc   = ones(1, N_pred + 1);
    R.u_trial_value = u_fixed;
    R.N_pred = N_pred;
    fprintf('  u_trial = %.4f A (constant), N_pred = %d\n', u_fixed, N_pred);

    % Reproduce the child-branch setup verbatim from DET_2.m:302-354.
    parent_finalOP = romOP;
    parent_last_state = [1 1];
    compl_branch_event_value = 0;

    simIn2 = Simulink.SimulationInput(name_rom);
    name_wks = dp.name_Wks;
    tsim = double(dp.tsim);

    simIn2 = setVariable(simIn2, 'parent_last_state', parent_last_state, 'Workspace', name_rom);
    simIn2 = setVariable(simIn2, 'initial_branch_state', parent_last_state, 'Workspace', name_rom);
    simIn2 = setVariable(simIn2, 'complementary_branch_flag', compl_branch_event_value, 'Workspace', name_rom);

    parent_last_time = safe_snaptime(parent_finalOP);
    simIn2 = setVariable(simIn2, 'parent_last_time', parent_last_time, 'Workspace', name_rom);
    simIn2 = setModelParameter(simIn2, 'StopTime', string(tsim));
    fprintf('  parent_last_time = %s, StopTime = %.1f\n', num2str(parent_last_time), tsim);

    % The mutation from DET_2.m:350 that is suspected of invalidating the OP.
    R.step_counter_mutated = false;
    try
        parent_finalOP.loggedStates{13}.Values.Data = 0;
        R.step_counter_mutated = true;
        fprintf('  loggedStates{13} (step_counter) set to 0  [as DET_2.m:350]\n');
    catch ME
        fprintf('  !! could not mutate loggedStates{13}: %s\n', ME.message);
        R.step_counter_error = ME.message;
    end

    % setInitialState — H3 fires here if FastRestart/accelerator rejects the OP.
    R.setInitialState_ok = false;
    try
        simIn2 = setInitialState(simIn2, parent_finalOP);
        R.setInitialState_ok = true;
        fprintf('  setInitialState: OK\n');
    catch ME
        R.setInitialState_error = ME.message;
        R.setInitialState_id = ME.identifier;
        fprintf('  >> H3 FIRES at setInitialState: [%s] %s\n', ME.identifier, ME.message);
    end

    % Inport data, verbatim from DET_2.m:408-438.
    try
        Y50 = dp.Y_denorm_50;
        outRoot = dp.outRoot; outrom = dp.outrom;
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

        R.len_time_ext = numel(time_ext);
        R.len_mdot_ext = numel(pred_mdot_ext);
        R.len_dtw_ext  = numel(pred_dtw_ext);
        R.len_u_ext    = numel(u_ext);
        R.len_d_ext    = numel(d_ext);
        fprintf('  inport lengths: time=%d mdot=%d dtw=%d u=%d d=%d\n', ...
            R.len_time_ext, R.len_mdot_ext, R.len_dtw_ext, R.len_u_ext, R.len_d_ext);
        L = [R.len_time_ext R.len_mdot_ext R.len_dtw_ext R.len_u_ext R.len_d_ext];
        if numel(unique(L)) > 1
            fprintf('  >> INPORT LENGTH MISMATCH: %s  (child branch cannot be built)\n', mat2str(L));
            R.inport_mismatch = true;
        else
            R.inport_mismatch = false;
        end

        inport_data = [time_ext', pred_mdot_ext', pred_dtw_ext', u_ext', d_ext'];
        simIn2 = simIn2.setExternalInput(inport_data);
        simIn2 = setVariable(simIn2, 'vel_ACS', double(dp.sampled_vel_ACS), 'Workspace', name_wks);
        R.inport_ok = true;
    catch ME
        R.inport_ok = false;
        R.inport_error = ME.message;
        fprintf('  >> inport construction FAILED: %s\n', ME.message);
    end

    % H1: run the sim that DET_2.m:450 has commented out.
    fprintf('\n  Running sim(simIn2)  [the call commented out at DET_2.m:450] ...\n');
    R.sim_ok = false;
    tsim_tic = tic;
    try
        out2 = sim(simIn2);
        R.sim_ok = true;
        R.sim_seconds = toc(tsim_tic);
        fprintf('  sim: OK in %.1f s\n', R.sim_seconds);

        % Can the child's tree_data entry actually be written?
        try
            finalOP_child = out2.finalOP;
            R.child_finalOP_ok = true;
            R.child_snapshotTime = safe_snaptime(finalOP_child);
            fprintf('  out2.finalOP: OK, snapshotTime = %s\n', num2str(R.child_snapshotTime));
            try
                R.child_queue_state = finalOP_child.loggedStates{12}.Values.Data;
                R.child_state_var_ps = finalOP_child.loggedStates{11}.Values.Data;
                fprintf('  loggedStates{12} (new_queue_state) = %s\n', mat2str(R.child_queue_state));
                fprintf('  loggedStates{11} (state_var_ps)    = %s\n', mat2str(R.child_state_var_ps));
            catch ME
                R.child_state_error = ME.message;
                fprintf('  >> loggedStates read FAILED: %s\n', ME.message);
            end
        catch ME
            R.child_finalOP_ok = false;
            R.child_finalOP_error = ME.message;
            fprintf('  >> H1 CONFIRMED downstream: out2.finalOP unavailable: %s\n', ME.message);
        end

        % Did the child produce a usable Tp segment past the suture?
        try
            tsig = out2.pred_tw_max_signal.Time(:);
            psig = out2.pred_tw_max_signal.Data(:);
            idxc = find(tsig == time_delayed(end));
            R.child_find_idx = idxc;
            R.child_t_range = [min(tsig) max(tsig)];
            fprintf('  child pred_tw_max Time = [%.2f, %.2f], n = %d\n', ...
                min(tsig), max(tsig), numel(tsig));
            fprintf('  find(t == %.1f) -> %s\n', time_delayed(end), mat2str(idxc(:)'));
            if isempty(idxc)
                fprintf('  >> H2 CONFIRMED: empty index -> child segment discarded.\n');
                R.H2_confirmed = true;
            elseif idxc(end) >= numel(tsig)
                fprintf('  >> H2 CONFIRMED: slice Tp(idx+1:end) is EMPTY.\n');
                R.H2_confirmed = true;
            else
                R.H2_confirmed = false;
                R.child_segment_len = numel(tsig) - idxc(end);
                R.child_peak_Tp = max(psig(idxc(end)+1:end));
                fprintf('  child segment = %d samples, peak Tp = %.4f K\n', ...
                    R.child_segment_len, R.child_peak_Tp);
                fprintf('  (1244 K limit -> margin %.4f K)\n', 1244 - R.child_peak_Tp);
            end
        catch ME
            R.child_signal_error = ME.message;
            fprintf('  >> child pred_tw_max_signal unavailable: %s\n', ME.message);
        end
    catch ME
        R.sim_seconds = toc(tsim_tic);
        R.sim_error = ME.message;
        R.sim_error_id = ME.identifier;
        fprintf('  >> sim FAILED after %.1f s\n', R.sim_seconds);
        fprintf('     identifier: %s\n', ME.identifier);
        fprintf('     message   : %s\n', ME.message);
        for k = 1:numel(ME.stack)
            fprintf('     at %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
        end
        if ~isempty(ME.cause)
            for k = 1:numel(ME.cause)
                fprintf('     cause %d: %s\n', k, ME.cause{k}.message);
            end
        end
    end

    %% ---- 6. Verdict ----
    fprintf('\n-- VERDICT --\n');
    R.verdict = {};
    if isfield(R,'setInitialState_error')
        R.verdict{end+1} = 'H3: setInitialState rejected the external OperatingPoint (FastRestart/accelerator conflict).';
    end
    if isfield(R,'sim_error')
        R.verdict{end+1} = sprintf('sim() failed: %s', R.sim_error_id);
    end
    if isfield(R,'H2_confirmed') && R.H2_confirmed
        R.verdict{end+1} = 'H2: snapshotTime/time_delayed collision empties the child Tp slice.';
    end
    if isfield(R,'inport_mismatch') && R.inport_mismatch
        R.verdict{end+1} = 'Inport length mismatch at the suture instant.';
    end
    if R.sim_ok && isfield(R,'H2_confirmed') && ~R.H2_confirmed
        R.verdict{end+1} = 'H1 ALONE: the child branch simulates and reconstructs fine once sim(simIn2) is restored at DET_2.m:450.';
    end
    if isempty(R.verdict)
        R.verdict{end+1} = 'Inconclusive — see fields above.';
    end
    for k = 1:numel(R.verdict)
        fprintf('  [%d] %s\n', k, R.verdict{k});
    end

    outfile = sprintf('diagnose_child_1170_%s.mat', datestr(now,'yyyymmdd_HHMMSS'));
    save(outfile, 'R');
    fprintf('\nSaved diagnosis to %s\n', outfile);
    fprintf('======================================================================\n\n');
end

%% ----------------------------------------------------------------------
function t = safe_snaptime(op)
    t = NaN;
    try
        t = op.snapshotTime;
    catch
    end
end

function tf = isequal_tol(a, b)
    tf = false;
    if isnumeric(a) && isnumeric(b) && isscalar(a) && isscalar(b)
        tf = abs(a - b) < 1e-9;
    end
end

function names = dump_logged_states(op, label)
    names = {};
    try
        ls = op.loggedStates;
        n = numel(ls);
        fprintf('  %s: %d logged states\n', label, n);
        watch = [11 12 13 23];
        for i = 1:n
            nm = '<no name>';
            try
                if isfield(ls{i}, 'Name') || isprop(ls{i}, 'Name')
                    nm = ls{i}.Name;
                end
            catch
            end
            names{i} = nm;
            mark = '';
            if ismember(i, watch)
                mark = '   <== accessed by DET_2';
            end
            if ismember(i, watch) || i <= 3
                sz = '';
                try
                    sz = mat2str(size(ls{i}.Values.Data));
                catch
                end
                fprintf('    {%2d} %-38s %s%s\n', i, nm, sz, mark);
            end
        end
        if n < 23
            fprintf('    >> ONLY %d states: DET_2 accesses {23} -> INDEX OUT OF BOUNDS.\n', n);
        end
    catch ME
        fprintf('  %s: could not read loggedStates (%s)\n', label, ME.message);
    end
end
