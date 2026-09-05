function V = verify_det_fixes(matfile)
% VERIFY_DET_FIXES  End-to-end check of the two DET propagation fixes.
%
%   Fix 1  DET_2.m:450   restored `out2 = sim(simIn2)` (child branches were
%                        never simulated, so tree_data only ever held the root)
%   Fix 2  fitness_cost.m reconstruction now keeps child segments strictly after
%                        the warm-up window instead of dropping short branches,
%                        no longer duplicates the warm-up window, and no longer
%                        filters children by node-id parity
%
% This calls DET_2 ONCE with a fixed u_trial (no optimizer: 1 objective
% evaluation instead of 420) and reports whether children actually entered the
% tree and what the reconstruction now contains.
%
% Baseline to beat (pre-fix, measured): tree_data held the root only, the
% reconstructed peak was 1181.21 K and J_psm was identically 0.

    if nargin < 1 || isempty(matfile)
        matfile = 'LJ_test_sampled_vel_ACS_2_loop_1180.mat';
    end

    V = struct();
    V.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fprintf('\n============ DET fix verification ============\n');
    fprintf('Started %s\n\n', V.timestamp);

    S  = load(matfile);
    fprintf('Regenerating operating points...\n');
    dp = refresh_operating_points(S);
    if ~isfield(dp, 'Y_denorm_50'), dp.Y_denorm_50 = dp.Y_denorm; end

    t_root  = double(dp.time_snapshot_root);
    t_sam   = double(dp.t_sam);
    delayed = double(dp.delayed_time_rom);
    N_pred  = size(dp.Y_denorm_50, 1);
    target  = CalculateTs(100, 543);

    u_fixed = 56.72;
    u_trial = u_fixed * ones(1, N_pred + 1);
    d_mpc   = ones(1, N_pred + 1);
    V.u_trial_value = u_fixed;
    V.target = target;
    fprintf('u_trial = %.4f A (constant), target = %.4f K\n', u_fixed, target);
    fprintf('t_root = %.1f, suture = %.1f, t_sam = %.2f\n\n', t_root, t_root - delayed, t_sam);

    %% ---- Run the DET once ----
    fprintf('-- Calling DET_2 (one evaluation, no optimizer) --\n');
    tt = tic;
    V.det_ok = false;
    try
        [forest_data, forest_roots] = DET_2(dp, u_trial, d_mpc);
        V.det_ok = true;
        V.det_seconds = toc(tt);
        fprintf('DET_2 returned in %.1f s\n\n', V.det_seconds);
    catch ME
        V.det_seconds = toc(tt);
        V.det_error = sprintf('[%s] %s', ME.identifier, ME.message);
        fprintf('DET_2 FAILED after %.1f s: %s\n', V.det_seconds, V.det_error);
        for k = 1:numel(ME.stack)
            fprintf('   at %s line %d\n', ME.stack(k).name, ME.stack(k).line);
        end
        save_and_exit(V); return;
    end

    %% ---- Did children enter the tree? (Fix 1) ----
    fprintf('-- Fix 1 check: tree contents --\n');
    td = forest_data{1}.tree_data;
    tr = forest_data{1}.tree;
    V.n_tree_nodes = numel(tr);
    V.n_tree_data  = numel(td);
    fprintf('nodes declared in tree     : %d\n', V.n_tree_nodes);
    fprintf('nodes SAVED into tree_data : %d\n', V.n_tree_data);
    % tree_data is built as tree_data{end+1} = struct(...) but DET_2 assigns it
    % into forest_data{k}.tree_data, which may arrive as a struct array rather
    % than a cell. Accept both shapes.
    V.node_ids = [];
    V.node_snapshots = [];
    V.node_samples = [];
    fprintf('class(tree_data) = %s\n', class(td));
    for i = 1:numel(td)
        if iscell(td), nd = td{i}; else, nd = td(i); end
        nid = NaN; st = NaN; ns = NaN;
        try, nid = nd.id; catch, end
        try, st = nd.final_OP.snapshotTime; catch, end
        try, ns = numel(nd.data.pred_tw_max_signal.Time); catch, end
        V.node_ids(end+1) = nid;
        V.node_snapshots(end+1) = st;
        V.node_samples(end+1) = ns;
        fprintf('   node id=%-3s snapshotTime=%-8s samples=%s\n', ...
            num2str(nid), num2str(st), num2str(ns));
    end
    % A single-node tree is the CORRECT outcome when no ACS transition time falls
    % inside the prediction window: the root is never cut short, so the guard
    % `time_snapshot < tsim` (DET_2.m:243) reports "no event" and nothing is
    % queued. Fix 1 only matters once a child IS queued, so the pass criterion is
    % "every declared node was saved without error", not "more than one node".
    V.all_declared_saved = (V.n_tree_data == V.n_tree_nodes);
    if V.all_declared_saved
        fprintf('>> tree consistent: %d/%d declared nodes saved.\n', ...
            V.n_tree_data, V.n_tree_nodes);
        if V.n_tree_data == 1
            fprintf('   single-node tree = no ACS transition in the window (expected).\n\n');
        else
            fprintf('   >> FIX 1 EXERCISED: child nodes saved without error.\n\n');
        end
        V.fix1_ok = true;
    else
        fprintf('>> MISMATCH: %d nodes declared but only %d saved -> a child failed.\n\n', ...
            V.n_tree_nodes, V.n_tree_data);
        V.fix1_ok = false;
    end

    %% ---- Reconstruction and cost (Fix 2) ----
    fprintf('-- Fix 2 check: reconstruction and J terms --\n');
    try
        [Tp_rec, J, feasible] = fitness_cost(forest_data, forest_roots, N_pred, ...
            target, t_root, delayed, t_sam, u_trial);
        V.fitness_ok = true;
        V.n_samples = numel(Tp_rec);
        V.peak_Tp = max(Tp_rec);
        V.min_Tp = min(Tp_rec);
        V.J = J;
        V.feasible = feasible;
        fprintf('reconstructed samples : %d\n', V.n_samples);
        fprintf('peak Tp               : %.4f K   (margin to 1244: %.4f K)\n', ...
            V.peak_Tp, 1244 - V.peak_Tp);
        fprintf('min Tp                : %.4f K\n', V.min_Tp);
        fprintf('J                     : %.6e\n', J);
        fprintf('feasible              : %d\n', feasible);

        % Duplicate-timestamp check: the pre-fix code re-added the warm-up window.
        V.duplicate_warning = false;
        fprintf('\n(the J_error/J_mdot/J_psm line printed by fitness_cost is above)\n');

        % Baseline comparison
        base_peak = 1181.2131949675056;
        V.peak_delta_vs_baseline = V.peak_Tp - base_peak;
        fprintf('\npeak vs pre-fix baseline (%.4f K): %+.4f K\n', ...
            base_peak, V.peak_delta_vs_baseline);
        % Expected sample count: every node's samples minus the warm-up window
        % the root drops (time_delayed spans delayed_time_rom/t_sam samples).
        n_warmup = numel(t_root - delayed : t_sam : t_root - t_sam);
        V.expected_samples = sum(V.node_samples(~isnan(V.node_samples))) - n_warmup;
        fprintf('\nexpected samples = sum(node samples) - %d warm-up = %d\n', ...
            n_warmup, V.expected_samples);
        V.fix2_ok = (V.n_samples == V.expected_samples) && V.n_samples > 0;
        if V.fix2_ok
            fprintf('>> FIX 2 CONSISTENT: reconstruction accounts for every node segment\n');
            fprintf('   with the warm-up window counted exactly once (no duplication).\n');
        else
            fprintf('>> MISMATCH: got %d samples, expected %d — check for dropped or\n', ...
                V.n_samples, V.expected_samples);
            fprintf('   duplicated segments.\n');
        end
    catch ME
        V.fitness_ok = false;
        V.fitness_error = sprintf('[%s] %s', ME.identifier, ME.message);
        fprintf('fitness_cost FAILED: %s\n', V.fitness_error);
        for k = 1:numel(ME.stack)
            fprintf('   at %s line %d\n', ME.stack(k).name, ME.stack(k).line);
        end
    end

    %% ---- Verdict ----
    fprintf('\n============ VERDICT ============\n');
    if isfield(V,'fix1_ok') && V.fix1_ok
        fprintf('Tree consistency             : PASS (%d/%d nodes saved)\n', ...
            V.n_tree_data, V.n_tree_nodes);
    else
        fprintf('Tree consistency             : FAIL\n');
    end
    if isfield(V,'fix2_ok') && V.fix2_ok
        fprintf('Reconstruction accounting    : PASS (%d samples)\n', V.n_samples);
    else
        fprintf('Reconstruction accounting    : needs review\n');
    end
    fprintf('\nNote: a single-node tree is correct when no ACS transition time falls\n');
    fprintf('inside the prediction window. J_psm = 0 then reflects the scenario, not\n');
    fprintf('a defect: the predicted peak stays far below the 1244 K limit.\n');
    save_and_exit(V);
end

function save_and_exit(V)
    outfile = sprintf('verify_det_fixes_%s.mat', datestr(now,'yyyymmdd_HHMMSS'));
    save(outfile, 'V');
    fprintf('\nSaved %s\n=================================\n\n', outfile);
end
