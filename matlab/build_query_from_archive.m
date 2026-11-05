function qfile = build_query_from_archive(t0)
% BUILD_QUERY_FROM_ARCHIVE  Rebuild an MPC query at instant t0 (< 1180 s) from the
% archived closed-loop realization stored in LJ_test_sampled_vel_ACS_2_loop_1180.mat.
%
% The archive holds the accident scenario (sampled_A, B, drift, vel_ACS), the
% currents actually applied by the loop from 1160 s on (control_series) and the
% HRM output up to 1180 s. A query at t0 needs:
%   - the HRM operating points at t0 and t0 - 10 s, with the controls applied up to
%     t0 - 0.5 s (control_series truncated);
%   - the ROM operating point at t0 - 10 s (measured leak, no control);
%   - a leak forecast over the horizon t0 + 0.5 ... t0 + 12.5 s.
% The BiLSTM forecast cannot be regenerated here (no torch, and mc_inference.py
% loads a .pth that is not in the repository), so the forecast is the leak the
% HRM actually produced in the archived realization (perfect forecast). For a
% comparison of optimizers this is immaterial: every method receives the same
% series. Column 2 (ΔPCT forecast) is set to 0; it is not connected in the ROM.
%
% Output: results/queries/query_<t0>.mat with det_params (DET_2 format), target,
% u_bounds, d_mpc and provenance notes. Run from the repository root:
%   matlab -batch "addpath('matlab'); addpath('models'); build_query_from_archive(1161)"

    repo = fileparts(fileparts(mfilename('fullpath')));
    outdir = fullfile(repo, 'results', 'queries');
    if ~isfolder(outdir), mkdir(outdir); end
    try, suppress_simulink_prompts(); catch, end

    S = load(fullfile(repo, 'models', 'LJ_test_sampled_vel_ACS_2_loop_1180.mat'));
    ts = double(S.t_sam); dly = double(S.delayed_time_rom);
    N_pred = double(S.N_pred);
    assert(t0 + N_pred*ts <= double(S.time_snapshot_root), ...
        'horizon of t0 = %g s goes past the archived HRM run (1180 s)', t0);
    snapshot_id = find(abs(S.time_snapshots - t0) < 1e-9);
    assert(~isempty(snapshot_id), 't0 = %g s is not a loop query instant', t0);

    name_hrm = 'SMDFR_hrm_MCDET_2022b_control_CMEX';
    name_rom = 'rom_MCDET_v2_2022b_control_CMEX';
    name_wks = 'simulation_core_dualfluid_FI3_sim_rom_extract_MCDET_v2';

    % controls applied before t0 (row 1 of control_series is a dummy row)
    cs = S.control_series;
    cs = cs([true; cs(2:end, 1) <= t0 - ts + 1e-9], :);

    %% HRM at t0 and t0 - 10 s
    load_system(name_hrm);
    simRoot = Simulink.SimulationInput(name_hrm);
    simRoot = setVariable(simRoot, 'B', 10, 'Workspace', name_hrm);
    simRoot = setVariable(simRoot, 'A', double(S.sampled_A), 'Workspace', name_hrm);
    simRoot = setVariable(simRoot, 'mdot_c_drift', 10, 'Workspace', name_hrm);
    simRoot = setVariable(simRoot, 'vel_ACS', double(S.sampled_vel_ACS), 'Workspace', name_hrm);
    simRoot = simRoot.setExternalInput(cs(2:end, :));
    sim1 = setModelParameter(simRoot, 'StopTime', string(t0), 'SaveFinalState', 'on', ...
        'SaveOperatingPoint', 'on', 'FinalStateName', 'finalOPRoot');
    outRoot = sim(sim1);
    sim2 = setModelParameter(simRoot, 'StopTime', string(t0 - dly), 'SaveState', 'off', ...
        'SaveFinalState', 'on', 'SaveOperatingPoint', 'on', 'FinalStateName', 'finalOPRoot2');
    outRoot2 = sim(sim2);
    close_system(name_hrm, 0);

    % consistency: the re-run HRM must match the archived one up to t0
    la = S.outRoot.mdot_c_leak_signal; ln = outRoot.mdot_c_leak_signal;
    ka = la.Time <= t0 + 1e-9;
    dleak = max(abs(double(la.Data(ka)) - double(ln.Data(:))));
    fprintf('HRM re-run vs archive up to t0: max |Δ leak| = %.3e kg/s\n', dleak);

    %% ROM at t0 - 10 s (measured leak, no control), as the loop does
    load_system([name_rom '.slx']);
    simrom = Simulink.SimulationInput(name_rom);
    simrom = setModelParameter(simrom, 'FastRestart', 'off', 'SaveFinalState', 'on', ...
        'SaveOperatingPoint', 'on', 'FinalStateName', 'finalOP', 'StopTime', string(t0 - dly));
    tt = (0:ts:t0)'; n = numel(tt);
    simrom = simrom.setExternalInput([tt, double(ln.Data(:)), zeros(n, 3)]);
    simrom = setVariable(simrom, 'vel_ACS', double(S.sampled_vel_ACS), 'Workspace', name_wks);
    outrom = sim(simrom);
    close_system(name_rom, 0);

    %% Perfect leak forecast from the archived HRM run
    tf = t0 + (1:N_pred)' * ts;
    fut = interp1(double(la.Time), double(la.Data(:)), tf);
    Yf = [fut, zeros(N_pred, 1)];

    %% det_params in DET_2 format
    dp = S.det_params;
    dp.finalOP_Raiz = outRoot.finalOPRoot;
    dp.finalOP2_Raiz = outRoot2.finalOPRoot2;
    dp.romOP = outrom.finalOP;
    dp.outRoot = outRoot;
    dp.outrom = outrom;
    dp.time_snapshot_root = t0;
    dp.tsim = t0 + N_pred*ts;
    dp.snapshot_id = snapshot_id;
    dp.control_series = cs;
    dp.Y_denorm_50 = Yf;
    dp.Y_denorm_95 = Yf;
    dp.Y_denorm = reshape(Yf, [1 N_pred 2]);
    for f = {'dsm_thr', 'lambda_dsm'}
        if isfield(dp, f{1}), dp = rmfield(dp, f{1}); end
    end

    Q = struct();
    Q.det_params = dp;
    Q.target = CalculateTs(100, 543);
    Q.u_bounds = repmat([56.69, 56.75], N_pred + 1, 1);
    Q.d_mpc = ones(1, N_pred + 1);
    Q.forecast = 'perfect: archived HRM leak (BiLSTM not reproducible locally)';
    Q.hrm_leak_check = dleak;
    Q.applied_before_t0 = cs(end-min(5, size(cs,1)-2):end, :);
    qfile = fullfile(outdir, sprintf('query_%g.mat', t0));
    save(qfile, '-struct', 'Q');
    fprintf('Saved %s (snapshot %d, leak at t0 = %.2f kg/s, forecast %.2f .. %.2f kg/s)\n', ...
        qfile, snapshot_id, ln.Data(end), min(fut), max(fut));
end
