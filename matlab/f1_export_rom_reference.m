function f1_export_rom_reference(u_csv, tag, query_file)
% F1 (monolithic formulation plan): reference data for the algebraic ROM replica.
%
% Rebuilds the ROOT-node ROM simulation of the archived query t_q = 1180 s exactly
% as DET_2 does it (initial operating point = ROM OP with HRM states overwritten,
% 10 s warm-up with measured leak and applied current, then the BiLSTM leak
% forecast and the archived optimal current u_opt), re-runs it in normal mode with
% full state logging, and writes everything the replica needs to
%   results/f1_rom_replica/rom_reference_1180.mat
%
% Also reports how far this re-run is from the archived prediction
% (forest_data_opt), which isolates the R2024a/R2024b gap from replica error.
%
% Optional: f1_export_rom_reference(u_csv, tag) replaces u_opt by the column
% named u_<tag> of the CSV file u_csv (26 values) and saves
% rom_reference_1180_<tag>.mat, printing the P0 cost of fitness_cost.m. Used to
% re-simulate the NLP action of julia/f2_monolithic_nlp.jl in Simulink (V3).
%
% Usage (from the repository root):
%   matlab -batch "addpath('matlab'); addpath('models'); f1_export_rom_reference"
%   matlab -batch "addpath('matlab'); addpath('models'); f1_export_rom_reference('results/f1_rom_replica/f2_nlp_vs_lj_1180.csv','NLP')"
%
% A third argument takes a query built by build_query_from_archive.m instead of
% the archived 1180 s workspace; the output is then rom_reference_<t_q>.mat:
%   matlab -batch "addpath('matlab'); addpath('models'); f1_export_rom_reference('','','results/queries/query_1161.mat')"

repo = fileparts(fileparts(mfilename('fullpath')));
outdir = fullfile(repo, 'results', 'f1_rom_replica');
if ~isfolder(outdir), mkdir(outdir); end

if nargin >= 3 && ~isempty(query_file)
    % query built by build_query_from_archive.m: no archived optimizer action,
    % the reference run uses the mid-range current
    A = load(query_file);
    nu = size(A.u_bounds, 1);
    A.u_opt = mean(A.u_bounds(1, :)) * ones(1, nu);
    fprintf('Query %s (%s)\n', query_file, A.forecast);
else
    A = load(fullfile(repo, 'models', 'LJ_test_sampled_vel_ACS_2_loop_1180.mat'), ...
        'det_params', 'u_opt', 'd_mpc', 'forest_data_opt', 'target', 'u_bounds');
end
P = A.det_params;
u_trial = A.u_opt;
d_mpc = A.d_mpc;
suffix = '';
if nargin >= 2 && ~isempty(tag)
    T = readtable(u_csv);
    u_trial = T.(['u_' tag])';
    assert(numel(u_trial) == numel(A.u_opt), 'expected %d actions', numel(A.u_opt));
    suffix = ['_' tag];
    fprintf('Using u_%s from %s\n', tag, u_csv);
end

name_rom = 'rom_MCDET_v2_2022b_control_CMEX';
load_system(name_rom);
set_param(name_rom, 'FastRestart', 'off');

t_q   = P.time_snapshot_root;
dly   = P.delayed_time_rom;
t_sam = P.t_sam;
tsim  = P.tsim;

outRoot = P.outRoot; outrom = P.outrom;

%% Initial operating point, as DET_2 (root node)
% DET_2 overwrites ALL 13 ROM states with HRM values at t_q - 10 s, so the ROM
% operating point only serves as a container. The archived one cannot be reused:
% it was saved in R2024a from a model then named SMDFR_rom_..., so its block
% paths do not resolve in models/rom_MCDET_v2_2022b_control_CMEX.slx. A fresh
% container is built by running the ROM 0 -> t_q - 10 s with the measured leak and
% no control (as run_phase2c_full_loop does), and the archived HRM values are
% copied into it, matching states by block path and occurrence order.
simOP = Simulink.SimulationInput(name_rom);
simOP = setModelParameter(simOP, 'StopTime', string(t_q - dly), ...
    'SaveFinalState', 'on', 'SaveOperatingPoint', 'on', 'FinalStateName', 'finalOP');
tt = (0:t_sam:t_q)';
simOP = simOP.setExternalInput([tt, outRoot.mdot_c_leak_signal.Data(1:numel(tt)), ...
    zeros(numel(tt), 3)]);
simOP = setVariable(simOP, 'vel_ACS', P.sampled_vel_ACS, 'Workspace', P.name_Wks);
outOP = sim(simOP);
romOP_temp = outOP.finalOP;

Raiz_to_Rom = [1,2,4,5,6,7,3,8,9,10,22,23,24];
tail = @(p) regexprep(p, '^[^/]*/', '');           % drop the model name
arch_keys = cell(1, P.romOP.loggedStates.numElements);
for i = 1:numel(arch_keys)
    arch_keys{i} = tail(char(P.romOP.loggedStates{i}.BlockPath.getBlock(1)));
end
nL = romOP_temp.loggedStates.numElements;
assert(nL == numel(arch_keys), 'ROM operating point has %d states, archive %d', nL, numel(arch_keys));
seen = containers.Map();
for j = 1:nL
    k = tail(char(romOP_temp.loggedStates{j}.BlockPath.getBlock(1)));
    if isKey(seen, k), seen(k) = seen(k) + 1; else, seen(k) = 1; end
    hits = find(strcmp(arch_keys, k));
    i = hits(seen(k));                              % same occurrence order
    v = P.finalOP2_Raiz.loggedStates{Raiz_to_Rom(i)}.Values.Data;
    if any(i == [3, 4, 5]), v = mean(v(1:12)); end  % 12 HRM nodes -> 0-D ROM
    old = romOP_temp.loggedStates{j}.Values.Data;
    assert(numel(old) == numel(v), 'size mismatch at %s', k);
    romOP_temp.loggedStates{j}.Values.Data = reshape(cast(v, 'like', old), size(old));
end

%% External inputs, exactly as DET_2 (root node)
% The archived run predates the switch to the median forecast and stores only
% the 95th percentile; any forecast serves the replica check as long as both
% models receive the same series.
if isfield(P, 'Y_denorm_50')
    Yf = P.Y_denorm_50; forecast_used = 'Y_denorm_50';
else
    Yf = P.Y_denorm_95; forecast_used = 'Y_denorm_95';
end
fprintf('Leak forecast used: %s\n', forecast_used);
pred_leak = cat(2, outRoot.mdot_c_leak_signal.Data(end,1), Yf(:,1)');
pred_dtw  = Yf(:,2)';
time_delayed = t_q-dly:t_sam:t_q-t_sam;
time_ext = cat(2, time_delayed, t_q:t_sam:tsim);
leak_ext = cat(2, outRoot.mdot_c_leak_signal.Data(end-length(time_delayed):end-1,1)', pred_leak);
dtw_delayed = outRoot.tw_max_signal.Data(end-length(time_delayed):end,1) - ...
              outrom.tw_max_signal.Data(end-length(time_delayed):end,1);
dtw_ext = cat(2, dtw_delayed', pred_dtw);
cs = P.control_series;
idx = arrayfun(@(t) find(cs(:,1) == t, 1), time_delayed);
u_ext = cat(2, cs(idx,2)', u_trial);
d_ext = cat(2, cs(idx,3)', d_mpc);
inport_data = [time_ext', leak_ext', dtw_ext', u_ext', d_ext'];

%% Simulate with full state logging
simIn = Simulink.SimulationInput(name_rom);
simIn = setInitialState(simIn, romOP_temp);
% The archived operating points were saved in R2024a; allow restoring them in a
% newer release (a warning, not an error). See docs/pilot_reproducibility.md.
simIn = setModelParameter(simIn, 'NonCurrentReleaseOperatingPointMsg', 'warning', ...
    'OperatingPointInterfaceChecksumMismatchMsg', 'warning');
simIn = setModelParameter(simIn, 'StopTime', string(tsim), ...
    'SaveState', 'on', 'StateSaveName', 'xout', 'SaveFormat', 'Dataset', ...
    'SaveFinalState', 'on', 'SaveOperatingPoint', 'on', 'FinalStateName', 'finalOP');
simIn = simIn.setExternalInput(inport_data);
simIn = setVariable(simIn, 'vel_ACS', P.sampled_vel_ACS, 'Workspace', P.name_Wks);
cell_trans = arrayfun(@(i) [0 P.Lambda_t(i); P.Mu_t(i) 0], 1:P.Ncomps, 'UniformOutput', false);
mat_trans = reshape(reshape(cell2mat(cell_trans), 2, 2, 5), 1, 20);
simIn = setVariable(simIn, 'mat_trans', mat_trans, 'Workspace', P.name_root);
simIn = setVariable(simIn, 'parent_last_time', P.time_snapshots(P.snapshot_id), 'Workspace', P.name_root);
simIn = setVariable(simIn, 'initial_branch_state', P.finalOP2_Raiz.loggedStates{23}.Values.Data, 'Workspace', P.name_root);
out = sim(simIn);

%% Collect states (initial values and trajectories), by block name
x0 = struct('block', {}, 'value', {});
for i = 1:romOP_temp.loggedStates.numElements
    s = romOP_temp.loggedStates{i};
    x0(i).block = char(s.BlockPath.getBlock(1));
    x0(i).value = double(s.Values.Data(:));
end
xs = struct('block', {}, 'time', {}, 'data', {});
for i = 1:out.xout.numElements
    e = out.xout.getElement(i);
    xs(i).block = char(e.BlockPath.getBlock(1));
    xs(i).time  = e.Values.Time(:);
    xs(i).data  = double(squeeze(e.Values.Data));
end

sig = @(n) [out.(n).Time(:), double(squeeze(out.(n).Data))];
ref = struct();
ref.t_q = t_q; ref.delayed_time_rom = dly; ref.t_sam = t_sam; ref.tsim = tsim;
ref.inport_data = inport_data;
ref.forecast_used = forecast_used;
ref.target = A.target;        % PCT setpoint y_ref of the MPC cost
ref.u_bounds = A.u_bounds;    % current bounds used by the archived LJ run
ref.u_opt = A.u_opt;          % archived Luus-Jaakola action (26 values, t_q .. t_q+12.5 s)
ref.x0 = x0; ref.xs = xs;
names = {'tw_signal','tf_signal','tc_signal','M_c_signal','mdot_c_signal', ...
         'totalpower_signal','rho_tot_signal','pred_tw_max_signal','tw_max_signal', ...
         'mdot_c_net_signal','mdot_c_leak_signal','Tin_c_signal','Te_c_signal','Te_f_signal'};
for k = 1:numel(names)
    try
        ref.(names{k}) = sig(names{k});
    catch
        warning('signal %s not logged', names{k});
    end
end

%% Parameters, as seen by the model (single source of truth for the replica)
W = get_param(name_rom, 'ModelWorkspace');
pn = {'Lambda','beta_i','beta','lambda_i','N_0','Pi_0','P_0','lambda_c','lambda_e', ...
      'alfa_f','alfa_c','rho_0','rhoex_value','C_i0','Ce_i0','A_fw','A_wc','mdot_f','mdot_c', ...
      'M_f','M_w','M_c','c_pf','c_pw','c_pc','h_fw','h_wc','Tc_0','Te_c0','Tf_0','Te_f0', ...
      'Tw_0','Tin_f0','Tin_c0','Tbar_f0','Tbar_c0','poly_vector'};
par = struct();
for k = 1:numel(pn)
    par.(pn{k}) = double(getVariable(W, pn{k}));
end
ref.par = par;

%% Re-run vs archived prediction (R2024a run stored in the .mat)
mine = out.pred_tw_max_signal;
if isfield(A, 'forest_data_opt')
    tr = A.forest_data_opt{1}(1).tree_data.data.pred_tw_max_signal;
    ta = tr.Time(:); ya = double(tr.Data(:));
    [tc, ia, ib] = intersect(round(ta*2)/2, round(mine.Time(:)*2)/2);
    d = double(mine.Data(ib)) - ya(ia);
    keep = tc > t_q - t_sam;
    ref.archived_pred_tw_max = [ta, ya];
    fprintf(['Re-run vs archived pred_tw_max on the horizon (%d samples): ' ...
        'diff in [%.4f, %.4f] K, std %.2e K\n'], nnz(keep), min(d(keep)), max(d(keep)), std(d(keep)));
end
fprintf('Initial states (%d logged):\n', numel(x0));
for i = 1:numel(x0)
    fprintf('  %2d  %-60s  n=%d  first=%.6g\n', i, x0(i).block, numel(x0(i).value), x0(i).value(1));
end

%% P0 cost of this action on the Simulink output (fitness_cost.m, default branch)
y = double(mine.Data(mine.Time > t_q - t_sam));
m = CurrentToMassFlow(u_trial);
J_track = sum((y - A.target).^2);
J_mdot  = 1e-4 * sum(max(0, abs(diff(m)) - 5).^2);
J_pct   = 10 * sum(max(0, y - 1244).^2);
fprintf('P0 cost on Simulink: J = %.6f (track %.6f, mdot %.2e, pct %.2e), max PCT %.4f K\n', ...
    J_track + J_mdot + J_pct, J_track, J_mdot, J_pct, max(y));

fn = fullfile(outdir, sprintf('rom_reference_%g%s.mat', t_q, suffix));
save(fn, '-struct', 'ref', '-v7');
fprintf('Saved %s\n', fn);
end
