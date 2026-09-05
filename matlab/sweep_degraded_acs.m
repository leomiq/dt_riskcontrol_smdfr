function D = sweep_degraded_acs(matfile, vel_list, n_u)
% SWEEP_DEGRADED_ACS  Make tracking and safety COMPETE by degrading the ACS.
%
% Motivation. At the nominal operating points the DSM penalty either has little
% leverage (t=1180: tracking has a sharp interior optimum) or is redundant
% (t=1161: the PCT is 2 K above setpoint, so tracking already demands maximum
% coolant and saturates at u_max). To show what the deterministic safety margin
% adds, the PCT must be able to climb toward the limit DESPITE high coolant
% demand -- i.e. the auxiliary cooling system must be degraded.
%
% The archive has no degraded run in the current parameterisation
% (LJ_test_sampled_vel_ACS_09.mat is an older experiment: u_bounds [-2,2],
% N_pred 50, PCT BELOW setpoint, so DSM is negative there and unusable).
% This generates the condition instead, by lowering sampled_vel_ACS in the
% current-series workspace and re-running the DET.
%
%   matfile  : current-series workspace (default: the t=1161 acute-phase one)
%   vel_list : ACS velocities to test (default [2 1 0.5 0.25 0.1])
%   n_u      : points on the current grid (default 7)
%
% For each ACS velocity it reports the reachable peak-DSM band and whether any
% (dsm_thr, lambda) pair reorients the optimum -- and, decisively, whether the
% PCT still approaches the limit at u_max (safety and tracking in conflict).

    if nargin < 1 || isempty(matfile),  matfile = 'LJ_test_sampled_vel_ACS_2_loop_1161.mat'; end
    if nargin < 2 || isempty(vel_list), vel_list = [2 1 0.5 0.25 0.1]; end
    if nargin < 3 || isempty(n_u),      n_u = 7; end

    D = struct();
    D.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fprintf('\n===== Degraded-ACS DSM sweep =====\n');
    fprintf('Started %s\n\n', D.timestamp);

    S = load(matfile);
    fprintf('Regenerating operating points from %s...\n', matfile);
    dp0 = refresh_operating_points(S);
    if ~isfield(dp0, 'Y_denorm_50'), dp0.Y_denorm_50 = dp0.Y_denorm; end

    t_root  = double(dp0.time_snapshot_root);
    t_sam   = double(dp0.t_sam);
    delayed = double(dp0.delayed_time_rom);
    N_pred  = size(dp0.Y_denorm_50, 1);
    target  = CalculateTs(100, 543);
    U_Y     = 1244;
    span    = U_Y - target;

    u_min = 56.69; u_max = 56.75;
    u_grid = linspace(u_min, u_max, n_u);
    d_mpc  = ones(1, N_pred + 1);

    D.matfile = matfile; D.t_root = t_root; D.target = target;
    D.U_Y = U_Y; D.span = span; D.u_grid = u_grid; D.vel_list = vel_list;
    fprintf('query state    : t_root = %.1f s\n', t_root);
    fprintf('nominal vel_ACS: %.5f\n', double(dp0.sampled_vel_ACS));
    fprintf('ACS velocities : %s\n', mat2str(vel_list));
    fprintf('u grid         : %d points in [%.5f, %.5f] A\n', n_u, u_min, u_max);
    fprintf('target %.4f K, U_y %.1f K, DSM span %.4f K\n\n', target, U_Y, span);

    nV = numel(vel_list);
    D.peak_Tp   = nan(nV, n_u);
    D.peak_DSM  = nan(nV, n_u);
    D.J_error   = nan(nV, n_u);
    D.traj      = cell(nV, n_u);
    D.mdot      = nan(1, n_u);

    % IMPORTANT (verified 2026-09-03): DET_2 sets vel_ACS with
    %   setVariable(..., 'Workspace', name_Wks)
    % where name_Wks = 'simulation_core_dualfluid_FI3_sim_rom_extract_MCDET_v2'
    % is a .m INITIALISATION SCRIPT, not a Simulink model. That call therefore
    % has no effect on what the ROM simulates: the ROM keeps the vel_ACS = 0.1
    % its own model workspace was populated with at load time (line 132 of that
    % script). A sweep that only changes det_params.sampled_vel_ACS produces
    % identical results for every value -- which is exactly what a first run of
    % this script showed. The degradation must be written into the ROM's own
    % model workspace instead.
    rom_name = 'rom_MCDET_v2_2022b_control_CMEX';
    if ~bdIsLoaded(rom_name), load_system([rom_name '.slx']); end
    rom_wks = get_param(rom_name, 'ModelWorkspace');
    D.vel_baseline = getVariable(rom_wks, 'vel_ACS');
    fprintf('ROM model-workspace vel_ACS before sweep: %s\n\n', mat2str(D.vel_baseline, 8));

    for v = 1:nV
        vel = vel_list(v);
        fprintf('---- vel_ACS = %.4f ----\n', vel);
        dp = dp0;
        dp.sampled_vel_ACS = vel;
        % write it where the ROM actually reads it
        assignin(rom_wks, 'vel_ACS', vel);
        chk = getVariable(rom_wks, 'vel_ACS');
        fprintf('   ROM workspace vel_ACS now = %s\n', mat2str(chk, 8));
        for k = 1:n_u
            uv = u_grid(k);
            u_trial = uv * ones(1, N_pred + 1);
            try
                [fd, fr] = DET_2(dp, u_trial, d_mpc);
                opts = struct('dsm_thr', 1.0, 'lambda_pct', 0, 'U_y', U_Y);
                [Tp, ~, ~, tm] = fitness_cost_dsm(fd, fr, N_pred, target, ...
                    t_root, delayed, t_sam, u_trial, opts);
                D.traj{v,k}     = Tp;
                D.peak_Tp(v,k)  = tm.peak_Tp;
                D.peak_DSM(v,k) = tm.peak_DSM;
                D.J_error(v,k)  = tm.J_error;
                if v == 1
                    try, m = CurrentToMassFlow(u_trial); D.mdot(k) = m(1); catch, end
                end
            catch ME
                fprintf('   u=%.5f FAILED: %s\n', uv, ME.message);
            end
        end
        ok = ~isnan(D.peak_Tp(v,:));
        if any(ok)
            fprintf('   peak Tp: %.4f (u_min) -> %.4f (u_max) K\n', ...
                D.peak_Tp(v,1), D.peak_Tp(v,end));
            fprintf('   peak DSM band: [%.6f, %.6f]\n', ...
                min(D.peak_DSM(v,ok)), max(D.peak_DSM(v,ok)));
            [~, im] = min(D.J_error(v,ok));
            vi = find(ok);
            fprintf('   unpenalised optimum: u = %.5f A%s\n', u_grid(vi(im)), ...
                ternary(vi(im)==n_u, '  <<< SATURATED at u_max', ''));
            fprintf('   margin to limit at u_max: %.4f K\n', U_Y - D.peak_Tp(v,end));
        end
        fprintf('\n');
    end

    %% ---- Which velocity makes safety and tracking compete? ----
    fprintf('===== COMPETITION ANALYSIS =====\n');
    fprintf('%10s %12s %12s %12s %10s %12s\n', 'vel_ACS', 'DSM@u_min', 'DSM@u_max', ...
        'u*_unpen', 'saturated', 'margin@u_max');
    fprintf('%s\n', repmat('-', 1, 74));
    thr_grid = [0.20 0.10 0.05 0.02];
    lam_grid = [0 1e1 1e2 1e3 1e4 1e5];
    D.moves = zeros(nV, 1);
    D.best = cell(nV, 1);
    for v = 1:nV
        ok = ~isnan(D.peak_Tp(v,:));
        if ~any(ok), continue; end
        vi = find(ok);
        [~, im] = min(D.J_error(v,ok));
        sat = (vi(im) == n_u);
        fprintf('%10.4f %12.6f %12.6f %12.5f %10s %12.4f\n', ...
            vel_list(v), D.peak_DSM(v,1), D.peak_DSM(v,end), u_grid(vi(im)), ...
            ternary(sat,'YES','no'), U_Y - D.peak_Tp(v,end));

        % does any (thr, lambda) move the optimum for this velocity?
        base_u = u_grid(vi(im));
        best_shift = 0; best = [];
        for a = 1:numel(thr_grid)
            for b = 1:numel(lam_grid)
                thr = thr_grid(a); lam = lam_grid(b);
                Jt = nan(1, n_u);
                for k = 1:n_u
                    if isempty(D.traj{v,k}), continue; end
                    dsm = (D.traj{v,k} - target)/span;
                    Jt(k) = D.J_error(v,k) + lam*sum(max(0, dsm-thr).^2);
                end
                [~, iq] = min(Jt);
                sh = u_grid(iq) - base_u;
                if abs(sh) > abs(best_shift)
                    best_shift = sh;
                    best = struct('thr', thr, 'lam', lam, 'shift', sh, ...
                        'u_new', u_grid(iq), 'dsm_new', D.peak_DSM(v,iq));
                end
            end
        end
        D.moves(v) = best_shift;
        D.best{v} = best;
    end

    fprintf('\n===== REORIENTATION BY ACS VELOCITY =====\n');
    for v = 1:nV
        if isempty(D.best{v})
            fprintf('vel_ACS %.4f : no data\n', vel_list(v));
        elseif abs(D.best{v}.shift) < 1e-9
            fprintf('vel_ACS %.4f : optimum unchanged by any (thr, lambda)\n', vel_list(v));
        else
            b = D.best{v};
            fprintf('vel_ACS %.4f : thr %.2f, lambda %.0e -> u* %+.5f A (to %.5f), DSM -> %.6f\n', ...
                vel_list(v), b.thr, b.lam, b.shift, b.u_new, b.dsm_new);
        end
    end

    %% ---- Verdict ----
    fprintf('\n===== VERDICT =====\n');
    dsm_max_all = max(D.peak_DSM(:, end));   % DSM at u_max, worst velocity
    [worst_dsm, wi] = max(D.peak_DSM(:, end));
    fprintf('Highest DSM still present at FULL coolant (u_max): %.6f (vel_ACS %.4f)\n', ...
        worst_dsm, vel_list(wi));
    if worst_dsm > 0.05
        D.conclusion = sprintf(['At vel_ACS = %.4f the PCT still consumes %.1f%% of the ' ...
            'safety margin even at maximum coolant demand: tracking cannot resolve the ' ...
            'excursion on its own, so the DSM penalty is no longer redundant. This is the ' ...
            'operating condition that demonstrates the added value of the deterministic ' ...
            'safety margin.'], vel_list(wi), 100*worst_dsm);
    else
        D.conclusion = sprintf(['Even the most degraded ACS tested (vel_ACS = %.4f) keeps the ' ...
            'peak DSM at %.6f under full coolant, i.e. below 5%% of the margin. Either the ' ...
            'degradation must go further, or this query instant is past the thermal peak.'], ...
            vel_list(end), D.peak_DSM(end, end));
    end
    fprintf('\n%s\n', D.conclusion);

    % Restore the ROM workspace so the model is left as it was found.
    try
        assignin(rom_wks, 'vel_ACS', D.vel_baseline);
        fprintf('\nROM workspace vel_ACS restored to %s\n', mat2str(getVariable(rom_wks,'vel_ACS'), 8));
    catch ME
        fprintf('\n!! could not restore vel_ACS: %s\n', ME.message);
    end

    outfile = sprintf('sweep_degraded_acs_%s.mat', datestr(now,'yyyymmdd_HHMMSS'));
    save(outfile, 'D');
    fprintf('\nSaved %s\n==================================\n\n', outfile);
end

function out = ternary(c, a, b)
    if c, out = a; else, out = b; end
end
