function run_retuned_lam1e4()
% RUN_RETUNED_LAM1E4  Closed-loop realization with the retuned LJ and a
% moderated penalty weight.
%
% WHY THESE SETTINGS. Two independent problems were found in the widened-bound
% run (u_max = 56.90, lambda = 1e6):
%
%   1. LJ was returning sub-optimal actions. Its pass reset was
%      r = |u_best - old_u_best|, which discarded the within-pass contraction
%      and let the radius grow again, so it never resolved the narrow interior
%      minimum (~0.05 A wide, cost x179 at its edges). The early-stop test was
%      absolute (improvement < 1) while late-horizon costs reach ~3000.
%      Measured gaps versus a grid reference: +161% to +18933%.
%      Fixed in luus_jaakola_det.m (geometric ceiling on the radius, relative
%      early stop) plus gamma 0.95 -> 0.90 and R 15 -> 24.
%
%   2. lambda = 1e6 was too large. With the retuned LJ the true optimum at the
%      late queries is u_max, i.e. maximum coolant even when the PCT is already
%      below setpoint -- so the over-cooling was NOT only an optimizer artefact.
%      Fixed-state check with the retuned LJ:
%          t=1180  lambda 1e4 -> u* 56.76173 (J_psm 0.3%)   INTERIOR
%                  lambda 1e5 -> u* 56.75133 (J_psm 7.7%)
%          t=1188  lambda 1e4 -> u* 56.69000 (J_psm 62.5%)  u_min, correct when
%                                                            already below setpoint
%                  lambda 1e5 -> u* 56.86340 (J_error 266!) pays a large tracking
%                                                            error for a negligible
%                                                            penalty
%      lambda = 1e4 keeps the penalty significant without forcing saturation.
%
% Expected outcome: coherent control (no 0.58 <-> 43.55 kg/s thrashing) and a PCT
% that converges toward the setpoint instead of drifting below it.

    fprintf('\n########## Retuned LJ, lambda = 1e4 ##########\n');
    fprintf('started %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf('u_max = 56.90, dsm_thr = 0.00237737, lambda = 1e4\n');
    fprintf('LJ: gamma = 0.90, R = 24, relative early stop (0.1%% of f_best)\n');
    fprintf('##############################################\n\n');

    try, suppress_simulink_prompts(); catch, end

    setenv('U_MAX_OVERRIDE', '56.90');
    setenv('DSM_THR', '0.00237737');
    setenv('DSM_LAMBDA', '1e4');

    t0 = tic;
    try
        K = run_phase2c_full_loop('LJ', 1, 'retuned_lam1e4_umax56p90');
        fprintf('\nfinished in %.2f h\n', toc(t0)/3600);
        save('KPI_retuned_lam1e4.mat', 'K');
    catch ME
        fprintf('\n*** FAILURE ***\n%s\n%s\n', ME.identifier, ME.message);
        for k = 1:numel(ME.stack)
            fprintf('  at %s line %d\n', ME.stack(k).name, ME.stack(k).line);
        end
        fprintf('(checkpoint kept -- rerun to resume)\n');
    end

    setenv('U_MAX_OVERRIDE', '');
    setenv('DSM_THR', '');
    setenv('DSM_LAMBDA', '');
    fprintf('\n########## done %s ##########\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
end
