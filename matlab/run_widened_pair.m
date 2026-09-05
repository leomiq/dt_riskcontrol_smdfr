function run_widened_pair(u_max_new, dsm_thr, lambda_pct, seed)
% RUN_WIDENED_PAIR  The two closed-loop realizations that isolate the DSM
% penalty's effect once the MHD pump is allowed to deliver more coolant.
%
%   (1) baseline  : lambda = 0,        u_max = u_max_new
%   (2) penalised : lambda = lambda_pct, u_max = u_max_new
%
% Both share solver (LJ), seed, accident stream and bounds, so the cost function
% is the only difference. The widened baseline is essential: without it a change
% in the PCT trajectory could not be attributed to the penalty rather than to
% simply having more range available.
%
% Why widen at all: with the archived bounds [56.69, 56.75] the control was
% SATURATED at u_max through most of the window where the penalty dominated the
% cost (J_psm/J_error up to 197 over 36 of 61 queries), so the safety term had
% already obtained everything it could ask for -- saturation, not the penalty,
% was the binding constraint. Fixed-state verification at u_max = 56.90 shows
% the penalised optimum becomes INTERIOR (56.7425 A, +3.58 kg/s, J_psm 70% of
% the cost), i.e. the term finally has room to act.
%
% Each realization takes ~4 h and checkpoints after every query, so a power loss
% costs one query rather than the run. Rerunning this function resumes.

    if nargin < 1 || isempty(u_max_new),  u_max_new = 56.90;      end
    if nargin < 2 || isempty(dsm_thr),    dsm_thr = 0.00237737;   end
    if nargin < 3 || isempty(lambda_pct), lambda_pct = 1e6;       end
    if nargin < 4 || isempty(seed),       seed = 1;               end

    fprintf('\n########## Widened-bound realization pair ##########\n');
    fprintf('started %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf('u_max = %.5f A (was 56.75), seed = %d\n', u_max_new, seed);
    fprintf('penalty: dsm_thr = %.8f, lambda = %.0e\n', dsm_thr, lambda_pct);
    fprintf('####################################################\n\n');

    try, suppress_simulink_prompts(); catch, end
    setenv('U_MAX_OVERRIDE', num2str(u_max_new, 12));

    tagw = strrep(sprintf('umax%.2f', u_max_new), '.', 'p');

    %% ---- (1) widened baseline, no penalty ----
    fprintf('\n===== (1/2) BASELINE  lambda = 0,  u_max = %.5f =====\n', u_max_new);
    setenv('DSM_THR', '');
    setenv('DSM_LAMBDA', '');
    tag1 = sprintf('widebase_%s', tagw);
    t1 = tic;
    try
        K1 = run_phase2c_full_loop('LJ', seed, tag1);
        fprintf('\nbaseline finished in %.2f h\n', toc(t1)/3600);
        save(sprintf('KPI_%s_LJ_seed%d.mat', tag1, seed), 'K1', 'u_max_new', 'seed');
    catch ME
        fprintf('\n*** BASELINE FAILED: %s\n%s\n', ME.identifier, ME.message);
        for k = 1:numel(ME.stack)
            fprintf('   at %s line %d\n', ME.stack(k).name, ME.stack(k).line);
        end
        fprintf('(its checkpoint is kept -- rerun to resume)\n');
    end

    %% ---- (2) widened + penalty ----
    fprintf('\n===== (2/2) PENALISED  lambda = %.0e,  u_max = %.5f =====\n', ...
        lambda_pct, u_max_new);
    setenv('DSM_THR', num2str(dsm_thr, 12));
    setenv('DSM_LAMBDA', num2str(lambda_pct, 12));
    tag2 = sprintf('widepen_%s', tagw);
    t2 = tic;
    try
        K2 = run_phase2c_full_loop('LJ', seed, tag2);
        fprintf('\npenalised finished in %.2f h\n', toc(t2)/3600);
        save(sprintf('KPI_%s_LJ_seed%d.mat', tag2, seed), 'K2', ...
            'u_max_new', 'dsm_thr', 'lambda_pct', 'seed');
    catch ME
        fprintf('\n*** PENALISED FAILED: %s\n%s\n', ME.identifier, ME.message);
        for k = 1:numel(ME.stack)
            fprintf('   at %s line %d\n', ME.stack(k).name, ME.stack(k).line);
        end
        fprintf('(its checkpoint is kept -- rerun to resume)\n');
    end

    %% ---- tidy up ----
    setenv('DSM_THR', '');
    setenv('DSM_LAMBDA', '');
    setenv('U_MAX_OVERRIDE', '');
    fprintf('\n########## pair complete %s ##########\n', ...
        datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf('outputs: %s_LJ_seed%d.mat  and  %s_LJ_seed%d.mat\n', tag1, seed, tag2, seed);
    fprintf('then run plot_penalty_comparison.py with these two files\n');
end
