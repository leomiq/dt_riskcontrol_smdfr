function R = run_penalised_trajectory(dsm_thr, lambda_pct, seed)
% RUN_PENALISED_TRAJECTORY  One full 61-query closed-loop realization with the
% DSM safety-margin penalty ACTIVE, to be compared against the existing
% unpenalised pilot (phase2c_fullloop_pilot_LJ_seed1.mat).
%
% Experimental design: identical solver (LJ), identical seed, identical accident
% stream and identical bounds. The ONLY difference is the cost function -- the
% DSM penalty term. That isolates the penalty's effect on the closed-loop PCT
% trajectory and on the coolant the controller injects.
%
%   dsm_thr    DSM threshold above which the penalty applies (default 0.00237737,
%              the value that reorients the decision at the t=1180 fixed state)
%   lambda_pct penalty weight (default 1e4 -- the setting where J_psm reaches
%              ~26% of the cost, i.e. comparable to J_error rather than
%              dominating it as lambda=1e6 does)
%   seed       optimizer seed (default 1, matching the baseline pilot)
%
% Cost: about 6 h (the baseline pilot recorded 22242.6 s = 6.18 h of optimizer
% time over its 61 queries). Run it in the background.
%
% Baseline for comparison (already on disk):
%   PCT first 1181.2132, last 1179.1560, min 1146.5177, max 1181.2132
%   max DSM 0.030002 at t=1160.0, last DSM -0.001780
%   u range [56.69000, 56.75000]

    if nargin < 1 || isempty(dsm_thr),    dsm_thr = 0.00237737; end
    if nargin < 2 || isempty(lambda_pct), lambda_pct = 1e4;     end
    if nargin < 3 || isempty(seed),       seed = 1;             end

    fprintf('\n===== Penalised closed-loop realization =====\n');
    fprintf('started %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf('dsm_thr = %.8f, lambda_pct = %.0e, seed = %d\n', dsm_thr, lambda_pct, seed);
    fprintf('baseline (lambda=0) already on disk: phase2c_fullloop_pilot_LJ_seed1.mat\n\n');

    try, suppress_simulink_prompts(); catch, end

    % Activate the DSM penalty. It travels through det_params -> DET_2 ->
    % tree_data root node -> fitness_cost, NOT through a global: the cost is
    % evaluated inside a parfor in luus_jaakola_det (line 71) and MATLAB does
    % not propagate client globals to parallel workers. A first attempt using a
    % global ran 50 min and produced J_psm = 0 on every evaluation.
    setenv('DSM_THR', num2str(dsm_thr, 12));
    setenv('DSM_LAMBDA', num2str(lambda_pct, 12));
    fprintf('DSM penalty settings staged for run_phase2c_full_loop\n');
    fprintf('   (read from the environment and written into det_params)\n\n');

    tag = sprintf('penalised_thr%.5f_lam%.0e', dsm_thr, lambda_pct);
    tag = strrep(tag, '.', 'p');
    tag = strrep(tag, '+', '');

    ttic = tic;
    KPI = run_phase2c_full_loop('LJ', seed, tag);
    elapsed = toc(ttic);

    R = struct();
    R.dsm_thr = dsm_thr;
    R.lambda_pct = lambda_pct;
    R.seed = seed;
    R.tag = tag;
    R.KPI = KPI;
    R.elapsed_hours = elapsed/3600;
    fprintf('\nrealization finished in %.2f h\n', R.elapsed_hours);

    outfile = sprintf('run_penalised_%s_%s.mat', tag, datestr(now,'yyyymmdd_HHMMSS'));
    save(outfile, 'R');
    fprintf('saved %s\n', outfile);

    % clear the staged settings so a later unpenalised run is unaffected
    setenv('DSM_THR', '');
    setenv('DSM_LAMBDA', '');
    fprintf('DSM penalty settings cleared\n');
    fprintf('compare with the baseline using compare_penalised_vs_baseline.py\n');
    fprintf('=============================================\n\n');
end
