function E = e2_compare_optimizers(seeds, solvers, nlp_csv, query_file)
% E2_COMPARE_OPTIMIZERS  Corrected Luus-Jaakola and GA versus the monolithic NLP
% on the archived query t_q = 1180 s, all costed by the SAME evaluator
% (DET_2 + fitness_cost.m, i.e. the cost the closed loop actually optimizes).
%
% Settings are those of the production loop (run_phase2c_full_loop.m):
%   bounds [56.69, 56.75] A;  LJ: R = 24, gamma = 0.90, M = 5, 30 passes,
%   r_init = (u_max-u_min)/2, relative stop;  GA: pop 20 x 20 generations,
%   crossover 0.8, mutation 0.3;  start u0 = mid-range;  optimizer stream
%   rng(seed) before each call.
% The NLP action (julia/f2_monolithic_nlp.jl) is read from nlp_csv (column u_NLP)
% and evaluated once with DET_2.
%
%   seeds    default 1:3
%   solvers  default {'LJ','GA'}
%   nlp_csv  default results/f1_rom_replica/f2_nlp_vs_lj_1180.csv
%
% Writes results/e2_optimizers_1180/e2_summary.csv and e2_actions.csv, and
% e2_results.mat. Run from the repository root:
%   matlab -batch "addpath('matlab'); addpath('models'); e2_compare_optimizers"

    if nargin < 1 || isempty(seeds),   seeds = 1:3; end
    if nargin < 2 || isempty(solvers), solvers = {'LJ', 'GA'}; end
    if isequal(solvers, 'none'), solvers = {}; end   % only cost the NLP and archived actions
    repo = fileparts(fileparts(mfilename('fullpath')));
    try, suppress_simulink_prompts(); catch, end
    if nargin >= 4 && ~isempty(query_file)
        % query from build_query_from_archive.m (operating points already current)
        S = load(query_file);
        dp = S.det_params;
        S.time_snapshot_root = dp.time_snapshot_root; S.t_sam = dp.t_sam;
        S.delayed_time_rom = dp.delayed_time_rom;
        S.u_opt = [];                                  % no archived optimizer action
    else
        S = load(fullfile(repo, 'models', 'LJ_test_sampled_vel_ACS_2_loop_1180.mat'));
        dp = refresh_operating_points(S);
    end
    tq_tag = sprintf('%g', double(S.time_snapshot_root));
    if nargin < 3 || isempty(nlp_csv)
        nlp_csv = fullfile(repo, 'results', 'f1_rom_replica', ['f2_nlp_vs_lj_' tq_tag '.csv']);
    end
    outdir = fullfile(repo, 'results', ['e2_optimizers_' tq_tag]);
    if ~isfolder(outdir), mkdir(outdir); end
    if ~isfield(dp, 'Y_denorm_50'), dp.Y_denorm_50 = dp.Y_denorm; end
    for f = {'dsm_thr', 'lambda_dsm'}               % historical P0 cost only
        if isfield(dp, f{1}), dp = rmfield(dp, f{1}); end
    end

    tq = double(S.time_snapshot_root); t_sam = double(S.t_sam);
    delayed = double(S.delayed_time_rom);
    N_pred = size(dp.Y_denorm_50, 1); nu = N_pred + 1;
    target = CalculateTs(100, 543);
    u_min = 56.69; u_max = 56.75;
    u_bounds = repmat([u_min, u_max], nu, 1);
    d_mpc = ones(1, nu);
    max_passes = 30; M = 5; R = 24; gamma = 0.90; r_init = 0.5*(u_max - u_min);
    u0 = (u_min + r_init) * ones(1, nu);

    evaluate = @(u) cost_of(dp, u, d_mpc, nu, target, tq, delayed, t_sam);

    rows = {};  acts = struct();
    % Resume: keep finished runs from a previous (interrupted) call
    prev = fullfile(outdir, 'e2_results.mat');
    if isfile(prev)
        P = load(prev);
        rows = table2cell(P.E); acts = P.acts;
        fprintf('Resuming: %d runs already in %s\n', size(rows, 1), prev);
    end
    done = @(name) isfield(acts, name);
    % --- NLP action and archived (pre-fix) LJ action, one evaluation each
    T = readtable(nlp_csv);
    named = {'NLP', T.u_NLP(:)'};
    if ~isempty(S.u_opt), named(end+1, :) = {'LJ_archived', double(S.u_opt(:))'}; end
    for i = 1:size(named, 1)
        if done(named{i, 1}), continue; end
        u = min(max(named{i, 2}, u_min), u_max);
        tt = tic; [J, pk] = evaluate(u);
        rows(end+1, :) = {named{i, 1}, NaN, J, pk, u(1), 1, toc(tt), NaN}; %#ok<AGROW>
        acts.(named{i, 1}) = u;
        fprintf('[%s] J = %.6f  peak = %.4f K  u0 = %.5f A\n', named{i, 1}, J, pk, u(1));
    end

    % --- the derivative-free optimizers of the loop
    for s = 1:numel(solvers)
        for seed = seeds
            name = sprintf('%s_seed%d', solvers{s}, seed);
            if done(name), fprintf('[%s] already done, skipped\n', name); continue; end
            rng(seed, 'twister');
            tt = tic;
            if strcmpi(solvers{s}, 'LJ')
                [~, u, f, ~, tvec] = luus_jaakola_det(u0, u_bounds, max_passes, M, ...
                    r_init, gamma, R, dp, target, d_mpc, tq, delayed, t_sam);
                n_eval = 1 + nnz(tvec) * M * R;
            else
                [~, u, f, ~, tvec] = ga_det(u0, u_bounds, 20, 20, 0.80, 0.30, ...
                    dp, target, d_mpc, tq, delayed, t_sam);
                n_eval = 20 * 21;
            end
            wall = toc(tt);
            [J, pk] = evaluate(u);                  % re-evaluate (sanity: J == f)
            rows(end+1, :) = {solvers{s}, seed, J, pk, u(1), n_eval, wall, f}; %#ok<AGROW>
            acts.(name) = u;
            fprintf('[%s] J = %.6f (optimizer reported %.6f)  peak = %.4f K  u0 = %.5f A  %d evals  %.0f s\n', ...
                name, J, f, pk, u(1), n_eval, wall);
            save_all(outdir, rows, acts, tq, t_sam, nu);   % checkpoint after each run
        end
    end
    E = save_all(outdir, rows, acts, tq, t_sam, nu);
    disp(E);
end

function [J, peak] = cost_of(dp, u, d_mpc, nu, target, tq, delayed, t_sam)
    [fd, fr] = DET_2(dp, u, d_mpc);
    [Tp, J] = fitness_cost(fd, fr, nu, target, tq, delayed, t_sam, u);
    peak = max(Tp(:, 1));
end

function E = save_all(outdir, rows, acts, tq, t_sam, nu)
    E = cell2table(rows, 'VariableNames', ...
        {'method', 'seed', 'J', 'peak_PCT', 'u0', 'n_evals', 'wall_s', 'J_reported'});
    writetable(E, fullfile(outdir, 'e2_summary.csv'));
    t = (tq + (0:nu-1) * t_sam)';
    A = table(t);
    fn = fieldnames(acts);
    for i = 1:numel(fn), A.(fn{i}) = acts.(fn{i})(:); end
    writetable(A, fullfile(outdir, 'e2_actions.csv'));
    save(fullfile(outdir, 'e2_results.mat'), 'E', 'acts');
end
