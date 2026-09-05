function [Y_best, u_best, f_best, history, t_vector] = ...
         luus_jaakola_det(u0, bounds, max_passes, M, r_init, gamma, R, det_params, target, d_mpc,...
         time_snapshot_root, delayed_time_rom, t_sam)

    % LUUS_JAAKOLA_DET_MULTIPASS_PAR:
    %   Illustrative multi-pass LJ optimization with parallelization in each iteration,
    %   and a display of best cost at the end of each pass.
    %
    % Inputs:
    %   u0         - Initial guess (1 x nu)
    %   bounds     - [nu x 2], [lower_bound upper_bound] for each variable
    %   max_passes - Maximum number of "passes"
    %   M          - Number of iterations per pass
    %   r_init     - Initial search radius (scalar or vector)
    %   gamma      - Radius decay factor (e.g. 0.95)
    %   R          - Number of random samples tested each iteration
    %   det_params - Struct with parameters for DET model
    %
    % Outputs:
    %   Y_best, u_best, f_best  - best solution found
    %   history                 - array of f_best after each pass
    %   t_vector                - array of elapsed time (in seconds) for each pass

    %% Basic setup
    nu = length(u0);   % dimension of decision variables
    u_best = u0;       % best solution so far

    % Evaluate initial guess
    [data0, root0] = DET_2(det_params, u_best, d_mpc);
    [Y_best, f_best, feasible0] = fitness_cost(data0, root0, nu, target,...
                                               time_snapshot_root, delayed_time_rom, t_sam, u_best);

    % ------ Display best cost found so far ------
    disp(['Pass ', num2str(0), ' complete. f_best0 = ', num2str(f_best)]);

    % Bookkeeping arrays
    history  = zeros(max_passes,1);
    t_vector = zeros(max_passes,1);

    % Tolerance logic
    %
    % RETUNED 2026-09-04. With the widened pump bound (u_max 56.75 -> 56.90) the
    % search range grew 3.5x (0.06 -> 0.21 A) and the ABSOLUTE threshold mu = 1
    % became meaningless: late-horizon costs reach J ~ 3000, so "no improvement
    % >= 1" fired while the solution was still hundreds of units from the
    % optimum. Measured consequence (verify_lj_suboptimality.m): at t_q = 1188
    % LJ returned u = 56.69 costing J = 3354 when u = 56.74 cost J = 18.7 --
    % 179x worse, with the optimum only 0.05 A away. Relative gaps of 161% to
    % 18932% were found across the drift region, which is what produced the
    % downward PCT drift.
    %
    % The threshold is now RELATIVE to the current best, so it scales with the
    % problem regardless of the bound width or the penalty weight.
    mu_rel           = 1e-3;   % stop when a pass improves f_best by < 0.1%
    mu_abs_floor     = 1e-6;   % guard for f_best ~ 0
    pass_no_improve  = 0;
    max_no_improve   = 5;      % stop if no improvement for 5 consecutive passes

    % Region size & pass-based re-initialization
    r = r_init;       % can be scalar or vector
    old_u_best = u_best;
    old_f_best = f_best;

    for pass = 1 : max_passes
        pass_tic = tic();  % measure time spent in this pass

        %% (A) Re-initialize region at start of pass
        %
        % RETUNED 2026-09-04. The original rule was r = |u_best - old_u_best|,
        % i.e. the radius was RESET every pass to how far the optimum had moved.
        % That makes gamma nearly irrelevant (the within-pass contraction at
        % line ~112 is discarded at every pass boundary) and, worse, it lets the
        % radius GROW again whenever a pass happens to take a large step. On the
        % widened range this kept the search coarse: the true optimum sits in a
        % narrow interior minimum (~0.05 A wide, cost rising 179x at its edges),
        % so a radius that never contracts below ~0.02 A cannot resolve it.
        %
        % Now the pass-to-pass reset is combined with a monotone geometric
        % ceiling: the radius may follow the optimum's motion, but never exceeds
        % r_init * gamma^(pass-1). With gamma = 0.90 over 30 passes the ceiling
        % reaches ~4e-3 A, an order of magnitude finer than the minimum's width.
        if pass > 1
            delta_move = abs(u_best - old_u_best);
            delta_move(delta_move < 1e-6) = 1e-6;     % avoid total collapse
            r_ceiling = r_init * gamma^(pass-1);
            r = min(delta_move, r_ceiling);
            % keep a floor so the search never becomes degenerate
            r = max(r, 1e-6);
        end

        %% (B) Do M iterations in this pass
        for iter = 1 : M
            improvedThisIter = false;

            % Allocate arrays to store trial results in parallel
            best_fvals = inf(R,1);
            best_uvals = zeros(R, nu);
            best_Yvals = cell(R,1);

            % ---------------- PARALLEL loop over R samples ----------------
            parfor trial = 1 : R
            %for trial = 1 : R
                % local copy so each worker has the same "u_best"
                u_trial = u_best;

                % random perturbation
                delta = (2*rand(1,nu) - 1) .* r;
                u_trial = u_trial + delta;

                % clamp to bounds
                for i = 1:nu
                    u_trial(i) = max(min(u_trial(i), bounds(i,2)), ...
                                     bounds(i,1));
                end

                % evaluate cost
                [forest_data, forest_roots] = DET_2(det_params, u_trial, d_mpc);
                [Y_trial, f_trial, feasible] = fitness_cost(forest_data, forest_roots, nu, target,...
                                              time_snapshot_root, delayed_time_rom, t_sam, u_trial);

                % store if feasible
                %if feasible
                    best_fvals(trial)   = f_trial;
                    best_uvals(trial,:) = u_trial;
                    best_Yvals{trial}   = Y_trial;
                %end
                % else remain inf
            end  % parfor

            % After parfor, find best of these R trials
            [f_trial_best, idx_best] = min(best_fvals);
            if f_trial_best < f_best
                % update global best
                f_best = f_trial_best;
                u_best = best_uvals(idx_best, :);
                Y_best = best_Yvals{idx_best};
                improvedThisIter = true;
            end

            % if no improvement, shrink region
            if ~improvedThisIter
                r = r * gamma;
            end
        end  % for iter = 1:M

        %% (C) End of pass: record best and check improvement
        history(pass)  = f_best;
        t_vector(pass) = toc(pass_tic);

        % ------ Display best cost found so far ------
        disp(['Pass ', num2str(pass), ' complete. f_best = ', num2str(f_best)]);

        % Relative improvement test (see the note at mu_rel above): an absolute
        % threshold cannot serve costs spanning 1e-2 to 1e4 across the horizon.
        mu_pass = max(mu_rel * abs(f_best), mu_abs_floor);
        if abs(old_f_best - f_best) < mu_pass
            pass_no_improve = pass_no_improve + 1;
        else
            pass_no_improve = 0;
        end

        if pass_no_improve >= max_no_improve
            disp(['Stopping early: improvement < ', num2str(mu_pass), ...
                  ' (', num2str(100*mu_rel), '% of f_best) for ', ...
                  num2str(max_no_improve),' passes.']);
            break;
        end

        % Save old best for next pass
        old_u_best = u_best;
        old_f_best = f_best;
    end % for pass

end
