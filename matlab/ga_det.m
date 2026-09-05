function [Y_best, u_best, f_best, history, t_vector] = ga_det(u0, bounds, num_generations, pop_size, crossover_frac, mutation_rate, det_params, target, d_mpc, time_snapshot_root, delayed_time_rom, t_sam)
% GA_DET: Genetic Algorithm optimization wrapper for DET-based fitness
%   Similar interface to luus_jaakola_det.m for use in DET_scheduler routines.
%
% Inputs:
%   u0              - Initial guess (1 x nu)
%   bounds          - [nu x 2] matrix of [lower_bound, upper_bound] for each variable
%   num_generations - Number of GA generations (analogous to max_passes)
%   pop_size        - Population size (analogous to R)
%   crossover_frac  - Fraction of population to create via crossover (0-1)
%   mutation_rate   - Probability of mutation per gene (0-1)
%   det_params      - Struct with parameters for DET model (as in DET_scheduler)
%   target          - Target for fitness cost
%   d_mpc           - Discrete MPC indicator vector
%   time_snapshot_root - Current time snapshot (for fitness_cost)
%   delayed_time_rom   - Delay for ROM (for fitness_cost)
%   t_sam           - Sampling time (for fitness_cost)
%
% Outputs:
%   Y_best    - Best output trajectory corresponding to u_best
%   u_best    - Best decision vector found
%   f_best    - Best fitness (cost) value
%   history   - Best cost per generation (num_generations x 1)
%   t_vector  - Elapsed time per generation (num_generations x 1)

%   Similar interface to luus_jaakola_det.m; uses parfor to speed up population evaluations.

    nu = length(u0);
    % Bounds
    lb = bounds(:,1)'; ub = bounds(:,2)';

    % Initialize parallel pool if needed
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('local', min(pop_size, feature('numcores')));
    end

    % Initialize population
    pop = repmat(u0, pop_size, 1) + (rand(pop_size,nu) - 0.5) .* repmat((ub-lb), pop_size, 1);
    pop = max(min(pop, ub), lb);

    %disp(pop');

    % Evaluate initial population in parallel
    fvals = inf(pop_size,1);
    Y_pool = cell(pop_size,1);
    parfor i = 1:pop_size
        [data_i, root_i] = DET_2(det_params, pop(i,:), d_mpc);
        [Y_i, f_i, feasible] = fitness_cost(data_i, root_i, nu, target, time_snapshot_root, delayed_time_rom, t_sam, pop(i,:));
        if feasible
            fvals(i) = f_i;
            Y_pool{i} = Y_i;
        end
    end

    % Set initial best
    [f_best, idx] = min(fvals);
    if isfinite(f_best)
        u_best = pop(idx,:);
        Y_best = Y_pool{idx};
    else
        u_best = u0;
        [data0, root0] = DET_2(det_params, u_best, d_mpc);
        [Y_best, f_best, ~] = fitness_cost(data0, root0, nu, target, time_snapshot_root, delayed_time_rom, t_sam, u_best);
    end

    history = zeros(num_generations,1);
    t_vector = zeros(num_generations,1);

    % Main GA loop
    for gen = 1:num_generations
        gen_tic = tic();

        % Selection (roulette wheel based on inverse cost)
        inv_f = 1 ./ (1 + fvals);
        probs = inv_f / sum(inv_f);
        cumprob = cumsum(probs);

        % Create next generation (sequential or parallel)
        new_pop = zeros(size(pop));
        for k = 1:pop_size
            % Parent selection
            r1 = rand(); parent1 = pop(find(cumprob>=r1,1),:);
            r2 = rand(); parent2 = pop(find(cumprob>=r2,1),:);
            % Crossover
            child = parent1;  % default: inherit parent1 (archived code left
                              % `child` undefined when crossover was skipped)
            if rand() < crossover_frac
                mask = rand(1, nu) < 0.5;  % Binary mask
                child = parent1.*mask + parent2.*(1-mask);
            end

            % Mutation
            mut_mask = (rand(1,nu) < mutation_rate);
            child(mut_mask) = lb(mut_mask) + rand(1,sum(mut_mask)) .* (ub(mut_mask)-lb(mut_mask));
            % Clamp
            new_pop(k,:) = max(min(child, ub), lb);
        end
        [~, best_idx] = min(fvals);
        new_pop(1, :) = pop(best_idx, :);  % Preserve best individual
        pop = new_pop;

        % Evaluate new population in parallel
        fvals = inf(pop_size,1);
        parfor i = 1:pop_size
            [data_i, root_i] = DET_2(det_params, pop(i,:), d_mpc);
            [Y_i, f_i, feasible] = fitness_cost(data_i, root_i, nu, target, time_snapshot_root, delayed_time_rom, t_sam, pop(i,:));
            if feasible
                fvals(i) = f_i;
                Y_pool{i} = Y_i;
            end
        end

        % Update best
        [gen_best, idx] = min(fvals);
        if gen_best < f_best
            f_best = gen_best;
            u_best = pop(idx,:);
            Y_best = Y_pool{idx};
        end
        history(gen) = f_best;
        t_vector(gen) = toc(gen_tic);
        disp(sprintf('Generation %d complete. Best fitness: %g', gen, f_best));
        %disp(u_best');
    end
end
