%% instantiate_ga_params.m
% Instantiate parameters for GA_DET (to integrate into DET_scheduler_v1_perc_online)

% --- Luus–Jaakola reference parameters ---
max_passes      = 30;                      % number of passes
M               = 2;                       % iterations per pass
R               = 15;                      % samples per iteration
r_init          = 0.5 * (u_max - u_min);  % initial search radius
gamma           = 0.95;                   % radius decay

% --- Derived GA parameters ---
num_generations = 20;               % total GA generations
pop_size        = 50;               % GA population size (15)
crossover_frac  = 0.80;             % fraction of children by crossover
mutation_rate   = 0.30;             % initial mutation probability (0.05)

% --- Problem setup (assumes these exist in workspace) ---
%   u_min, u_max, N_pred, det_params,
%   target, d_mpc, time_snapshot_root,
%   delayed_time_rom, t_sam

% Initial guess and bounds for decision vector
u0     = (u_min + r_init) * ones(1, N_pred+1);
bounds = repmat([u_min, u_max], N_pred+1, 1);

