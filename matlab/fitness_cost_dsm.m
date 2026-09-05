function [Tp_reconstructed, J, feasible, terms] = fitness_cost_dsm(forest_data, forest_roots, nu, target, ...
                                            time_snapshot_root, delayed_time_rom, t_sam, u_trial, opts)
% FITNESS_COST_DSM  Parameterised variant of fitness_cost with the safety-margin
% penalty written on the DSM itself instead of on (Tp - U_y) in Kelvin.
%
% Rationale: the deterministic safety margin
%       DSM = (Tp - target) / (U_y - target)
% is the dimensionless quantity the constraint is meant to act on. Writing the
% penalty as
%       J_psm = lambda_pct * sum(max(0, DSM - dsm_thr).^2)
% makes lambda dimensionless too, so the same lambda means the same thing across
% operating points and thresholds. The production fitness_cost.m keeps its
% Kelvin form (equivalent to dsm_thr = 1 up to a scale factor); this variant
% exists for the (dsm_thr, lambda) sweep and does not replace it.
%
% opts fields (all optional):
%   dsm_thr      DSM threshold above which the penalty applies   (default 1.0)
%   lambda_pct   weight on the DSM excess                        (default 10)
%   lambda_mdot  weight on the mass-flow-rate soft constraint    (default 1e-4)
%   U_y          absolute safety limit [K]                       (default 1244)
%   verbose      print the J breakdown                           (default false)
%
% terms returns the decomposition so a sweep can report each contribution.

    if nargin < 9 || isempty(opts), opts = struct(); end
    if ~isfield(opts, 'dsm_thr'),     opts.dsm_thr = 1.0;     end
    if ~isfield(opts, 'lambda_pct'),  opts.lambda_pct = 10;   end
    if ~isfield(opts, 'lambda_mdot'), opts.lambda_mdot = 1e-4; end
    if ~isfield(opts, 'U_y'),         opts.U_y = 1244;        end
    if ~isfield(opts, 'verbose'),     opts.verbose = false;   end

    single_tree = forest_data{1};
    time_delayed = time_snapshot_root-delayed_time_rom:t_sam:time_snapshot_root-t_sam;

    Tp_reconstructed = [];
    t_reconstructed = [];
    feasible = true;
    terms = struct('J_error', NaN, 'J_mdot', NaN, 'J_psm', NaN, ...
                   'peak_Tp', NaN, 'peak_DSM', NaN, 'n_violating', 0, ...
                   'n_samples', 0, 'penalty_raw', 0);

    %% ---- Reconstruction (identical to the fixed fitness_cost.m) ----
    % Future horizon only: every node's segment is taken strictly after the
    % warm-up window, so past instants never enter the cost.
    root_init_state = single_tree(1).tree_data.init_state;
    for i = 1:length(single_tree)
        node_data_list = single_tree(i).tree_data;
        node_id = single_tree(i).tree.id;
        Tp_signal = node_data_list.data.pred_tw_max_signal.Data(:);
        t_signal  = node_data_list.data.pred_tw_max_signal.Time(:);
        if isempty(t_signal), continue; end

        if node_id == 0
            if isequal(root_init_state, [1,1])
                keep = t_signal > time_delayed(end);
            else
                keep = true(size(t_signal));
            end
        else
            keep = t_signal > time_delayed(end);
            if ~any(keep), keep = t_signal > t_signal(1); end
        end
        Tp_reconstructed = cat(1, Tp_reconstructed, Tp_signal(keep));
        t_reconstructed  = cat(1, t_reconstructed,  t_signal(keep));
    end

    if isempty(Tp_reconstructed)
        J = 0; feasible = false; return;
    end
    terms.n_samples = numel(Tp_reconstructed);

    %% ---- Tracking error ----
    nominal_value = target;
    err = Tp_reconstructed(:,1) - nominal_value;
    J_error = sum(err.^2);
    terms.J_error = J_error;
    J = J_error;

    %% ---- Soft constraint 1: mass-flow rate ----
    m_seq = CurrentToMassFlow(u_trial);
    dm_seq = diff(m_seq);
    delta_limit = 5;
    penalty_mdot = sum(max(0, abs(dm_seq) - delta_limit).^2);
    J_mdot = opts.lambda_mdot * penalty_mdot;
    terms.J_mdot = J_mdot;
    J = J + J_mdot;

    %% ---- Soft constraint 2: DSM safety-margin penalty ----
    span = opts.U_y - nominal_value;
    DSM = (Tp_reconstructed - nominal_value) ./ span;
    dsm_excess = max(0, DSM - opts.dsm_thr);
    penalty_raw = sum(dsm_excess.^2);
    J_psm = opts.lambda_pct * penalty_raw;

    terms.peak_Tp = max(Tp_reconstructed);
    terms.peak_DSM = max(DSM);
    terms.n_violating = sum(DSM > opts.dsm_thr);
    terms.penalty_raw = penalty_raw;
    terms.J_psm = J_psm;
    J = J + J_psm;

    %% ---- Hard constraint: DSM must not exceed 1 (Tp at the limit) ----
    if any(abs(DSM) > 1)
        feasible = false;
    end

    if opts.verbose
        fprintf('J_error: %.6e, J_mdot: %.6e, J_psm: %.6e  (peak DSM %.6f, %d viol)\n', ...
            J_error, J_mdot, J_psm, terms.peak_DSM, terms.n_violating);
    end
end
