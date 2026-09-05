function [Tp_reconstructed, J, feasible] = fitness_cost_penalised(forest_data, forest_roots, nu, target, ...
                                            time_snapshot_root, delayed_time_rom, t_sam, u_trial)
% FITNESS_COST_PENALISED  Drop-in replacement for fitness_cost with the DSM
% safety-margin penalty ACTIVE, for the penalised closed-loop realization.
%
% Signature is identical to fitness_cost so luus_jaakola_det / ga_det can call
% it unchanged; the penalty parameters are read from a global so the closed-loop
% driver can set them once per run.
%
%   global DSM_PENALTY_OPTS
%   DSM_PENALTY_OPTS = struct('dsm_thr', 0.00237737, 'lambda_pct', 1e4);
%
% Defaults (if the global is unset) reproduce the unpenalised behaviour, so an
% accidental call cannot silently change results.
%
% Penalty form, written on the DSM itself so lambda is dimensionless:
%       J_psm = lambda_pct * sum( max(0, DSM - dsm_thr)^2 )
%       DSM   = (Tp - target) / (U_y - target),  U_y = 1244 K
%
% The reconstruction is the fixed one: every node's segment is taken strictly
% after the warm-up window, so ONLY future (post-decision) instants enter the
% cost -- verified by audit_psm_future_only.m.

    global DSM_PENALTY_OPTS %#ok<GVMIS>
    if isempty(DSM_PENALTY_OPTS)
        dsm_thr = inf;          % never binds -> unpenalised
        lambda_pct = 0;
    else
        dsm_thr    = DSM_PENALTY_OPTS.dsm_thr;
        lambda_pct = DSM_PENALTY_OPTS.lambda_pct;
    end
    U_y = 1244;

    single_tree = forest_data{1};
    time_delayed = time_snapshot_root-delayed_time_rom:t_sam:time_snapshot_root-t_sam;

    Tp_reconstructed = [];
    t_reconstructed = [];
    feasible = true;

    %% ---- Reconstruction (future horizon only) ----
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

    %% ---- Tracking error ----
    nominal_value = target;
    err = Tp_reconstructed(:,1) - nominal_value;
    J_error = sum(err.^2);
    J = J_error;

    %% ---- Soft constraint 1: mass-flow rate ----
    m_seq = CurrentToMassFlow(u_trial);
    dm_seq = diff(m_seq);
    delta_limit = 5;
    lambda_mdot = 1E-4;
    penalty_mdot = sum(max(0, abs(dm_seq) - delta_limit).^2);
    J_mdot = lambda_mdot * penalty_mdot;
    J = J + J_mdot;

    %% ---- Soft constraint 2: DSM safety-margin penalty (ACTIVE) ----
    span = U_y - nominal_value;
    DSM = (Tp_reconstructed - nominal_value) ./ span;
    dsm_excess = max(0, DSM - dsm_thr);
    penalty_pct = sum(dsm_excess.^2);
    J_psm = lambda_pct * penalty_pct;
    J = J + J_psm;

    %% ---- Hard constraint ----
    if any(abs(DSM) > 1)
        feasible = false;
    end

    fprintf('J_error: %.6e, J_mdot: %.6e, J_psm: %.6e (peak DSM %.6f, %d viol)\n', ...
        J_error, J_mdot, J_psm, max(DSM), sum(DSM > dsm_thr));
end
