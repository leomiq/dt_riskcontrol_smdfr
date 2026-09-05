% fitness_cost.m
function [Tp_reconstructed, J, feasible] = fitness_cost(forest_data, forest_roots, nu, target,...
                                            time_snapshot_root, delayed_time_rom, t_sam, u_trial)
% FITNESS_CODE Reconstructs 'Tp_signal' from forest_data{1}'s sequential nodes 
%              and computes a cost measure J as the sum of squared errors from 1144.
%              Also checks feasibility if at each t, the deviation is within ±0.3.
%
% Inputs:
%   forest_data  - Cell array of forest structures; forest_data{1} holds a single tree
%   forest_roots - Typically the initial states for each forest (unused here)
%
% Outputs:
%   Tp_reconstructed - Two-column matrix: [time, Tp_signal_values]
%   J                - Scalar cost (sum of squared deviations from 1144)
%   feasible         - Boolean indicating if margin constraints are satisfied

    % We assume forest_data{1} is an array, each element having .tree_data
    single_tree = forest_data{1};

    %adjusting root node
    time_delayed = time_snapshot_root-delayed_time_rom:t_sam:time_snapshot_root-t_sam;

    % Initialize outputs
    Tp_reconstructed = [];
    t_reconstructed = [];
    feasible = true;

    % 1) Reconstruct the 'Tp_signal'

    % FIX 2026-09-03: reconstruction of the child-node segments.
    %
    % DET branching is driven by the ACS decision: whether the auxiliary cooling
    % system must be activated, and once activated whether it is available or
    % failed -- i.e. the [demand_ps, ava_ps] pair tracked by
    % msfcn_schedule_rom_v2. A branch cut short by such an event carries a SHORT
    % segment that need not reach time_delayed(end) (measured: a complementary
    % node halting one t_sam past the suture at t=1170.5). The previous code lost
    % those segments three ways:
    %   (a) it searched find(t_signal == time_delayed(end)) on nodes that never
    %       reach that instant -> empty index -> Tp_signal(idx+1:end) empty;
    %   (b) child nodes concatenated their FULL signal, re-adding the
    %       time_delayed warm-up window the root had already contributed
    %       -> duplicated samples between 1170 and 1179.5;
    %   (c) with init_state_root == [1 1] only even node_ids were accepted, so
    %       odd-id children (the transition branch, id 1) were dropped silently.
    %
    % Every node's segment is now taken strictly after the branch instant, which
    % is the node's own start time for children and time_delayed(end) for the
    % root warm-up. Nodes are no longer filtered by id parity.
    %
    % NOTE: when no transition time falls inside the prediction window, the tree
    % legitimately stays at the root alone and this loop reduces to the original
    % single-node behaviour. That is the expected outcome, not a failure.
    root_init_state = single_tree(1).tree_data.init_state;
    for i = 1:length(single_tree)
        % 'tree_data' is a cell array with node-data structures
        node_data_list = single_tree(i).tree_data;
        node_id = single_tree(i).tree.id;

        Tp_signal = node_data_list.data.pred_tw_max_signal.Data(:);
        t_signal  = node_data_list.data.pred_tw_max_signal.Time(:);
        if isempty(t_signal)
            continue;
        end

        if node_id == 0
            % Root: drop the ROM warm-up window [time_delayed(1) .. time_delayed(end)].
            % Kept conditional on the root state, as in the original code.
            if isequal(root_init_state, [1,1])
                keep = t_signal > time_delayed(end);
            else
                keep = true(size(t_signal));
            end
        else
            % Child: keep everything strictly after the warm-up window. A branch
            % that stopped early simply contributes fewer samples instead of none.
            keep = t_signal > time_delayed(end);
            if ~any(keep)
                % Short branch that halted inside/at the warm-up window: keep the
                % samples past its own branching instant so the node still counts.
                keep = t_signal > t_signal(1);
            end
        end

        Tp_reconstructed = cat(1, Tp_reconstructed, Tp_signal(keep));
        t_reconstructed  = cat(1, t_reconstructed,  t_signal(keep));
    end

    % NOTE: LENGTH CAN BE SHORTER THAN PREDICTION HORIZON!!!
    %Tp_reconstructed = Tp_reconstructed(1:nu);

    %% === Objective cost (J) (deviation from nominal) ===
    nominal_value = target;
    if ~isempty(Tp_reconstructed)
        err = Tp_reconstructed(:,1) - nominal_value;
        % No weighting
        %J   = sum(err.^2);
        % Early weighting
        time_horizon = length(Tp_reconstructed);
        %weights = exp(-0.2*(0:time_horizon-1))'; 
        weights = 1;
        J_error = sum(weights.*(err).^2);
    else
        % If no data was found, we can set cost to 0 or large
        J_error = 0;
        feasible = false;  % Or handle differently
        return;
    end
    J = J_error;
    %% === SOFT CONSTRAINT 1: ∆mass flow rate (J_mdot) ===
    m_seq = CurrentToMassFlow(u_trial);       % <-- Analytical model
    dm_seq = diff(m_seq);
    delta_limit =  5;  % kg/s per step, customary 5
    lambda_mdot = 1E-4; % penalty weight, customary 0
    penalty_mdot = sum(max(0, abs(dm_seq) - delta_limit).^2);
    J_mdot = lambda_mdot * penalty_mdot;
    J = J + J_mdot;

    %% === SOFT CONSTRAINT 2: PCT margin soft penalty (J_psm) ===
    % Default (historical) form: penalise only an ABSOLUTE limit violation,
    % Tp > U_y. With the setpoint 64.7 K below U_y this term is identically zero
    % in every scenario tested, so lambda_pct has no effect (lambda * 0 = 0).
    %
    % OPTIONAL DSM form, enabled by setting the global DSM_PENALTY_OPTS:
    %     global DSM_PENALTY_OPTS
    %     DSM_PENALTY_OPTS = struct('dsm_thr', 0.00237737, 'lambda_pct', 1e4);
    % which instead penalises the deterministic safety margin exceeding a
    % threshold,  J_psm = lambda_pct * sum(max(0, DSM - dsm_thr)^2),  with
    % DSM = (Tp - target)/(U_y - target). Writing it on the DSM keeps lambda
    % dimensionless so the same value means the same thing across operating
    % points. Leaving the global unset reproduces the historical behaviour
    % exactly, so existing results are unaffected.
    % The settings travel inside forest_data{1}(1).tree_data (attached by DET_2
    % from det_params) rather than through a global, because luus_jaakola_det
    % evaluates the cost inside a `parfor` (line 71) and MATLAB does NOT
    % propagate client globals to parallel workers -- a global here silently
    % yields J_psm = 0 on every worker.
    U_y = 1244;  % Absolute safety limit
    dsm_thr = [];
    lambda_dsm = [];
    try
        rootnode = single_tree(1).tree_data;
        if isfield(rootnode, 'dsm_thr') && isfield(rootnode, 'lambda_dsm')
            dsm_thr = rootnode.dsm_thr;
            lambda_dsm = rootnode.lambda_dsm;
        end
    catch
    end

    if isempty(dsm_thr) || isempty(lambda_dsm)
        pct_violation = max(0, Tp_reconstructed - U_y);
        lambda_pct = 10; % customary value, 100
        penalty_pct = sum(pct_violation.^2);
    else
        span_dsm = U_y - nominal_value;
        DSM_pen = (Tp_reconstructed - nominal_value) ./ span_dsm;
        lambda_pct = lambda_dsm;
        penalty_pct = sum(max(0, DSM_pen - dsm_thr).^2);
    end
    J_psm = lambda_pct * penalty_pct;
    J = J + J_psm;

    %% HARD CONSTRAINT. Check feasibility: if absolute deviation from nominal exceeds safety margin, it's infeasible
    U_y = 1244*ones(size(Tp_reconstructed)); 
    nominal_vector = nominal_value*ones(size(Tp_reconstructed));
    DSM = (Tp_reconstructed - nominal_vector) ./ (U_y - nominal_vector);
    abs_margin = 1;
    if any(abs(DSM) > abs_margin)
        feasible = false;
    end
    J_disp = sprintf('J_error: %d, J_mdot: %d, J_psm: %d.', J_error, J_mdot, J_psm);
    disp(J_disp);
end

