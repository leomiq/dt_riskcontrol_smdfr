function [lb, terms] = certified_cost_lower_bound_box(T_box, m_box, target, opts)
% CERTIFIED_COST_LOWER_BOUND_BOX Lower bound for the current MPC cost.
% T_box: N-by-2 certified enclosures [min,max] of EVERY future temperature
% sample reconstructed by fitness_cost for a fixed DET topology.
% m_box: K-by-2 certified enclosures of mass flow for each control step.
% The caller, not this function, must certify enclosure validity and topology.
% A singleton box recovers the exact cost, up to floating-point roundoff.
%
% Defaults match fitness_cost.m. Set opts.dsm_thr and opts.lambda_dsm
% together to select the active DSM penalty. This function does not certify
% the plant, scenario probabilities, feasibility, or global optimality alone.
    if nargin < 4 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'U_y'), opts.U_y = 1244; end
    if ~isfield(opts,'delta_limit'), opts.delta_limit = 5; end
    if ~isfield(opts,'lambda_mdot'), opts.lambda_mdot = 1e-4; end
    has_dsm = isfield(opts,'dsm_thr') && isfield(opts,'lambda_dsm');
    if ~has_dsm && (isfield(opts,'dsm_thr') || isfield(opts,'lambda_dsm'))
        error('Provide dsm_thr and lambda_dsm together.');
    end
    if has_dsm
        if ~isfield(opts,'lambda_pct'), opts.lambda_pct = 10; end
        safety_weight = opts.lambda_dsm;
    else
        if ~isfield(opts,'lambda_pct'), opts.lambda_pct = 10; end
        safety_weight = opts.lambda_pct;
    end
    if ~isnumeric(T_box) || size(T_box,2) ~= 2 || isempty(T_box) || ...
            any(~isfinite(T_box(:))) || any(T_box(:,1) > T_box(:,2))
        error('T_box must be a nonempty finite N-by-2 interval array.');
    end
    if ~isnumeric(m_box) || size(m_box,2) ~= 2 || isempty(m_box) || ...
            any(~isfinite(m_box(:))) || any(m_box(:,1) > m_box(:,2))
        error('m_box must be a nonempty finite K-by-2 interval array.');
    end
    if ~isscalar(target) || ~isfinite(target) || target >= opts.U_y || ...
            ~isscalar(opts.U_y) || ~isfinite(opts.U_y)
        error('Require finite target < U_y.');
    end
    weights = [opts.delta_limit opts.lambda_mdot safety_weight];
    if any(~isfinite(weights)) || any(weights < 0) || ...
            ~isscalar(opts.delta_limit) || ~isscalar(opts.lambda_mdot) || ...
            ~isscalar(safety_weight)
        error('Limits and penalty weights must be finite nonnegative scalars.');
    end
    if has_dsm && (~isscalar(opts.dsm_thr) || ~isfinite(opts.dsm_thr))
        error('dsm_thr must be a finite scalar.');
    end

    % Distance from the target to the temperature interval.
    tracking_distance = max([T_box(:,1)-target, target-T_box(:,2), ...
                             zeros(size(T_box,1),1)],[],2);
    terms.tracking = sum(tracking_distance.^2);

    % Any true consecutive flow difference lies in [d_lo,d_hi].
    if size(m_box,1) < 2
        terms.flow = 0;
    else
        d_lo = m_box(2:end,1) - m_box(1:end-1,2);
        d_hi = m_box(2:end,2) - m_box(1:end-1,1);
        excess = max([d_lo-opts.delta_limit, -opts.delta_limit-d_hi, ...
                      zeros(numel(d_lo),1)],[],2);
        terms.flow = opts.lambda_mdot * sum(excess.^2);
    end

    % Both safety hinges increase with temperature, so the lower endpoint
    % minimizes each one over the enclosing interval.
    if has_dsm
        dsm_lo = (T_box(:,1)-target)/(opts.U_y-target);
        excess_T = max(0,dsm_lo-opts.dsm_thr);
        terms.safety_mode = 'dsm';
    else
        excess_T = max(0,T_box(:,1)-opts.U_y);
        terms.safety_mode = 'absolute';
    end
    terms.safety = safety_weight * sum(excess_T.^2);
    lb = terms.tracking + terms.flow + terms.safety;
end
