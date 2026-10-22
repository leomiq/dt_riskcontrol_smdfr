function result = m1_branch_and_bound(domain, bound_fn, cost_fn, opts)
% Conditional branch-and-bound. bound_fn(box) must cover every control and
% every reachable DET topology in the box. This requirement is external.
% cost_fn(u) must return [cost,feasible] for the SAME model and objective.
    if nargin < 4, opts = struct(); end
    if ~isfield(opts,'abs_gap'), opts.abs_gap = 1e-3; end
    if ~isfield(opts,'max_splits'), opts.max_splits = 2000; end
    if ~isfield(opts,'initial_u'), opts.initial_u = []; end
    if ~isnumeric(domain) || ~isreal(domain) || size(domain,2) ~= 2 || isempty(domain) || ...
       any(~isfinite(domain(:))) || any(domain(:,1)>=domain(:,2)) || ...
       ~isscalar(opts.abs_gap) || ~isfinite(opts.abs_gap) || opts.abs_gap < 0 || ...
       ~isscalar(opts.max_splits) || ~isfinite(opts.max_splits) || opts.max_splits < 0 || ...
       opts.max_splits ~= floor(opts.max_splits)
        error('Invalid domain or options.');
    end
    leaves = struct('box',domain,'lb',checked_bound(bound_fn,domain));
    U = inf; best_u = []; evaluations = 0;
    if ~isempty(opts.initial_u)
        u = opts.initial_u(:);
        if numel(u) ~= size(domain,1) || any(u<domain(:,1)) || any(u>domain(:,2))
            error('Initial control outside domain.');
        end
        [U,best_u] = check_candidate(cost_fn,u,U,best_u);
        evaluations = evaluations+1;
    end
    [U,best_u] = check_candidate(cost_fn,mean(domain,2),U,best_u);
    evaluations = evaluations+1;
    splits = 0;
    while splits < opts.max_splits
        [L,idx] = min([leaves.lb]);
        if isfinite(U) && U-L <= opts.abs_gap, break; end
        box = leaves(idx).box;
        [width,dim] = max(box(:,2)-box(:,1));
        mid = box(dim,1)+width/2;
        if mid <= box(dim,1) || mid >= box(dim,2), break; end
        left = box; right = box;
        left(dim,2) = mid; right(dim,1) = mid;
        a = struct('box',left,'lb',checked_bound(bound_fn,left));
        b = struct('box',right,'lb',checked_bound(bound_fn,right));
        leaves(idx) = a; leaves(end+1) = b; %#ok<AGROW>
        [U,best_u] = check_candidate(cost_fn,mean(left,2),U,best_u);
        [U,best_u] = check_candidate(cost_fn,mean(right,2),U,best_u);
        evaluations = evaluations+2; splits = splits+1;
    end
    L = min([leaves.lb]);
    if isfinite(U) && L > U + 1e-10*max(1,abs(U))
        error('Lower bound exceeds an evaluated feasible cost. Check the bound oracle.');
    end
    result = struct('lower_bound',L,'upper_bound',U,'gap',U-L, ...
        'best_u',best_u,'splits',splits,'leaves',numel(leaves), ...
        'evaluations',evaluations,'tolerance_met', ...
        isfinite(U) && U-L<=opts.abs_gap, ...
        'status','CONDITIONAL: verified enclosures and exact arithmetic required', ...
        'pending_boxes',leaves);
end

function lb = checked_bound(fn,box)
    lb = fn(box);
    if ~isscalar(lb) || ~isreal(lb) || ~isfinite(lb) || lb < 0
        error('Invalid lower bound.');
    end
end

function [U,best_u] = check_candidate(fn,u,U,best_u)
    [J,feasible] = fn(u);
    if ~isscalar(feasible) || ~ismember(feasible,[0 1])
        error('Invalid feasibility flag.');
    end
    if feasible
        if ~isscalar(J) || ~isreal(J) || ~isfinite(J) || J < 0
            error('Invalid candidate cost.');
        end
        if J < U, U = J; best_u = u; end
    end
end
