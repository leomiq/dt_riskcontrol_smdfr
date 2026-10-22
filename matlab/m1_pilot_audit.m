function A = m1_pilot_audit()
% Audits the saved R2024b candidate without rerunning DET.
% Flow-only global bound assumes exact monotonic CurrentToMassFlow on bounds.
% Thermal bound remains zero: sampled trajectories cannot enclose all controls.
    S=load(fullfile('models','LJ_test_sampled_vel_ACS_2_loop_1180.mat'), ...
        'u_bounds','target');
    P=load(fullfile('results','pilot_r2024b_20261022.mat'));
    u=P.u_opt(:); bounds=S.u_bounds; target=double(S.target);
    T=P.Tp(:);
    if size(bounds,1)~=numel(u) || any(u<bounds(:,1)-1e-8) || ...
            any(u>bounds(:,2)+1e-8) || isempty(T)
        error('Candidate outside domain or missing temperatures.');
    end
    m=arrayfun(@CurrentToMassFlow,u);
    flow_lo=arrayfun(@CurrentToMassFlow,bounds(:,1));
    flow_hi=arrayfun(@CurrentToMassFlow,bounds(:,2));
    if any(flow_lo>flow_hi), error('Pump enclosure is not monotone.'); end
    A.J_tracking=sum((T-target).^2);
    A.J_flow=1e-4*sum(max(0,abs(diff(m))-5).^2);
    A.J_safety=10*sum(max(0,T-1244).^2);
    A.recomputed_cost=A.J_tracking+A.J_flow+A.J_safety;
    A.saved_cost=P.f_opt;
    A.consistent=abs(A.recomputed_cost-A.saved_cost) <= ...
        1e-8*max(1,abs(A.saved_cost));
    A.feasible=all(abs((T-target)/(1244-target))<=1);
    if ~A.consistent || ~A.feasible
        error('Saved candidate is not a valid upper bound for this audit.');
    end
    A.upper_bound=A.recomputed_cost;
    [A.lower_bound_flow,term]=certified_cost_lower_bound_box( ...
        [target target],[flow_lo flow_hi],target);
    A.lower_bound_flow=term.flow;
    A.lower_bound_thermal=0;
    A.lower_bound=A.lower_bound_flow;
    A.gap=A.upper_bound-A.lower_bound;
    A.status='CONDITIONAL flow bound; no validated ROM/DET thermal enclosure';
    fprintf('M1 pilot: Jtrack=%.12g Jflow=%.12g Jsafety=%.12g\n', ...
        A.J_tracking,A.J_flow,A.J_safety);
    fprintf('M1 pilot: Lflow=%.12g U=%.12g gap=%.12g feasible=%d\n', ...
        A.lower_bound,A.upper_bound,A.gap,A.feasible);
    save(fullfile('results','m1_pilot_audit.mat'),'A');
end
