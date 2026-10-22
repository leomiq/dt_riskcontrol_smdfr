function test_certified_cost_lower_bound_box()
% Synthetic checks for the algebraic bound; no Simulink licence required.
    target = 1179.2712;
    T = [1178.9; 1180.1; 1245.0];
    m = [2; 8; 15];
    point_T = [T T];
    point_m = [m m];
    [lb, terms] = certified_cost_lower_bound_box(point_T,point_m,target);
    exact = sum((T-target).^2) + 1e-4*sum(max(0,abs(diff(m))-5).^2) ...
        + 10*sum(max(0,T-1244).^2);
    assert(abs(lb-exact) <= 1e-10*max(1,exact));
    assert(strcmp(terms.safety_mode,'absolute'));

    opts = struct('dsm_thr',0.01,'lambda_dsm',1e4);
    [lb_dsm, terms_dsm] = certified_cost_lower_bound_box(point_T,point_m,target,opts);
    exact_dsm = sum((T-target).^2) + 1e-4*sum(max(0,abs(diff(m))-5).^2) ...
        + opts.lambda_dsm*sum(max(0,(T-target)/(1244-target)-opts.dsm_thr).^2);
    assert(abs(lb_dsm-exact_dsm) <= 1e-10*max(1,exact_dsm));
    assert(strcmp(terms_dsm.safety_mode,'dsm'));

    % The box covers the point, so its lower bound cannot exceed point cost.
    [lb_box, ~] = certified_cost_lower_bound_box([T-2,T+3],[m-1,m+2],target,opts);
    assert(lb_box <= exact_dsm + 1e-10);

    % A box wholly above target gives a strictly positive tracking bound.
    [lb_positive,~] = certified_cost_lower_bound_box( ...
        [target+1,target+2],[0,0;10,10],target);
    assert(lb_positive >= 1 + 1e-4*25 - 1e-12);

    did_fail = false;
    try
        certified_cost_lower_bound_box([NaN 1],[0 1],target);
    catch
        did_fail = true;
    end
    assert(did_fail);
    fprintf('PASS: exact singleton costs, enclosing-box bound, positive bound, invalid input.\n');
end
