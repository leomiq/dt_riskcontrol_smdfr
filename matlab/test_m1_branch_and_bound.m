function test_m1_branch_and_bound()
% Exact toy optimum J*=0.5 at u=0.5.
    opts=struct('abs_gap',1e-3,'max_splits',2000);
    R=m1_branch_and_bound([0 1],@m1_toy_det_bound,@m1_toy_det_cost,opts);
    assert(R.tolerance_met);
    assert(R.lower_bound<=0.5+1e-12);
    assert(R.upper_bound>=0.5-1e-12);
    assert(abs(R.upper_bound-0.5)<=1e-12);
    fprintf('M1 toy: L=%.9f U=%.9f gap=%.9f splits=%d evaluations=%d\n', ...
        R.lower_bound,R.upper_bound,R.gap,R.splits,R.evaluations);
    bad=false;
    try
        m1_branch_and_bound([0 1],@(box) -1,@m1_toy_det_cost,opts);
    catch
        bad=true;
    end
    assert(bad);
    fprintf('PASS: hybrid topology, global gap, invalid bound.\n');
end
