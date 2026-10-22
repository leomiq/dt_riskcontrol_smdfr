function lb = m1_toy_det_bound(box)
% Analytic enclosure for both topologies:
% A (u<=0.5): T=[u,1-u]; B (u>0.5): T=[u,1.2-u].
    if ~isequal(size(box),[1 2]) || any(~isfinite(box)) || ...
            box(1)<0 || box(2)>1 || box(1)>box(2)
        error('Toy box outside [0,1].');
    end
    candidates = [];
    if box(1)<=0.5
        lo=box(1); hi=min(box(2),0.5);
        T_box=[lo hi;1-hi 1-lo];
        candidates(end+1)=certified_cost_lower_bound_box(T_box,[0 0],0); %#ok<AGROW>
    end
    if box(2)>0.5
        lo=max(box(1),0.5); hi=box(2);
        T_box=[lo hi;1.2-hi 1.2-lo];
        candidates(end+1)=certified_cost_lower_bound_box(T_box,[0 0],0); %#ok<AGROW>
    end
    lb=min(candidates);
end
