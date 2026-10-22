function [J,feasible] = m1_toy_det_cost(u)
    if numel(u)~=1 || ~isfinite(u) || u<0 || u>1
        error('Toy control outside [0,1].');
    end
    if u<=0.5, T=[u;1-u]; else, T=[u;1.2-u]; end
    J=sum(T.^2); feasible=true;
end
