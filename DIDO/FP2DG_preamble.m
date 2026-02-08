function [r, v, m, aT, t] = FP2DG_preamble(primal)

r = primal.states(1:3, :);
v = primal.states(4:6, :);
m = primal.states(7, :);

aT = primal.controls(1:3, :);

t = primal.time;


end