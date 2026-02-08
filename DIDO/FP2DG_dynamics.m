function dxdt = FP2DG_dynamics(primal)
[r, v, m, aT, ~] = FP2DG_preamble(primal);

rMoonND = primal.constants.nonDim.rMoonND;
ispND = primal.constants.nonDim.ispND;

rMag = vecnorm(r);
g = -(rMoonND^2) .* r ./ (rMag.^3);

drdt = v;
dvdt = aT + g;

aTmag = vecnorm(aT);
dmdt = -(aTmag .* m) ./ispND;

dxdt = [drdt; dvdt; dmdt];
end