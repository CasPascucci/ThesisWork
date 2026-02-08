function [eCost, fCost] = FP2DG_cost(primal)
    [~, ~, ~, aT, ~] = FP2DG_preamble(primal);
    beta = primal.constants.beta;

    eCost = 0;
    fCost = beta * vecnorm(aT) + (1 - beta) .* dot(aT, aT);
end
