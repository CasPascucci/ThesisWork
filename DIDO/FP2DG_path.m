function h = FP2DG_path(primal)
    [~, ~, m, aT, ~] = FP2DG_preamble(primal);

    h = m .* vecnorm(aT); % thrust magnitude
end
