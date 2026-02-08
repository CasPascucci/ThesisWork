function efun = FP2DG_events(primal)
    r0 = primal.initial.states(1:3);
    v0 = primal.initial.states(4:6);
    m0 = primal.initial.states(7);

    rf = primal.final.states(1:3);
    vf = primal.final.states(4:6);

    efun = [r0; v0; m0; rf; vf];
end
