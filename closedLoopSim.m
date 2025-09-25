function [tTraj, stateTraj] = closedLoopSim(gamma,kr,tgo0, problemParams, nonDimParams, refVals)

    r0 = nonDimParams.r0ND;
    v0 = nonDimParams.v0ND;
    m0 = nonDimParams.m0ND;
    
    rMoonND = nonDimParams.rMoonND;
    isp = nonDimParams.ispND;

    rfStar = nonDimParams.rfStarND;
    vfStar = nonDimParams.vfStarND;
    afStar = nonDimParams.afStarND;

    X0 = [r0; v0; m0];
    odeoptions = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
     [tTraj, stateTraj] = ode45(@(t, X) trajectory(t, X, gamma, kr, ...
                          tgo0, isp, rMoonND, rfStar, vfStar, afStar), [0, tgo0], X0, odeoptions);

     rState = stateTraj(:,1:3);
     vState = stateTraj(:,4:6);
     mState = stateTraj(:,7);
end



function dXdt = trajectory(t, X, gamma, kr, tgo0, isp, rMoonND, rfStar, vfStar, afStar)
    r    = X(1:3);
    v    = X(4:6); 
    mass = X(7);

    g = -(rMoonND^2) * r / (norm(r)^3);
    

    tgo  = tgo0 - t;
    aT1 = afStar + (( (gamma*kr) / (2*(gamma + 2)) ) - gamma - 1)*(afStar + g);
    aT2 = ((gamma + 1) / tgo)*(1 - kr/(gamma + 2))*(vfStar - v);
    aT3 = (rfStar - r - v*tgo)* kr / tgo^2;
    aT = aT1 + aT2 + aT3;
    F_mag = norm(aT) * mass;

    dm_dt = -F_mag / (isp);
    dXdt = [v; aT + g; dm_dt];
end