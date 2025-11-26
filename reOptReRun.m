function outputSingle = reOptReRun(idx, ICStates, optHistory, betaParam, problemParams, nonDimParams, optimizationParams, refVals, delta_tND, verboseOutput)
       ICStates = table2array(ICStates);
       optHistory = table2array(optHistory);

       IC = ICStates(idx,:);
       nonDimParams.r0ND = IC(1:3)';
       nonDimParams.v0ND = IC(4:6)';
       nonDimParams.m0ND = IC(7);
       X0 = [nonDimParams.r0ND; nonDimParams.v0ND; nonDimParams.m0ND];
       paramsX0 = optHistory(idx,[2,3,5]);
       gamma = paramsX0(1);
       kr = optHistory(idx,4);
       tgo = optHistory(idx,5);
       updateFreqND = optimizationParams.updateFreq / refVals.T_ref;
       tStart = optHistory(idx,1);
       tEnd = tStart + updateFreqND;

       rfStar = nonDimParams.rfStarND;
       vfStar = nonDimParams.vfStarND;
       afStar = nonDimParams.afStarND;
       isp = nonDimParams.ispND;
       rMoonND = nonDimParams.rMoonND;
       gGuidance = nonDimParams.gConst;

       [optParams, optCost, aTOptim, mOptim, rdOptim, vdOptim, exitflag] = ...
        optimizationLoop(paramsX0, betaParam, problemParams, nonDimParams, optimizationParams, refVals, delta_tND, verboseOutput, false);
        

       odeoptions = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
       [tSeg, stateSeg] = ode45(@(t, X) trajectorySegmentODE(t, X, tStart, gamma, kr, tgo, ...
            isp, rMoonND, rfStar, vfStar, afStar, gGuidance, nonDimParams), ...
            [tStart, tEnd], X0, odeoptions);

       % Do optimization and sim from here
       % Issue might be thrust getting pushed up in optimization, but it
       % might be outside of that segment's own time

end

function [aTi, limited] = computeAccel(r, v, m, gamma, kr, tgo, rfStar, vfStar, afStar, gGuidance, nonDimParams)
    % Compute thrust acceleration with limits
    minAccel = nonDimParams.minThrustND / m;
    maxAccel = nonDimParams.maxThrustND / m;
    
    % Guidance law
    aT1 = gamma * (kr/(2*gamma + 4) - 1) * afStar;
    aT2 = (gamma*kr/(2*gamma+4) - gamma - 1) * gGuidance;
    aT3 = ((gamma+1)/tgo) * (1 - kr/(gamma+2)) * (vfStar - v);
    aT4 = (kr/tgo^2) * (rfStar - r - v*tgo);
    aTi = aT1 + aT2 + aT3 + aT4;
    
    % Apply thrust limits
    normAT = norm(aTi);
    if normAT > maxAccel
        aTi = aTi / normAT * maxAccel;
        limited = true;
    elseif normAT < minAccel
        aTi = aTi / normAT * minAccel;
        limited = true;
    else
        limited = false;
    end
end

function dXdt = trajectorySegmentODE(t, X, t_reopt_start, gamma, kr, tgo_at_reopt, isp, rMoonND, rfStar, vfStar, afStar, gGuidance, nonDimParams)
    % ODE function for trajectory segment
    % t - current time (elapsed from simulation start)
    % t_reopt_start - time when current parameters were set
    % tgo_at_reopt - time-to-go when current parameters were set
    
    r = X(1:3);
    v = X(4:6);
    m = X(7);
    
    % Calculate current time-to-go
    t_since_reopt = t - t_reopt_start;
    tgo = tgo_at_reopt - t_since_reopt;
    
    % Ensure tgo doesn't go negative
    if tgo < 0.001
        tgo = 0.001;
    end
    
    % Gravitational acceleration
    g = -(rMoonND^2) * r / (norm(r)^3);
    
    % Compute thrust acceleration
    [aT, ~] = computeAccel(r, v, m, gamma, kr, tgo, rfStar, vfStar, afStar, gGuidance, nonDimParams);
    
    % Thrust magnitude and mass flow rate
    F_mag = norm(aT) * m;
    dm_dt = -F_mag / isp;
    
    % State derivatives
    dXdt = [v; aT + g; dm_dt];
end