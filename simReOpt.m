function [tTraj, stateTraj, aTList, flag_thrustGotLimited, optHistory, exitFlags] = simReOpt(gamma0,kr0,tgo0, problemParams, nonDimParams, refVals, delta_t, optimizationParams, betaParam, verboseOutput)

% Original IC's
r0 = nonDimParams.r0ND;
v0 = nonDimParams.v0ND;
m0 = nonDimParams.m0ND;
rfStar = nonDimParams.rfStarND;
vfStar = nonDimParams.vfStarND;
afStar = nonDimParams.afStarND;
isp = nonDimParams.ispND;
rMoonND = nonDimParams.rMoonND;
gConst = nonDimParams.gConst;

X0 = [r0; v0; m0];
t_curr = 0;
odeoptions = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

% Setup Results
gamma = gamma0;
kr = kr0;
tgo = tgo0;

tTraj = [];
stateTraj = [];
aTList = [];
flag_thrustGotLimited = false;
optHistory = [t_curr, gamma, kr, tgo];
exitFlags = [];

% Main loop
while (tgo * refVals.T_ref) > optimizationParams.updateStop
    t_end = t_curr + (optimizationParams.updateFreq/refVals.T_ref);
    [tSeg, stateSeg] = ode45(@(t, X) trajectorySegment(t, X, gamma, kr, tgo, isp, rMoonND, rfStar, vfStar, afStar, gConst, nonDimParams.minThrustND, nonDimParams.maxThrustND), [t_curr, t_end], X0, odeoptions);
    tTraj = [tTraj; tSeg];
    stateTraj = [stateTraj; stateSeg];

    segAT = zeros(3, numel(tSeg));
    for i = 1:numel(tSeg)
        r = stateSeg(i,1:3)';
        v = stateSeg(i,4:6)';
        m = stateSeg(i,7);
        tgo_now = tgo - (tSeg(i) - t_curr);
        [aTi, limited] = computeAccel(r, v, m, gamma, kr, tgo_now, rfStar, vfStar, afStar, gConst, nonDimParams);
        segAT(:,i) = aTi;
        flag_thrustGotLimited = flag_thrustGotLimited || limited;
    end
    aTList = [aTList, segAT];
    
    % Update state conditions for reOpt
    X0 = stateSeg(end,:)';
    t_curr = tSeg(end);

    newNonDimParams = nonDimParams;
    newNonDimParams.r0ND = X0(1:3);
    newNonDimParams.v0ND = X0(4:6);
    newNonDimParams.m0ND = X0(7);

    % ReOpt
    paramsX0 = [gamma, kr, tgo];
    [optParams, optCost, aTOptim, mOptim, rdOptim, vdOptim, exitflag] = optimizationLoop(paramsX0, betaParam, problemParams, newNonDimParams, optimizationParams, refVals, delta_t, verboseOutput, false);
    gamma = optParams(1);
    kr = optParams(2);
    tgo = optParams(3);

    optHistory = [optHistory; [t_curr, gamma, kr, tgo]];
    exitFlags = [exitFlags; exitflag];

    if (tgo*refVals.T_ref) <= optimizationParams.updateStop
        break;
    end
end
end


function [aTi, limited] = computeAccel(r, v, m, gamma, kr, tgo, rfStar, vfStar, afStar, gGuidance, nonDimParams)
    minAccel = nonDimParams.minThrustND / m;
    maxAccel = nonDimParams.maxThrustND / m;

    aT1 = gamma*(kr/(2*gamma +4) -1)*afStar;
    aT2 = (gamma*kr/(2*gamma+4)-gamma-1)*gGuidance;
    aT3 = ((gamma+1)/tgo)*(1-kr/(gamma+2))*(vfStar-v);
    aT4 = (kr/tgo^2)*(rfStar-r-v*tgo);
    aTi = aT1 + aT2 + aT3 + aT4;
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

function dXdt = trajectorySegment(t, X, gamma, kr, tgo0, isp, rMoonND, rfStar, vfStar, afStar, gConst, minThrust, maxThrust)
    r = X(1:3);
    v = X(4:6);
    m = X(7);
    tgo = tgo0 - t;

    g = -(rMoonND^2) * r / (norm(r)^3);
    [aT,~] = computeAccel(r, v, m, gamma, kr, tgo, rfStar, vfStar, afStar, gConst, struct('minThrustND', minThrust, 'maxThrustND', maxThrust));
    F_mag = norm(aT) * m;
    dm_dt = -F_mag / isp;
    dXdt = [v; aT + g; dm_dt];
end