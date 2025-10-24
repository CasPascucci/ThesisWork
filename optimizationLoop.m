function [optParams, optCost, aTOptim, mOptim, rdOptim, vdOptim] = optimizationLoop(paramsX0, betaParam, problemParams, nonDimParams, optimParams, refVals, delta_t)

    r0 = nonDimParams.r0ND;
    v0 = nonDimParams.v0ND;

    rfStar = nonDimParams.rfStarND;
    vfStar = nonDimParams.vfStarND;
    afStar = nonDimParams.afStarND;

    gConst = nonDimParams.gConst;

    maxThrust = nonDimParams.maxThrustND;
    minThrust = nonDimParams.minThrustND;
    isp = nonDimParams.ispND;

    % Fmincon Constraints
    Aineq = [-1  0  0;      % -gamma <= -1e-6 ---- gamma >= 1e-6
              2 -1  0;      % 2gamma - kr <= -4-1e-4 ---- kr >= 2gamma +4+1e-4 ---- kr >= 2*(gamma + 2) + 1e-4
              0  0 -1];     % tgo >= 0.01
    bineq = [-1e-6; -4-1e-4; -0.01];

    lb = [0, 0, 3];
    ub = [8, 30, 11];

    fminconOptions = optimoptions('fmincon', 'Display', 'iter-detailed', 'MaxFunctionEvaluations', 10000, ...
        'FiniteDifferenceType','central','FiniteDifferenceStepSize', 1e-6,'MaxIterations', 1000, ...
        'Algorithm','interior-point','OptimalityTolerance', 1e-8, 'EnableFeasibilityMode',true, ...
        'ConstraintTolerance',0);
    nodeCount = optimParams.nodeCount;
    glideSlopeFlag = optimParams.glideSlopeEnabled;
    pointingFlag = optimParams.pointingEnabled;

    %[rfStar, vfStar, afStar, tgoVirt] = computeBeyondTerminationTargeting(r0, v0, paramsX0(1), paramsX0(2), rfStar, vfStar, afStar, delta_t, paramsX0(3), gConst);
    %paramsX0(3) = tgoVirt

    obj = @(params) objectiveFunction(params, betaParam, afStar, rfStar, r0, vfStar, v0, gConst, nonDimParams, nodeCount);
    nonlincon = @(params) nonLinearLimits(params, r0, v0, rfStar, vfStar, afStar, gConst, isp, minThrust, maxThrust, nodeCount, glideSlopeFlag, pointingFlag, problemParams, refVals);

    [optParams, optCost] = fmincon(obj, paramsX0, Aineq, bineq, [], [], lb, ub, nonlincon, fminconOptions);

    gamma1 = optParams(1);
    gamma2 = optParams(2)/(optParams(1)+2) - 2;

    [c1, c2] = calculateCoeffs(r0, v0, optParams(3), gamma1, gamma2, afStar, rfStar, vfStar, gConst);

    tgospan = linspace(0,optParams(3),nodeCount);
    %tspan = optParams(3) - tgospan;

    aTOptim = afStar + c1*tgospan.^gamma1 + c2*tgospan.^gamma2;

    phi1hat = (tgospan.^(gamma1+2))./((gamma1+1)*(gamma1+2));
    phi2hat = (tgospan.^(gamma2+2))./((gamma2+1)*(gamma2+2));
    phi1bar = (tgospan.^(gamma1+1))./(gamma1+1);
    phi2bar = (tgospan.^(gamma2+1))./(gamma2+1);

    rdOptim = rfStar + c1*phi1hat + c2*phi2hat - vfStar.*tgospan + 0.5*(gConst+afStar).*tgospan.^2;
    vdOptim = vfStar + c1*phi1bar + c2*phi2bar -(gConst+afStar).*tgospan;

    aTNorm = vecnorm(aTOptim,2,1);

    Q = cumtrapz(tgospan,aTNorm./isp);
    Q = Q(end) - Q;

    mOptim = 1 .* exp(-Q);

    
end
%% Functions
function cost = objectiveFunction(params, betaParam, afStar, rfStar, r, vfStar, v, gConst, nonDimParams, nodeCount)
    gamma  = params(1);
    kr     = params(2);
    tgo   = params(3);
    

    gamma1 = gamma;
    gamma2 = kr/(gamma+2) - 2;

    
    [c1, c2] = calculateCoeffs(r, v, tgo, gamma1, gamma2, afStar, rfStar, vfStar, gConst);

    tspan = linspace(0,tgo,nodeCount);

    aT = afStar + c1*tspan.^gamma1 + c2*tspan.^gamma2;
    aTmag = vecnorm(aT,2,1);

    simpson1 = simpsonComp13Integral(tspan,aTmag);
    simpson2 = simpsonComp13Integral(tspan,dot(aT,aT));

    cost = betaParam*simpson1 + (1-betaParam)*simpson2;

end

function [c, ceq] = nonLinearLimits(params, r0, v0, rfStar, vfStar, afStar, gConst, isp, minThrust, maxThrust, nodeCount, glideSlopeFlag, pointingFlag, problemParams, refVals)
    if pointingFlag || glideSlopeFlag
        c = zeros(4*nodeCount,1)-1;
    end
    gamma  = params(1);
    kr     = params(2);
    tgo0   = params(3);
    m0 = 1;

    gamma1 = gamma;
    gamma2 = kr/(gamma+2) - 2;



    [c1, c2] = calculateCoeffs(r0, v0, tgo0, gamma1, gamma2, afStar, rfStar, vfStar, gConst);
    tgospan = linspace(0,tgo0,nodeCount);

    aT = afStar + c1*tgospan.^gamma1 + c2*tgospan.^gamma2;
    aTmag = vecnorm(aT,2,1);

    Q = cumtrapz(tgospan,aTmag./isp);
    Q = Q(end) - Q;

    m = m0 .* exp(-Q);

    thrust = m .*(aTmag);

    upper = thrust - (maxThrust);
    lower = minThrust - thrust;
    if pointingFlag || glideSlopeFlag
        c(1:2*nodeCount) = [upper(:); lower(:)];
    else
        c = [upper(:); lower(:)];
    end
    ceq = [];
    if glideSlopeFlag
        phi1hat = (tgospan.^(gamma1+2))./((gamma1+1)*(gamma1+2));
        phi2hat = (tgospan.^(gamma2+2))./((gamma2+1)*(gamma2+2));
        rdOptim = rfStar + c1*phi1hat + c2*phi2hat - vfStar.*tgospan + 0.5*(gConst+afStar).*tgospan.^2;
        rdOptimTOPO = MCMF2ENU(rdOptim,problemParams.landingLatDeg,problemParams.landingLonDeg,true,false);
        rUnitVec = rdOptimTOPO./vecnorm(rdOptimTOPO,2,1);
        if norm(rdOptimTOPO(:,1)) < 1e-10
            rUnitVec(:,1) = [0;0;1];
        end
        vertUnitVec = [0;0;1];
        
        rMoon = problemParams.rMoon;
        rdOptimDim = rdOptim*refVals.L_ref;
        alt_opt = vecnorm(rdOptimDim,2,1) - rMoon;

        altMask = find(alt_opt <= 500,1,'last');
        altActive = alt_opt(1:altMask);
        frac = (altActive - 250)/250;
        frac(frac < 0) = 0;
        theta = (frac*45) + 45;
        cGlide = zeros(altMask,1);
        for idx = 1:altMask
            cGlide(idx) = cosd(theta(idx)) - dot(rUnitVec(:,idx),vertUnitVec);
        end
        c(2*nodeCount+1:2*nodeCount+altMask) = cGlide;
    end
end