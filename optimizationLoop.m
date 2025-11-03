function [optParams, optCost, aTOptim, mOptim, rdOptim, vdOptim, exitflag] = optimizationLoop(paramsX0, betaParam, problemParams, nonDimParams, optimParams, refVals, delta_t, verboseOutput, dispersion)
    if nargin < 9
        dispersion = false;
    end
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

    lb = [1e-6, 0, 0.01];
    ub = [8, 30, 11];

    if optimParams.pointingEnabled
        if dispersion
            fminconOptions = optimoptions('fmincon', 'Display','none', 'MaxFunctionEvaluations', 10000, ...
            'FiniteDifferenceType','central','FiniteDifferenceStepSize', 1e-4,'MaxIterations', 1000, ...
            'Algorithm','active-set', 'EnableFeasibilityMode',true, ...
            'HessianApproximation','lbfgs');
        else
            fminconOptions = optimoptions('fmincon', 'Display','final', 'MaxFunctionEvaluations', 10000, ...
            'FiniteDifferenceType','central','FiniteDifferenceStepSize', 1e-4,'MaxIterations', 1000, ...
            'Algorithm','active-set', 'EnableFeasibilityMode',true, ...
            'HessianApproximation','lbfgs');
        end
    else
        if dispersion
            fminconOptions = optimoptions('fmincon', 'Display', 'none', 'MaxFunctionEvaluations', 10000, ...
            'FiniteDifferenceType','central','FiniteDifferenceStepSize', 1e-4,'MaxIterations', 1000, ...
            'Algorithm','interior-point', 'EnableFeasibilityMode',true, ...
            'HessianApproximation','lbfgs');
        else
            fminconOptions = optimoptions('fmincon', 'Display', 'final', 'MaxFunctionEvaluations', 10000, ...
            'FiniteDifferenceType','central','FiniteDifferenceStepSize', 1e-4,'MaxIterations', 1000, ...
            'Algorithm','interior-point', 'EnableFeasibilityMode',true, ...
            'HessianApproximation','lbfgs');
        end

    end
    nodeCount = optimParams.nodeCount;


    %[rfStar, vfStar, afStar, tgoVirt] = computeBeyondTerminationTargeting(r0, v0, paramsX0(1), paramsX0(2), rfStar, vfStar, afStar, delta_t, paramsX0(3), gConst);
    %paramsX0(3) = tgoVirt

    obj = @(params) objectiveFunction(params, betaParam, afStar, rfStar, r0, vfStar, v0, gConst, nonDimParams, optimParams);
    nonlincon = @(params) nonLinearLimits(params, r0, v0, rfStar, vfStar, afStar, gConst, isp, minThrust, maxThrust, optimParams, problemParams, refVals);

    
    [optParams, optCost, exitflag, output] = fmincon(obj, paramsX0, Aineq, bineq, [], [], lb, ub, nonlincon, fminconOptions);

    if verboseOutput
        fprintf('=== Optimization Results ===\n');
        fprintf('Exit flag: %d\n', exitflag);
        fprintf('Exit message: %s\n', output.message);
        fprintf('Final cost: %.6e\n', optCost);
        fprintf('Iterations: %d\n', output.iterations);
        fprintf('First-order optimality: %.6e (tolerance: %.1e)\n', output.firstorderopt, fminconOptions.OptimalityTolerance);
        fprintf('\nOptimal parameters:\n');
        fprintf('  gamma = %.6f\n', optParams(1));
        fprintf('  kr    = %.6f\n', optParams(2));
        fprintf('  tgo   = %.6f (s)\n', optParams(3)* refVals.T_ref);
    
        [c, ~] = nonlincon(optParams);
        activeIneq = find(c >= 0); % Nearly active
        fprintf('\nActive inequality constraints: %d/%d\n', length(activeIneq), length(c));
        [maxViolation, idx] = max(c);
        fprintf('Max violation: %.3f at constraint %d\n', maxViolation, idx);
    end

    


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
function cost = objectiveFunction(params, betaParam, afStar, rfStar, r, vfStar, v, gConst, nonDimParams, optimParams)
    gamma  = params(1);
    kr     = params(2);
    tgo   = params(3);
    nodeCount = optimParams.nodeCount;

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

function [c, ceq] = nonLinearLimits(params, r0, v0, rfStar, vfStar, afStar, gConst, isp, minThrust, maxThrust, optimParams, problemParams, refVals)
    nodeCount = optimParams.nodeCount;
    glideSlopeFlag = optimParams.glideSlopeEnabled;
    %pointingFlag = optimParams.pointingEnabled;
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

    c = [upper(:); lower(:)];
    ceq = [];
    if glideSlopeFlag
        glideNodes = floor(nodeCount*0.20);
        freeGlideNodes = floor(nodeCount*0.01) + 1;
        tgospanGlide = tgospan(1:glideNodes);
        phi1hat = (tgospanGlide.^(gamma1+2))./((gamma1+1)*(gamma1+2));
        phi2hat = (tgospanGlide.^(gamma2+2))./((gamma2+1)*(gamma2+2));
        rdOptim = rfStar + c1*phi1hat + c2*phi2hat - vfStar.*tgospanGlide + 0.5*(gConst+afStar).*tgospanGlide.^2;
        rdOptimTOPO = MCMF2ENU(rdOptim,problemParams.landingLatDeg,problemParams.landingLonDeg,true,false);
        rUnitVec = rdOptimTOPO./vecnorm(rdOptimTOPO,2,1);


        if norm(rdOptimTOPO(:,1)) < 1e-10
            rUnitVec(:,1) = [0;0;1];
        end
        vertUnitVec = [0;0;1];

        rMoon = problemParams.rMoon;
        rdOptimDim = rdOptim*refVals.L_ref;
        alt_opt = vecnorm(rdOptimDim,2,1) - rMoon;
        frac = (alt_opt - 250)/250;
        frac(frac < 0) = 0;
        frac(frac > 1) = 1;
        theta = (frac*45) + 45;
        theta(1:freeGlideNodes) = 90;
        cGlide = zeros(glideNodes,1);
        for idx = 1:glideNodes
            cGlide(idx) = cosd(theta(idx)) - dot(rUnitVec(:,idx),vertUnitVec);
        end
        c = [c; cGlide];
    end
    % if glideSlopeFlag
    %     glideNodes = max(1, floor(0.20*nodeCount));   % first 20% of indices (end of flight)
    %     freeGlideNodes = max(1, floor(0.02*nodeCount));   % first 2% of indices
    %     idx = 1:glideNodes;
    % 
    %     % calc TOPO radius, tgospan works from landing backwards
    %     tgospanGlide = tgospan(idx);
    %     phi1hat = (tgospanGlide.^(gamma1+2))./((gamma1+1)*(gamma1+2));
    %     phi2hat = (tgospanGlide.^(gamma2+2))./((gamma2+1)*(gamma2+2));
    %     rdOptim = rfStar + c1*phi1hat + c2*phi2hat - vfStar.*tgospanGlide + 0.5*(gConst+afStar).*tgospanGlide.^2;
    %     rdOptimTOPO = MCMF2ENU(rdOptim, problemParams.landingLatDeg, problemParams.landingLonDeg, true, false);
    % 
    %     % unit vectors
    %     rdNorms = vecnorm(rdOptimTOPO,2,1);
    %     rUnitVec = rdOptimTOPO ./ max(rdNorms, 1e-12);
    %     vertUnitVec = [0;0;1];
    % 
    %     % calc altitude
    %     rMoon     = problemParams.rMoon;
    %     rdDim     = rdOptim * refVals.L_ref;
    %     alt_opt   = vecnorm(rdDim,2,1) - rMoon;
    % 
    %     %fraction for how far through linear theta transition
    %     frac = (alt_opt - 250)/250;     % 0 below 250 m, 1 at 500 m and above
    %     frac(frac < 0) = 0;
    %     frac(frac > 1) = 1;
    %     theta = (frac*45) + 45;         % degrees from vertical
    % 
    %     % first 2% of indices are the very end of flight -> set to 90 deg
    %     theta(1:freeGlideNodes) = 90;
    % 
    %     % constraint: angle(rhat, zhat) <= theta  <=>  rhat·zhat >= cos(theta)
    %     % so cGlide = cos(theta) - (rhat·zhat) <= 0
    %     cGlide = cosd(theta) - sum(rUnitVec .* vertUnitVec, 1);
    % 
    %     % append
    %     c = [c; cGlide(:)];
    % end

    % if pointingFlag
    %     aTTOPO = MCMF2ENU(aT,problemParams.landingLatDeg,problemParams.landingLonDeg,false,false);
    %     aTmagTOPO = vecnorm(aTTOPO,2,1);
    %     thrustUnitVec = aTTOPO./aTmagTOPO;
    %     vertUnitVec = [0;0;1];
    %     dotProducts = thrustUnitVec'*vertUnitVec;
    %     dotProducts = min(1,max(-1,dotProducts));
    %     phi = acosd(dotProducts); % Values at idx 1 are landing, idx end are PDI, Deg
    %     phiR = (optimParams.maxTiltRate * refVals.T_ref)* tgospan; % DEG
    %     phiA = 0.5 * (optimParams.maxTiltAccel * refVals.T_ref^2) * (tgospan.^2); % DEG
    %     phiLimit = min(min(phiR,phiA),90)';
    %     %phiLimit = max(phiLimit,0.1); % Adds a 0.1 degree floor, useful for final nodes to touchdown when constraint is basically 0
    % 
    %     cPointing = phi - (phiLimit); % phi - phiLim <= 0
    %     c(3*nodeCount+4:4*nodeCount) = cPointing(4:end); % Currently not imposing in last few nodes before landing
    % end
end