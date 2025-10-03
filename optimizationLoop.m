function [optParams, optCost] = optimizationLoop(paramsX0, problemParams, nonDimParams, refVals)
    gammaTest = paramsX0(1);
    krTest = paramsX0(2);
    tgoTest = paramsX0(3);

    r0 = nonDimParams.r0ND;
    v0 = nonDimParams.v0ND;

    rfStar = nonDimParams.rfStarND;
    vfStar = nonDimParams.vfStarND;
    afStar = nonDimParams.afStarND;

    gConst = nonDimParams.gConst;

    % Fmincon Constraints
    Aineq = [-1  0  0;      % -gamma <= -1e-6 ---- gamma >= 1e-6
              2 -1  0;      % 2gamma - kr <= -4-1e-4 ---- kr >= 2gamma +4+1e-4 ---- kr >= 2*(gamma + 2) + 1e-4
              0  0 -1];     % tgo >= 0.001
    bineq = [-1e-3; -4-1e-4; -0.01];
    lb = [1.0, 6, 8];%9.70249498737309
    ub = [1.0, 6.01, 11];

    fminconOptions = optimoptions('fmincon', 'Display', 'iter-detailed', 'MaxFunctionEvaluations', 5000, ...
        'FiniteDifferenceType', 'central', 'FiniteDifferenceStepSize', 1e-6, ...
        'Algorithm', 'sqp');%, 'HessianApproximation', 'lbfgs');


    obj = @(params) objectiveFunction(params, afStar, rfStar, r0, vfStar, v0, gConst, problemParams.landingLatDeg, problemParams.landingLonDeg);
    [optParams, optCost] = fmincon(obj, paramsX0, Aineq, bineq, [], [], lb, ub, [], fminconOptions);
    
end
%% Functions
function cost = objectiveFunction(params, afStar, rfStar, r, vfStar, v, gConst, landingLat, landingLon)
    gamma  = params(1);
    kr     = params(2);
    tgo   = params(3);

    gamma1 = gamma;
    gamma2 = kr/(gamma+2) - 2;

    rENU = MCMF2ENU(r,landingLat,landingLon,false);
    vENU = MCMF2ENU(v,landingLat,landingLon,true);
    rfStarENU = MCMF2ENU(rfStar,landingLat,landingLon,false);
    vfStarENU = MCMF2ENU(vfStar,landingLat,landingLon,true);
    gConstENU = MCMF2ENU(gConst,landingLat, landingLon, true); % This only applies as long as landing at south pole
    afStarENU = MCMF2ENU(afStar, landingLat, landingLon, true); %This only applies as long as landing at south pole

    [c1, c2] = calculateCoeffs(rENU, vENU, tgo, gamma1, gamma2, afStarENU, rfStarENU, vfStarENU, gConstENU);

    % Simpson Composite 1/3 Rule
    tspan = linspace(0,tgo,997);
    aT = afStarENU + c1*tspan.^gamma1 + c2*tspan.^gamma2;
    simpson = simpsonComp13Integral(tspan,dot(aT,aT));
    
    cost = simpson;

end