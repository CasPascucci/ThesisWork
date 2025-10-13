function [optParams, optCost] = optimizationLoop(paramsX0, problemParams, nonDimParams, refVals)

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
              0  0 -1];     % tgo >= 0.001
    bineq = [-1e-6; -4; -0.01];

    %Specific Test
        %lb = [1.0, 6, 9.70848233443312];%9.70848233443312
        %ub = [1.0, 6, 9.70848233443312];
    % E-Guidance Test
        lb = [0, 0, 3];
        ub = [1, 18, 11];

    fminconOptions = optimoptions('fmincon', 'Display', 'iter-detailed', 'MaxFunctionEvaluations', 5000, ...
        'FiniteDifferenceType','central','FiniteDifferenceStepSize', 1e-6, ...
        'Algorithm','sqp','OptimalityTolerance', 1e-6, 'EnableFeasibilityMode',true);


    obj = @(params) objectiveFunction(params, afStar, rfStar, r0, vfStar, v0, gConst, nonDimParams);
    nonlincon = @(params) thrustLimits(params, r0, v0, rfStar, vfStar, afStar, gConst, isp, minThrust, maxThrust);

    [optParams, optCost] = fmincon(obj, paramsX0, Aineq, bineq, [], [], lb, ub, nonlincon, fminconOptions);
    
end
%% Functions
function cost = objectiveFunction(params, afStar, rfStar, r, vfStar, v, gConst, nonDimParams)
    gamma  = params(1);
    kr     = params(2);
    tgo   = params(3);
    

    gamma1 = gamma;
    gamma2 = kr/(gamma+2) - 2;


    [c1, c2] = calculateCoeffs(r, v, tgo, gamma1, gamma2, afStar, rfStar, vfStar, gConst);

    % Simpson Composite 1/3 Rule
    tspan = linspace(0,tgo,997);

    aT = afStar + c1*tspan.^gamma1 + c2*tspan.^gamma2;
    %simpson = simpsonComp13Integral(tspan,dot(aT,aT));

    %% Fuel Opt
    m0 = 1;
    aTmag = vecnorm(aT,2,1);
    isp = nonDimParams.ispND;

    simpson = simpsonComp13Integral(tspan,aTmag);

    mf = m0 * exp(-simpson/isp);
    cost = m0-mf;

end

function [c, ceq] = thrustLimits(params, r0, v0, rfStar, vfStar, afStar, gConst, isp, minThrust, maxThrust)
    gamma  = params(1);
    kr     = params(2);
    tgo   = params(3);
    m0 = 1;

    gamma1 = gamma;
    gamma2 = kr/(gamma+2) - 2;


    [c1, c2] = calculateCoeffs(r0, v0, tgo, gamma1, gamma2, afStar, rfStar, vfStar, gConst);

    tspan = linspace(0,tgo,997);

    aT = afStar + c1*tspan.^gamma1 + c2*tspan.^gamma2;
    aTmag = vecnorm(aT,2,1);

    Q = cumtrapz(tspan,aTmag./isp);

    m = m0 .* exp(-Q);
    m = flip(m);

    thrust = m .*aTmag;

    upper = thrust - maxThrust;
    lower = minThrust - thrust;

    c = [upper(:), lower(:)];
    ceq = [];
end