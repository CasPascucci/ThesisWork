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
    bineq = [-1e-6; -4-1e-4; -0.001];
    lb = [1, 6, 3];
    ub = [1, 25, 20];

    fminconOptions = optimoptions('fmincon', 'Display', 'none', 'MaxFunctionEvaluations', 1000, ...
        'FiniteDifferenceType', 'central', 'FiniteDifferenceStepSize', 1e-6, ...
        'Algorithm', 'sqp', 'HessianApproximation', 'lbfgs');


    obj = @(params) objectiveFunction(params, afStar, rfStar, r0, vfStar, v0, gConst);
    [optParams, optCost] = fmincon(obj, paramsX0, Aineq, bineq, [], [], lb, ub, [], fminconOptions);
    
end

function cost = objectiveFunction(params, afStar, rfStar, r, vfStar, v, gConst)
    gamma  = params(1);
    kr     = params(2);
    tgo   = params(3);

    gamma1 = gamma;
    gamma2 = kr/(gamma+2) - 2;

    [c1, c2] = calculateCoeffs(r, v, tgo, gamma1, gamma2, afStar, rfStar, vfStar, gConst);

    %% Simpson 1/3 Rule, goes fast (< 1microsecond per loop) but tgo is inaccurate due to only having 3 data points
    % h = (tgo - 0)/2;
    % aT_0 = afStar + c1*0^gamma1 + c2*0^gamma2;
    % aT_h = afStar + c1*h^gamma1 + c2*h^gamma2;
    % aT_tgo = afStar + c1*tgo^gamma1 + c2*tgo^gamma2;

    % aT or aT^2 for cost, set value to 1 or 2
    % tic
    % aTexp = 2;
    % 
    % if aTexp == 1
    %     aT_0 = norm(aT_0);
    %     aT_h = norm(aT_h);
    %     aT_tgo = norm(aT_tgo);
    % elseif aTexp == 2
    %     aT_0 = dot(aT_0,aT_0);
    %     aT_h = dot(aT_h,aT_h);
    %     aT_tgo = dot(aT_tgo,aT_tgo);
    % else
    %     error("Cost function not set to aT or aT^2");
    % end
    % 
    % simpson = h*(aT_0 + 4*aT_h + aT_tgo)/3;
    % toc
    %% Simpson Composite 1/3 Rule
    tspan = linspace(0,tgo,997);
    aT = afStar + c1*tspan.^gamma1 + c2*tspan.^gamma2;
    simpson = simpsonComp13Integral(tspan,dot(aT,aT));
    %%
    cost = simpson;

end