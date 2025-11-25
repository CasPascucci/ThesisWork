function [gammaOpt, gamma2Opt, krOpt, tgoOpt, aTOptim, exitflag, optFuelCost, simFuelCost, aTSim, finalPosSim, optHistory, exitFlags] = ...
    getParams(PDIState, planetaryParams, targetState, vehicleParams, optimizationParams, betaParam, doPlots, verboseOutput, dispersion, runSimulation)

    addpath([pwd, '/CoordinateFunctions']);

    % Default Constants (not user defined)
    gEarth = 9.81;
    L_ref = 10000;
    T_ref = sqrt(L_ref/planetaryParams.gPlanet);
    A_ref = planetaryParams.gPlanet;
    V_ref = L_ref / T_ref;
    M_ref = 15103.0;

    if nargin < 9
        dispersion = false;
        runSimulation = true;
    end
    if nargin < 10
        runSimulation = true;
    end

    % PDI Breakout:
    altitude_km        = PDIState.altitude / 1000;
    lonInitDeg         = PDIState.lonInitDeg;
    latInitDeg         = PDIState.latInitDeg;
    inertialVelocity   = PDIState.inertialVelocity;
    flightPathAngleDeg = PDIState.flightPathAngleDeg;
    azimuth            = PDIState.azimuth * pi / 180; % conver to rad for processing

    % Target State Breakout:
    landingLonDeg      = targetState.landingLonDeg;
    landingLatDeg      = targetState.landingLatDeg;
    rfLanding          = targetState.rfLanding;
    vfLanding          = targetState.vfLanding;
    afLanding          = targetState.afLanding;
    delta_t            = targetState.delta_t;

    % Vehicle Params Breakout:
    massInit           = vehicleParams.massInit;
    dryMass            = vehicleParams.dryMass;
    isp                = vehicleParams.isp;
    maxThrust          = vehicleParams.maxThrust;
    minThrust          = vehicleParams.minThrust;

    % Planetary Params Breakout:
    rPlanet            = planetaryParams.rPlanet;
    gPlanet            = planetaryParams.gPlanet;

    %% Setup and NonDimensionalization

    %Convert PDI initial conditions into the MCMF frame of the problem
    [r0Dim, v0Dim] = PDI2MCMF(altitude_km, lonInitDeg, latInitDeg, ...
                                       landingLonDeg, landingLatDeg, ...
                                       inertialVelocity, flightPathAngleDeg, azimuth, rPlanet);

    rfDim = 10000*ENU2MCMF(rfLanding/10000,landingLatDeg,landingLonDeg,true);
    vfDim = ENU2MCMF(vfLanding,landingLatDeg,landingLonDeg,false);
    afDim = ENU2MCMF(afLanding,landingLatDeg,landingLonDeg,false);

    rPlanetND     = rPlanet / L_ref;
    r0ND        = r0Dim / L_ref;
    v0ND        = v0Dim / V_ref;
    rfStarND    = rfDim / L_ref;
    vfStarND    = vfDim / V_ref;
    afStarND    = afDim / A_ref;
    gConst      = -(rPlanetND^2) * rfStarND / (norm(rfStarND)^3); % Constant landing site gravity for optimization loop
    massInitND        = massInit / M_ref;
    dryMassND      = dryMass / M_ref;
    ispND       = isp * gEarth / (V_ref);
    maxThrustND = maxThrust/(M_ref*A_ref);
    minThrustND = minThrust/(M_ref*A_ref);
    delta_tND = delta_t / T_ref;

    %% Create Structs to be Passed Along

    % Problem Parameters Struct
    problemParams = struct;
    problemParams.r0Dim = r0Dim; % m
    problemParams.v0Dim = v0Dim; % m
    problemParams.rMoon = rPlanet; % m
    problemParams.gMoon = gPlanet; % m/s^2
    problemParams.g0 = gEarth; % m/s^2
    problemParams.rfDim = rfDim; % m
    problemParams.vfDim = vfDim; % m/s
    problemParams.afDim = afDim; % m/s^2
    problemParams.massInitDim = massInit; % kg
    problemParams.dryMassDim = dryMass; % kg
    problemParams.ispDim = isp; % s
    problemParams.maxThrustDim = maxThrust; % N
    problemParams.minThrustDim = minThrust; % N
    problemParams.landingLatDeg = landingLatDeg;
    problemParams.landingLonDeg = landingLonDeg;

    % Reference Values for Non-Dim
    refVals = struct;
    refVals.L_ref = L_ref;
    refVals.T_ref = T_ref;
    refVals.A_ref = A_ref;
    refVals.V_ref = V_ref;
    refVals.M_ref = M_ref;

    % Non dimensional parameters
    nonDimParams = struct;
    nonDimParams.rMoonND = rPlanetND;
    nonDimParams.r0ND = r0ND;
    nonDimParams.v0ND = v0ND;
    nonDimParams.rfStarND = rfStarND;
    nonDimParams.vfStarND = vfStarND;
    nonDimParams.afStarND = afStarND;
    nonDimParams.gConst = gConst;
    nonDimParams.m0ND = massInitND;
    nonDimParams.mMinND = dryMassND;
    nonDimParams.ispND = ispND;
    nonDimParams.maxThrustND = maxThrustND;
    nonDimParams.minThrustND = minThrustND;

    %% Optimization
    paramsX0 = [0.3, 0.4, 6];
    reopt = optimizationParams.updateOpt;

    [optParams, optCost, aTOptim, mOptim, rdOptim, vdOptim, exitflag] = ...
        optimizationLoop(paramsX0, betaParam, problemParams, nonDimParams, optimizationParams, refVals, delta_tND, verboseOutput, dispersion);

    gammaOpt = optParams(1);
    gamma2Opt = optParams(2);
    krOpt = (gamma2Opt+2)*(gammaOpt+2);
    tgoOpt = optParams(3) * T_ref;
    optFuelCost = (mOptim(end)-mOptim(1))*M_ref;

    %Determine what needs to be run
    needsUnConOpt = doPlots && ~dispersion;
    needsSim = runSimulation;
    needsPlotting = doPlots;


    unconstrained = struct();
    if needsUnConOpt % If plotting, sho unconstrained for comparison
        fprintf("-----------------------------------------------\n");
        fprintf("Unconstrained OPT for Plotting\n");
        optimizationParamsUC = optimizationParams;
        optimizationParamsUC.glideSlopeEnabled = false;
        optimizationParamsUC.pointingEnabled = false;
        [optParamsUC, optCostUC, aTOptimUC, mOptimUC, rdOptimUC, vdOptimUC] = optimizationLoop(paramsX0, betaParam, problemParams, nonDimParams, optimizationParamsUC, refVals, delta_tND, false, false);
        ucOptFuelCost = (mOptimUC(end)-mOptimUC(1))*M_ref;
        gammaOptUC = optParamsUC(1);
        gamma2OptUC = optParamsUC(2);
        tgoOptUC = optParamsUC(3) * T_ref;
        [tTrajUC, stateTrajUC, aTListUC, flag_thrustGotLimitedUC] = closedLoopSim(gammaOptUC, gamma2OptUC, tgoOptUC/T_ref, problemParams, nonDimParams, refVals, delta_tND);
        ucSimFuelCost = M_ref*(stateTrajUC(1,7) - stateTrajUC(end,7));
        unconstrained = struct('tTraj', tTrajUC, 'stateTraj', stateTrajUC, 'optParams', optParamsUC, ...
                              'optCost', optCostUC, 'aTOptim', aTOptimUC, 'mOptim', mOptimUC, ...
                              'rdOptim', rdOptimUC, 'vdOptim', vdOptimUC, 'aTList', aTListUC, ...
                              'flag_thrustGotLimited', flag_thrustGotLimitedUC);
        fprintf("-----------------------------------------------\n");
    end

    % If doing simulation:
    tTraj = [];
    stateTraj = [];
    aTSim = [];
    flag_thrustGotLimited = false;
    simFuelCost = [];
    finalPosSim = [];
    optHistory = [];
    exitFlags = [];

    if needsSim
        if ~reopt
            [tTraj, stateTraj, aTSim, flag_thrustGotLimited] = closedLoopSim(gammaOpt, gamma2Opt, tgoOpt/T_ref, problemParams, nonDimParams, refVals, delta_tND);
        else
            [tTraj, stateTraj, aTSim, flag_thrustGotLimited, optHistory, exitFlags] = simReOpt(gammaOpt, gamma2Opt, tgoOpt/T_ref, problemParams, nonDimParams, refVals, delta_tND, optimizationParams, betaParam, verboseOutput);
            % Format optHistory as table
            if ~isempty(optHistory)
                optHistory = array2table(optHistory);
                optHistory.Properties.VariableNames(1:5) = {'t_elapsedND','gamma1','gamma2','kr','tgoND'};
            end
        end
        simFuelCost = M_ref * (stateTraj(1,7) - stateTraj(end,7));
        finalPosSim = MCMF2ENU(stateTraj(end,1:3)' * L_ref, landingLatDeg, landingLonDeg, true, true);
    end

    %If doing Plots
    if needsPlotting
        plotting(tTraj, stateTraj, optParams, optCost, aTOptim, mOptim, rdOptim, vdOptim, aTSim, refVals, problemParams, nonDimParams, optimizationParams, flag_thrustGotLimited, unconstrained);
    end
end