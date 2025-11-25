function [gammaOpt, gamma2Opt, krOpt, tgoOpt, aTOptim, exitflag, optFuelCost, simFuelCost, aTSim, finalPosSim, optHistory, exitFlags] = ...
    getParams(PDIState, planetaryParams, targetState, vehicleParams, optimizationParams, simulationParams, constants)
%getParams Main Function to Run for FP2PDG Optimization Single or Stat
%   This function takes in the initial state of the vehicle, planetary
%   parameters, target state, vehicle parameters, optimization parameters,
%   and simulation parameters and returns the optimal guidance parameters
%   and simulation results.
%
%   Inputs:
%       PDIState - Struct containing the initial state of the vehicle
%       planetaryParams - Struct containing planetary parameters
%       targetState - Struct containing the target state
%       vehicleParams - Struct containing vehicle parameters
%       optimizationParams - Struct containing optimization parameters
%       simulationParams - Struct containing simulation parameters
%       constants - Struct containing physical and simulation constants
%
%   Outputs:
%       gammaOpt, gamma2Opt, krOpt, tgoOpt - Optimal guidance parameters
%       aTOptim - Optimal thrust acceleration profile
%       exitflag - Exit flag from the optimization
%       optFuelCost - Fuel cost from the optimization
%       simFuelCost - Fuel cost from the simulation
%       aTSim - Thrust acceleration profile from the simulation
%       finalPosSim - Final position from the simulation
%       optHistory - History of the optimization parameters
%       exitFlags - Exit flags from the re-optimization steps

    addpath([pwd, '/CoordinateFunctions']);

    betaParam = simulationParams.beta;
    doPlots = simulationParams.doPlotting;
    verboseOutput = simulationParams.verboseOutput;
    runSimulation = simulationParams.runSimulation;

    altitude_km        = PDIState.altitude / 1000;
    lonInitDeg         = PDIState.lonInitDeg;
    latInitDeg         = PDIState.latInitDeg;
    inertialVelocity   = PDIState.inertialVelocity;
    flightPathAngleDeg = PDIState.flightPathAngleDeg;
    azimuth            = PDIState.azimuth * pi / 180;

    landingLonDeg      = targetState.landingLonDeg;
    landingLatDeg      = targetState.landingLatDeg;
    rfLanding          = targetState.rfLanding;
    vfLanding          = targetState.vfLanding;
    afLanding          = targetState.afLanding;
    delta_t            = targetState.delta_t;

    massInit           = vehicleParams.massInit;
    dryMass            = vehicleParams.dryMass;
    isp                = vehicleParams.isp;
    maxThrust          = vehicleParams.maxThrust;
    minThrust          = vehicleParams.minThrust;

    rPlanet            = planetaryParams.rPlanet;
    gPlanet            = planetaryParams.gPlanet;

    gEarth             = constants.gEarth;
    L_ref              = constants.L_ref;
    T_ref              = constants.T_ref;
    A_ref              = constants.A_ref;
    V_ref              = constants.V_ref;
    M_ref              = constants.M_ref;

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
    gConst      = -(rPlanetND^2) * rfStarND / (norm(rfStarND)^3);
    massInitND        = massInit / M_ref;
    dryMassND      = dryMass / M_ref;
    ispND       = isp * gEarth / (V_ref);
    maxThrustND = maxThrust/(M_ref*A_ref);
    minThrustND = minThrust/(M_ref*A_ref);
    delta_tND = delta_t / T_ref;

    problemParams = struct;
    problemParams.r0Dim = r0Dim;
    problemParams.v0Dim = v0Dim;
    problemParams.rMoon = rPlanet;
    problemParams.gMoon = gPlanet;
    problemParams.g0 = gEarth;
    problemParams.rfDim = rfDim;
    problemParams.vfDim = vfDim;
    problemParams.afDim = afDim;
    problemParams.massInitDim = massInit;
    problemParams.dryMassDim = dryMass;
    problemParams.ispDim = isp;
    problemParams.maxThrustDim = maxThrust;
    problemParams.minThrustDim = minThrust;
    problemParams.landingLatDeg = landingLatDeg;
    problemParams.landingLonDeg = landingLonDeg;

    refVals = struct;
    refVals.L_ref = L_ref;
    refVals.T_ref = T_ref;
    refVals.A_ref = A_ref;
    refVals.V_ref = V_ref;
    refVals.M_ref = M_ref;

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

    paramsX0 = optimizationParams.initialGuess;
    reopt = optimizationParams.updateOpt;

    [optParams, optCost, aTOptim, mOptim, rdOptim, vdOptim, exitflag] = ...
        optimizationLoop(paramsX0, betaParam, problemParams, nonDimParams, optimizationParams, refVals, delta_tND, verboseOutput, false);

    gammaOpt = optParams(1);
    gamma2Opt = optParams(2);
    krOpt = (gamma2Opt+2)*(gammaOpt+2);
    tgoOpt = optParams(3) * T_ref;
    optFuelCost = (mOptim(end)-mOptim(1))*M_ref;

    needsUnConOpt = doPlots;
    needsSim = runSimulation;
    needsPlotting = doPlots;

    unconstrained = struct();
    if needsUnConOpt
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
        [tTrajUC, stateTrajUC, aTListUC, flag_thrustGotLimitedUC] = closedLoopSim(gammaOptUC, gamma2OptUC, tgoOptUC/T_ref, problemParams, nonDimParams, refVals, delta_tND, simulationParams);
        ucSimFuelCost = M_ref*(stateTrajUC(1,7) - stateTrajUC(end,7));
        unconstrained = struct('tTraj', tTrajUC, 'stateTraj', stateTrajUC, 'optParams', optParamsUC, ...
                              'optCost', optCostUC, 'aTOptim', aTOptimUC, 'mOptim', mOptimUC, ...
                              'rdOptim', rdOptimUC, 'vdOptim', vdOptimUC, 'aTList', aTListUC, ...
                              'flag_thrustGotLimited', flag_thrustGotLimitedUC);
        fprintf("-----------------------------------------------\n");
    end

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
            [tTraj, stateTraj, aTSim, flag_thrustGotLimited] = closedLoopSim(gammaOpt, gamma2Opt, tgoOpt/T_ref, problemParams, nonDimParams, refVals, delta_tND, simulationParams);
        else
            [tTraj, stateTraj, aTSim, flag_thrustGotLimited, optHistory, exitFlags] = simReOpt(gammaOpt, gamma2Opt, tgoOpt/T_ref, problemParams, nonDimParams, refVals, delta_tND, optimizationParams, betaParam, verboseOutput, simulationParams);
            if ~isempty(optHistory)
                optHistory = array2table(optHistory);
                optHistory.Properties.VariableNames(1:5) = {'t_elapsedND','gamma1','gamma2','kr','tgoND'};
            end
        end
        simFuelCost = M_ref * (stateTraj(1,7) - stateTraj(end,7));
        finalPosSim = MCMF2ENU(stateTraj(end,1:3)' * L_ref, landingLatDeg, landingLonDeg, true, true);
    end

    if needsPlotting
        plotting(tTraj, stateTraj, optParams, optCost, aTOptim, mOptim, rdOptim, vdOptim, aTSim, refVals, problemParams, nonDimParams, optimizationParams, flag_thrustGotLimited, unconstrained);
    end
end
