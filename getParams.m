function [gammaOpt, gamma2Opt, krOpt, tgoOpt, aTOptim, exitflag, optFuelCost, simFuelCost, aTSim, finalPosSim, optHistory, ICstates, exitFlags, problemParams, nonDimParams, refVals] = ...
    getParams(PDIState, planetaryParams, targetState, vehicleParams, optimizationParams, betaParam, doPlots, verboseOutput, dispersion, runSimulation)

    addpath([pwd, '/CoordinateFunctions']);

    %% 1. Initialization & Defaults

    % Planetary Constants
    rPlanet = planetaryParams.rPlanet;
    gPlanet = planetaryParams.gPlanet;

    % Reference & Scaling Constants
    gEarth = 9.81;
    L_ref = 10000;
    T_ref = sqrt(L_ref/gPlanet);
    A_ref = gPlanet;
    V_ref = L_ref / T_ref;
    M_ref = 15103.0;

    % Extract State & Vehicle Parameters
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

    %% 2. Coordinate Transformation & Non-Dimensionalization
    [r0Dim, v0Dim] = PDI2MCMF(altitude_km, lonInitDeg, latInitDeg, ...
                                       landingLonDeg, landingLatDeg, ...
                                       inertialVelocity, flightPathAngleDeg, azimuth, rPlanet);

    rfDim = 10000 * ENU2MCMF(rfLanding/10000, landingLatDeg, landingLonDeg, true);
    vfDim = ENU2MCMF(vfLanding, landingLatDeg, landingLonDeg, false);
    afDim = ENU2MCMF(afLanding, landingLatDeg, landingLonDeg, false);

    rPlanetND   = rPlanet / L_ref;
    r0ND        = r0Dim / L_ref;
    v0ND        = v0Dim / V_ref;
    rfStarND    = rfDim / L_ref;
    vfStarND    = vfDim / V_ref;
    afStarND    = afDim / A_ref;
    gConst      = -(rPlanetND^2) * rfStarND / (norm(rfStarND)^3); 
    massInitND  = massInit / M_ref;
    dryMassND   = dryMass / M_ref;
    ispND       = isp * gEarth / (V_ref);
    maxThrustND = maxThrust / (M_ref * A_ref);
    minThrustND = minThrust / (M_ref * A_ref);
    delta_tND   = delta_t / T_ref;

    %% 3. Structs Packing
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

    %% 4. Optimization
    
    paramsX0 = [0.3, 0.4, 6];
    reopt = optimizationParams.updateOpt;

    [optParams, optCost, aTOptim, mOptim, rdOptim, vdOptim, exitflag] = ...
        optimizationLoop(paramsX0, betaParam, problemParams, nonDimParams, optimizationParams, refVals, delta_tND, verboseOutput, dispersion);
    if exitflag ~= 1
        fprintf("\n First Optimization Converged to flag ~=2, rerunning optimization starting from first rounds parameters:\n");
        [optParams, optCost, aTOptim, mOptim, rdOptim, vdOptim, exitflag] = ...
        optimizationLoop(optParams, betaParam, problemParams, nonDimParams, optimizationParams, refVals, delta_tND, verboseOutput, dispersion);
    end

    gammaOpt = optParams(1);
    gamma2Opt = optParams(2);
    krOpt = (gamma2Opt+2)*(gammaOpt+2);
    tgoOpt = optParams(3) * T_ref;
    optFuelCost = (mOptim(end)-mOptim(1))*M_ref;

    %% 5. Comparison/Secondary Data Generation (Static or Unconstrained)
    % Generates a secondary dataset for plotting comparisons.
    needsSecondary = doPlots && ~dispersion;
    secondaryData = struct();

    if needsSecondary
        if reopt
            % Mode 1: Re-Opt is Enabled. 
            % Secondary Data = Static Baseline (Same params, no updates).
            if verboseOutput; fprintf("\nGenerating Static Baseline for Comparison...\n"); end
            
            optParamsSec = optParams;
            optCostSec   = optCost;
            aTOptimSec   = aTOptim;
            mOptimSec    = mOptim;
            rdOptimSec   = rdOptim;
            vdOptimSec   = vdOptim;
            
            [tTrajSec, stateTrajSec, aTListSec, flag_thrustGotLimitedSec] = ...
                closedLoopSim(gammaOpt, gamma2Opt, tgoOpt/T_ref, problemParams, nonDimParams, refVals, delta_tND);
                
        else
            % Mode 2: Re-Opt is Disabled (Single Run).
            % Secondary Data = Unconstrained Optimization (No glideslope/pointing).
            if verboseOutput; fprintf("Generating Unconstrained Trajectory for Comparison...\n"); end
            
            optimizationParamsUC = optimizationParams;
            optimizationParamsUC.glideSlopeEnabled = false;
            optimizationParamsUC.pointingEnabled = false;
            
            [optParamsSec, optCostSec, aTOptimSec, mOptimSec, rdOptimSec, vdOptimSec] = ...
                optimizationLoop(paramsX0, betaParam, problemParams, nonDimParams, optimizationParamsUC, refVals, delta_tND, false, false);
                
            gammaOptSec = optParamsSec(1);
            gamma2OptSec = optParamsSec(2);
            tgoOptSec = optParamsSec(3) * T_ref;
            
            [tTrajSec, stateTrajSec, aTListSec, flag_thrustGotLimitedSec] = ...
                closedLoopSim(gammaOptSec, gamma2OptSec, tgoOptSec/T_ref, problemParams, nonDimParams, refVals, delta_tND);
        end
        
        % Pack Secondary Data Structure
        secondaryData = struct('tTraj', tTrajSec, 'stateTraj', stateTrajSec, 'optParams', optParamsSec, ...
                              'optCost', optCostSec, 'aTOptim', aTOptimSec, 'mOptim', mOptimSec, ...
                              'rdOptim', rdOptimSec, 'vdOptim', vdOptimSec, 'aTList', aTListSec, ...
                              'flag_thrustGotLimited', flag_thrustGotLimitedSec);
    end

    %% 6. Simulation
    tTraj = [];
    stateTraj = [];
    aTSim = [];
    flag_thrustGotLimited = false;
    simFuelCost = [];
    finalPosSim = [];
    optHistory = [];
    ICstates = [];
    exitFlags = [];

    if runSimulation
        if ~reopt
            % Static Simulation
            [tTraj, stateTraj, aTSim, flag_thrustGotLimited] = ...
                closedLoopSim(gammaOpt, gamma2Opt, tgoOpt/T_ref, problemParams, nonDimParams, refVals, delta_tND);
        else
            % Re-Optimization Simulation
            [tTraj, stateTraj, aTSim, flag_thrustGotLimited, optHistory, ICstates, exitFlags] = ...
                simReOpt(gammaOpt, gamma2Opt, tgoOpt/T_ref, problemParams, nonDimParams, refVals, delta_tND, optimizationParams, betaParam, verboseOutput);
            
            % Format History Output
            if ~isempty(optHistory)
                optHistory = array2table(optHistory);
                optHistory.Properties.VariableNames(1:5) = {'t_elapsedND','gamma1','gamma2','kr','tgoND'};
                ICstates = array2table(ICstates');
                ICstates.Properties.VariableNames(1:7) = {'r0_X', 'r0_Y', 'r0_Z', 'v0_X', 'v0_Y', 'v0_Z', 'm0'};
            end
        end
        
        simFuelCost = M_ref * (stateTraj(1,7) - stateTraj(end,7));
        finalPosSim = MCMF2ENU(stateTraj(end,1:3)' * L_ref, landingLatDeg, landingLonDeg, true, true);
    end

    %% 7. Plotting
    if doPlots
        plotting(tTraj, stateTraj, optParams, optCost, aTOptim, mOptim, rdOptim, vdOptim, aTSim, ...
            refVals, problemParams, nonDimParams, optimizationParams, flag_thrustGotLimited, ...
            secondaryData, optHistory, ICstates);
    end

end