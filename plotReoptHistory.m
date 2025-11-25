function plotReoptHistory(optHistory, tTraj, stateTraj, aTList, t_plot, PDIState, planetaryParams, targetState, vehicleParams, optimizationParams, betaParam)
    % plotReoptHistory Reconstructs and plots the optimization state at a specific time in the reoptimization history.
    %
    % Inputs:
    %   optHistory: Table or matrix of optimization history from simReOpt.
    %               Columns/Fields: t_elapsed, gamma, gamma2, kr, tgo.
    %   tTraj: Vector of time points from simulation (ND or Dim? simReOpt returns ND usually).
    %          NOTE: This function assumes tTraj is in Nondimensional Time if it comes from simReOpt.
    %   stateTraj: Matrix of state trajectory from simReOpt.
    %   aTList: Matrix of thrust acceleration history (3xN) from simReOpt.
    %   t_plot: Time (in seconds, dimensional) to plot. The function will find the closest reoptimization point.
    %   PDIState: Struct defining the initial state.
    %   planetaryParams: Struct defining planetary parameters.
    %   targetState: Struct defining target state.
    %   vehicleParams: Struct defining vehicle parameters.
    %   optimizationParams: Struct defining optimization parameters.
    %   betaParam: Scalar optimization weight.

    addpath([pwd, '/CoordinateFunctions']);

    %% 1. Recreate Parameters (Copied from getParams.m)

    % PDI Breakout:
    altitude_km        = PDIState.altitude / 1000;
    lonInitDeg         = PDIState.lonInitDeg;
    latInitDeg         = PDIState.latInitDeg;
    inertialVelocity   = PDIState.inertialVelocity;
    flightPathAngleDeg = PDIState.flightPathAngleDeg;
    azimuth            = PDIState.azimuth * pi / 180;

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

    % Default Constants
    gEarth = 9.81;
    L_ref = 10000;
    T_ref = sqrt(L_ref/gPlanet);
    A_ref = gPlanet;
    V_ref = L_ref / T_ref;
    M_ref = 15103.0;

    %% Setup and NonDimensionalization

    % Convert PDI initial conditions into the MCMF frame
    [r0Dim, v0Dim] = PDI2MCMF(altitude_km, lonInitDeg, latInitDeg, ...
                                       landingLonDeg, landingLatDeg, ...
                                       inertialVelocity, flightPathAngleDeg, azimuth, rPlanet);

    rfDim = 10000*ENU2MCMF(rfLanding/10000,landingLatDeg,landingLonDeg,true);
    vfDim = ENU2MCMF(vfLanding,landingLatDeg,landingLonDeg,false);
    afDim = ENU2MCMF(afLanding,landingLatDeg,landingLonDeg,false);

    rPlanetND     = rPlanet / L_ref;
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

    % Problem Parameters Struct
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

    % Reference Values
    refVals = struct;
    refVals.L_ref = L_ref;
    refVals.T_ref = T_ref;
    refVals.A_ref = A_ref;
    refVals.V_ref = V_ref;
    refVals.M_ref = M_ref;

    % Non dimensional parameters (base structure)
    nonDimParams = struct;
    nonDimParams.rMoonND = rPlanetND;
    nonDimParams.rfStarND = rfStarND;
    nonDimParams.vfStarND = vfStarND;
    nonDimParams.afStarND = afStarND;
    nonDimParams.gConst = gConst;
    nonDimParams.mMinND = dryMassND;
    nonDimParams.ispND = ispND;
    nonDimParams.maxThrustND = maxThrustND;
    nonDimParams.minThrustND = minThrustND;

    %% 2. Locate Data Point

    % Handle optHistory format (Table or Matrix)
    if istable(optHistory)
        t_hist = optHistory.t_elapsedND * T_ref;
        gamma_hist = optHistory.gamma1;
        gamma2_hist = optHistory.gamma2;
        kr_hist = optHistory.kr;
        tgo_hist = optHistory.tgoND;
    else
        % Matrix [t_elapsed, gamma, gamma2, kr, tgo]
        t_hist = optHistory(:,1) * T_ref;
        gamma_hist = optHistory(:,2);
        gamma2_hist = optHistory(:,3);
        kr_hist = optHistory(:,4);
        tgo_hist = optHistory(:,5);
    end

    % Find closest index
    [min_diff, idx] = min(abs(t_hist - t_plot));

    if min_diff > 1.0 % Warning if more than 1 second away
        fprintf('Warning: Closest reoptimization point is %.2f s away from requested time %.2f s.\n', min_diff, t_plot);
    end

    t_selected = t_hist(idx);
    gamma_sel = gamma_hist(idx);
    gamma2_sel = gamma2_hist(idx);
    kr_sel = kr_hist(idx);
    tgo_sel = tgo_hist(idx); % This is ND

    fprintf('Plotting for Reoptimization at t = %.2f s (Index %d)\n', t_selected, idx);
    fprintf('Selected Params: gamma=%.4f, gamma2=%.4f, tgo=%.2f s\n', gamma_sel, gamma2_sel, tgo_sel * T_ref);

    %% 3. Extract State at that time
    % Interpolate to find state at t_selected (sim time)
    % tTraj is assumed to be ND if from simReOpt

    t_selected_ND = t_selected / T_ref;

    % Ensure unique points for interpolation
    [tTrajUnique, uniqueIdx] = unique(tTraj);
    stateTrajUnique = stateTraj(uniqueIdx, :);

    if isempty(tTrajUnique)
         error('tTraj is empty. Cannot interpolate state.');
    end

    if t_selected_ND <= tTrajUnique(1)
        %warning('Requested time is before simulation start. Using first point.');
        X0 = stateTrajUnique(1, :)';
    elseif t_selected_ND >= tTrajUnique(end)
        %warning('Requested time is after simulation end. Using last point.');
        X0 = stateTrajUnique(end, :)';
    else
        X0 = interp1(tTrajUnique, stateTrajUnique, t_selected_ND, 'linear')';
    end

    r0_curr = X0(1:3);
    v0_curr = X0(4:6);
    m0_curr = X0(7);

    nonDimParams.r0ND = r0_curr;
    nonDimParams.v0ND = v0_curr;
    nonDimParams.m0ND = m0_curr;

    %% 4. Reconstruct Optimal Trajectory (Prediction)

    nodeCount = optimizationParams.nodeCount;

    [c1, c2] = calculateCoeffs(r0_curr, v0_curr, tgo_sel, gamma_sel, gamma2_sel, afStarND, rfStarND, vfStarND, gConst);

    tgospan = linspace(0, tgo_sel, nodeCount);

    % Generate Trajectory Arrays
    aTOptim = afStarND + c1*tgospan.^gamma_sel + c2*tgospan.^gamma2_sel;

    phi1hat = (tgospan.^(gamma_sel+2))./((gamma_sel+1)*(gamma_sel+2));
    phi2hat = (tgospan.^(gamma2_sel+2))./((gamma2_sel+1)*(gamma2_sel+2));
    phi1bar = (tgospan.^(gamma_sel+1))./(gamma_sel+1);
    phi2bar = (tgospan.^(gamma2_sel+1))./(gamma2_sel+1);

    rdOptim = rfStarND + c1*phi1hat + c2*phi2hat - vfStarND.*tgospan + 0.5*(gConst+afStarND).*tgospan.^2;
    vdOptim = vfStarND + c1*phi1bar + c2*phi2bar -(gConst+afStarND).*tgospan;

    aTNorm = vecnorm(aTOptim, 2, 1);
    Q = cumtrapz(tgospan, aTNorm./ispND);
    Q = Q(end) - Q;
    mOptim = m0_curr .* exp(-Q);

    %% 5. Prepare Data for Plotting

    optParams = [gamma_sel, gamma2_sel, tgo_sel];
    optCost = NaN;

    flag_thrustGotLimited = false;
    unconstrained = [];

    % Call plotting
    % Note: tTraj passed to plotting should match stateTraj.
    plotting(tTraj, stateTraj, optParams, optCost, aTOptim, mOptim, rdOptim, vdOptim, aTList, refVals, problemParams, nonDimParams, optimizationParams, flag_thrustGotLimited, unconstrained);

end