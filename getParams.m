function [gammaOpt, krOpt, tgoOpt, aTOptim, exitflag, optFuelCost, simFuelCost, aTList] = getParams(PDIState, planetaryParams, targetState, vehicleParams, optimizationParams, betaParam, doPlots, verboseOutput, dispersion)
    addpath([pwd, '/CoordinateFunctions']);
%% Main Function to Run for FP2PDG Optimization Single or Stat
% Inputs: All given in dimensional values unless stated otherwise
%   PDIState:
%       .altitude_km: Altitude of PDI, km
%       .lonInitDeg: Longitude of PDI, deg
%       .latInitDeg: Latitude of PDI, deg
%       .inertialVelocity: Inertial frame velocity at PDI, m/s
%       .flightPathAngleDeg: The angle between the flight path, and the
%                            horizontal plane, deg
%       .azimuth: Heading angle, measured clockwise from North, rad
%
%   planetaryParams:
%       .rPlanet: radius of the planet, m
%       .gPlanet: average surface gravity of the planet, m/s
%       .gEarth: value used for Earth gravity, used in
%                nondimensionalization of isp
%
%   targetState:
%       .landingLonDeg: Longitude of Landing
%       .landingLatDeg: Latitude of Landing
%       .rfLanding: target position vector at landing, topocentric frame
%       .vfLanding: target velocity vector at landing, topocentric frame
%       .afLanding: target acceleration vector at landing, topocentric
%                   frame
%
%   vehicleParams:
%       .massInit: Mass of lander at PDI, kg
%       .dryMass: Mass of lander without any fuel, kg
%       .isp: isp of the engine(s) aboard, given in seconds
%       .maxThrust: maximum total thrust of the lander, Newtons
%       .minThrust: minimum total thrust of the lander, Newtons
%
%   targetState:
%       .landingLonDeg: Longitude of Landing Site, degrees
%       .landingLatDeg: Latitude of Landing Site, degrees
%       .rfLanding: Topocentric coordinates of landing site, frame centered
%                   at the landing lat & lon at planetary radius, meters
%       .vfLanding: Topocentric frame desired final velocity, m/s
%       .afLanding: Topocentric frame desired final acceleration, m/s^2
%       .delta_t: Delta T value used for Beyond Termination Targeting,
%                 DISABLED, seconds
%
%   optimParams:
%       .nodeCount: number of nodes to use in optimization
%
%   beta: weight value between 0 and 1, weighing the objective function
%         towards smoothing the acceleration curve, or minimzing fuel cost
%         respecitively. Default value of 0.80
%
% Outputs:
%   gammaOpt, krOpt, tgoOpt: optimal parameters for FP2PDG from the given
%                            conditions
%   fuelCost: Estimated fuel cost from simulation of flight
if nargin < 6
    betaParam = 0.65;
elseif nargin < 7
    doPlots = false;
end
if nargin < 9
    dispersion = false;
end
    


% PDI Breakout:
altitude_km        = PDIState.altitude_km;
lonInitDeg         = PDIState.lonInitDeg;
latInitDeg         = PDIState.latInitDeg;
inertialVelocity   = PDIState.inertialVelocity;
flightPathAngleDeg = PDIState.flightPathAngleDeg;
azimuth            = PDIState.azimuth;

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

% Default Constants (not user defined)
gEarth = 9.81;
L_ref = 10000;
T_ref = sqrt(L_ref/gPlanet);
A_ref = gPlanet;
V_ref = L_ref / T_ref;
M_ref = 15103.0;

%% Setup and NonDimensionalization

%Convert PDI initial conditions into the MCMF frame of the problem
[r0Dim, v0Dim] = PDI2MCMF(altitude_km, lonInitDeg, latInitDeg, ...
                                   landingLonDeg, landingLatDeg, ...
                                   inertialVelocity, flightPathAngleDeg, azimuth, rPlanet);

% r0Dim = [411608.123492225;
%          368665.998391003;
%          -1659817.39432084];
% 
% v0Dim = [-1195.85903735503;
%          -1073.07538948840;
%          -536.193540835963];

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
paramsX0 = [1, 6.5, 6];

[optParams, optCost, aTOptim, mOptim, rdOptim, vdOptim, exitflag] = optimizationLoop(paramsX0, betaParam, problemParams, nonDimParams, optimizationParams, refVals, delta_tND, verboseOutput, dispersion);
gammaOpt = optParams(1);
krOpt = optParams(2);
tgoOpt = optParams(3) * T_ref;
optFuelCost = (mOptim(end)-mOptim(1))*M_ref;

% Run Unconstrained Opt and Sim for Plotting
if ~dispersion
    fprintf("-----------------------------------------------\n");
    fprintf("Unconstrained OPT for Plotting");
    optimizationParamsUC = optimizationParams;
    optimizationParamsUC.glideSlopeEnabled = false;
    optimizationParamsUC.pointingEnabled = false;
    [optParamsUC,optCostUC, aTOptimUC, mOptimUC, rdOptimUC, vdOptimUC] = optimizationLoop(paramsX0, betaParam, problemParams, nonDimParams, optimizationParamsUC, refVals, delta_tND, false);
    ucOptFuelCost = (mOptimUC(end)-mOptimUC(1))*M_ref;
    gammaOptUC = optParamsUC(1);
    krOptUC = optParamsUC(2);
    tgoOptUC = optParamsUC(3) * T_ref;
    [tTrajUC, stateTrajUC, aTListUC, flag_thrustGotLimitedUC] = closedLoopSim(gammaOptUC, krOptUC, tgoOptUC/T_ref, problemParams, nonDimParams, refVals, delta_tND);
    ucSimFuelCost = M_ref*(stateTrajUC(1,7) - stateTrajUC(end,7));
    unconstrained = struct('tTraj',tTrajUC,'stateTraj',stateTrajUC,'optParams',optParamsUC,'optCost',optCostUC,'aTOptim',aTOptimUC,'mOptim',mOptimUC,'rdOptim',rdOptimUC,'vdOptim',vdOptimUC,'aTList',aTListUC,'flag_thrustGotLimited',flag_thrustGotLimitedUC);
    fprintf("-----------------------------------------------");
end
% Plotting Handling
    if nargout > 5 || doPlots
        [tTraj, stateTraj, aTList, flag_thrustGotLimited] = closedLoopSim(gammaOpt, krOpt, tgoOpt/T_ref, problemParams, nonDimParams, refVals, delta_tND);
        simFuelCost = M_ref*(stateTraj(1,7) - stateTraj(end,7));
        if doPlots
            plotting(tTraj, stateTraj, optParams, optCost, aTOptim, mOptim, rdOptim, vdOptim, aTList, refVals, problemParams, nonDimParams, optimizationParams, flag_thrustGotLimited, unconstrained);
        end
    end
end