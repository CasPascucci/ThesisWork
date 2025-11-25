clear all; close all; clc; format short
%% Function to use for Single Runs
% All values are dimensional
PDIState = struct;
PDIState.altitude           = 15240; % m
PDIState.lonInitDeg         = 41.85; % deg
PDIState.latInitDeg         = -71.59; % deg
PDIState.inertialVelocity   = 1693.8; % m/s
PDIState.flightPathAngleDeg = 0; % deg
PDIState.azimuth            = 180; % deg

planetaryParams = struct;
planetaryParams.rPlanet = 1736.01428 * 1000; % m
planetaryParams.gPlanet = 1.622; % m/s^2
planetaryParams.gEarth = 9.81; % m/s^2

vehicleParams = struct;
vehicleParams.massInit = 15103.0; % kg
vehicleParams.dryMass = vehicleParams.massInit - 8248; % kg
vehicleParams.isp = 311; % s
vehicleParams.maxThrust = 45000; % N
vehicleParams.minThrust = 4500; % N

targetState = struct;
targetState.landingLonDeg      = 41.85; % deg
targetState.landingLatDeg      = -90.0; % deg
targetState.rfLanding = [0;0;0]; % m, TOPO frame centered at landing site lat/lon, and lunar radius
targetState.vfLanding = [0;0;-1]; % m/s
targetState.afLanding = [0;0;2*planetaryParams.gPlanet]; % m/s^2
%targetState.afLanding = [0;0;vehicleParams.maxThrust/(6710)]; % This
                                                               %afStar seems to be too high for a solvable path
targetState.delta_t   = 5; % seconds dim, for btt, not implemented

optimizationParams = struct;
optimizationParams.nodeCount = 301; %Count must be odd for Simpson
optimizationParams.glideSlopeFinalTheta = 45; %deg
optimizationParams.glideSlopeHigh = 500; %m
optimizationParams.glideSlopeLow = 250; %m
optimizationParams.glideSlopeCutoff = 50; %m
optimizationParams.glideSlopeEnabled = false;
optimizationParams.pointingEnabled = true;
optimizationParams.maxTiltAccel = 2; % deg/s^2
optimizationParams.minPointing = 10; %deg, floor for pointing constraint
optimizationParams.updateFreq = 10;
optimizationParams.updateStop = 60;
optimizationParams.updateOpt = false; % Only Applies if Sim is also set to turn on
optimizationParams.initialGuess = [0.3, 0.4, 6];
optimizationParams.Aineq = [-1  0  0; 1 -1  0; 0  0 -1];
optimizationParams.bineq = [0; -1e-4; -0.01];
optimizationParams.lb = [0, 0, 0.01];
optimizationParams.ub = [10, 10, 15];
optimizationParams.marginTop = 0.95;

simulationParams = struct;
simulationParams.beta = 0.5;
simulationParams.runSimulation = true;
simulationParams.doPlotting = true;
simulationParams.verboseOutput = true;
simulationParams.minNodeCount = 21;
simulationParams.gammaChangeThreshold = 0.5;
simulationParams.minTime = 0.2;
simulationParams.odeRelTol = 1e-6;
simulationParams.odeAbsTol = 1e-6;

constants = struct;
constants.gEarth = 9.81;
constants.L_ref = 10000;
constants.T_ref = sqrt(constants.L_ref/planetaryParams.gPlanet);
constants.A_ref = planetaryParams.gPlanet;
constants.V_ref = constants.L_ref / constants.T_ref;
constants.M_ref = 15103.0;

tic
[gammaOpt, gamma2Opt, krOpt, tgoOptSec,~,~, optFuelCost, simFuelCost, aTSim,finalPosSim, optHistory, exitFlags] = getParams(PDIState,...
    planetaryParams, targetState, vehicleParams, optimizationParams, simulationParams, constants);
toc
% tgoOpt returned in seconds
gammaOpt
gamma2Opt
krOpt
tgoOptSec
optFuelCost
simFuelCost