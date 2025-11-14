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
optimizationParams.glideSlopeEnabled = false;
optimizationParams.pointingEnabled = false;
optimizationParams.maxTiltAccel = 2; % deg/s^2
optimizationParams.minPointing = 10; %deg, floor for pointing constraint
optimizationParams.updateFreq = 10;
optimizationParams.updateStop = 30;
optimizationParams.updateOpt = true;

beta = 0.85;
doPlotting = true; % disable this to not plot results
verboseOutput = true;
tic
[gammaOpt, krOpt, tgoOpt,~,~, optFuelCost, simFuelCost, aTSim,finalPosSim] = getParams(PDIState, planetaryParams, targetState, vehicleParams, optimizationParams, beta, doPlotting, verboseOutput);
toc
% tgoOpt returned in seconds
gammaOpt
krOpt
tgoOpt
gamma2Opt = krOpt/(gammaOpt+2) - 2
optFuelCost
simFuelCost