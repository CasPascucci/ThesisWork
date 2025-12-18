clear all; close all; clc; format short

%% 1. Initial State Definitions
% PDI (Powered Descent Initiation) State - Dimensional
PDIState = struct;
PDIState.altitude           = 15240;   % m
PDIState.lonInitDeg         = 41.85;   % deg
PDIState.latInitDeg         = -71.59;  % deg
PDIState.inertialVelocity   = 1693.8;  % m/s
PDIState.flightPathAngleDeg = 0;       % deg
PDIState.azimuth            = 180;     % deg

% Planetary Constants (Moon)
planetaryParams = struct;
planetaryParams.rPlanet = 1736.01428 * 1000; % m
planetaryParams.gPlanet = 1.622;             % m/s^2
planetaryParams.gEarth  = 9.81;              % m/s^2

% Vehicle Parameters
vehicleParams = struct;
vehicleParams.massInit  = 15103.0; % kg
vehicleParams.dryMass   = vehicleParams.massInit - 8248; % kg
vehicleParams.isp       = 311;     % s
vehicleParams.maxThrust = 45000;   % N
vehicleParams.minThrust = 4500;    % N

% Target Landing State
targetState = struct;
targetState.landingLonDeg = 41.85; % deg
targetState.landingLatDeg = -90.0; % deg
targetState.rfLanding     = [0; 0; 0]; % m (TOPO)
targetState.vfLanding     = [0; 0; -1]; % m/s
targetState.afLanding     = [0; 0; 2*planetaryParams.gPlanet]; % m/s^2
targetState.delta_t       = 5; % s (BTT parameter, not yet implemented)

%% 2. Optimization Configuration
optimizationParams = struct;
optimizationParams.nodeCount = 301; % Must be odd for Simpson's rule

% Glideslope Constraints
optimizationParams.glideSlopeEnabled    = false;
optimizationParams.glideSlopeFinalTheta = 45;  % deg
optimizationParams.glideSlopeHigh       = 500; % m
optimizationParams.glideSlopeLow        = 250; % m
optimizationParams.glideSlopeCutoff     = 50;  % m

% Pointing Constraints
optimizationParams.pointingEnabled = false;
optimizationParams.maxTiltAccel    = 2;  % deg/s^2
optimizationParams.minPointing     = 10; % deg

% Re-Optimization Settings
optimizationParams.updateOpt  = true; 
optimizationParams.updateFreq = 10;   % s
optimizationParams.updateStop = 30;   % s (Time before landing to stop updates)

% Tolerances
optimizationParams.gamma1eps = 1e-2;
optimizationParams.gamma2eps = 1e-2;

% Divert
targetState.divertEnabled = false; % Also requires reopt to be enabled to truly work
% divertDistances = [0, 50, 100, 150, 200, 300];
% divert1E = divertDistances';
% divert1N = zeros(length(divertDistances),1);
% divert2E = zeros(length(divertDistances),1);
% divert2N = divertDistances';
% divert3E = [divertDistances.*cosd(45)]';
% divert3N = divert3E;
% targetState.divertPoints = [divert1E, divert1N, zeros(length(divertDistances),1);
%                             divert2E, divert2N, zeros(length(divertDistances),1);
%                             divert3E, divert3N, zeros(length(divertDistances),1)];
targetState.divertPoints = [-150, -300, 0]; % m % Single case, above block is multiple divert plots
targetState.altDivert = 2000; % m

%% 3. Execution Flags & Run
beta          = 0.8;  % Weighting: 1.0 = Fuel Optimal, 0.0 = Smoothest Throttle
runSimulation = true;
doPlotting    = true; 
verboseOutput = true;

%tic
[gammaOpt, gamma2Opt, krOpt, tgoOptSec, ~, ~, optFuelCost, simFuelCost, ...
 aTSim, finalPosSim, optHistory, ICstates, exitFlags, problemParams, ...
 nonDimParams, refVals] = getParams(PDIState, planetaryParams, targetState, ...
    vehicleParams, optimizationParams, beta, doPlotting, verboseOutput, false, runSimulation);
%toc

%% 4. Results Display
if optimizationParams.updateOpt
    fprintf("\n--- Initial Solution ---\n")
else
    fprintf('\n--- Final Results ---\n');
end
fprintf('Gamma1:      %.4f\n', gammaOpt);
fprintf('Gamma2:      %.4f\n', gamma2Opt);
fprintf('Kr:          %.4f\n', krOpt);
fprintf('Tgo (sec):   %.2f\n', tgoOptSec);
fprintf('Opt Cost:    %.2f kg\n', optFuelCost);
fprintf('Sim Cost:    %.2f kg\n', simFuelCost);

%% 5. Single Segment Re-Run
% Use this block to isolate and troubleshoot specific re-optimization segments
if optimizationParams.updateOpt
    % segmentIdx = 63;
    % outputSingle = reOptReRun(segmentIdx, ICstates, optHistory, beta, ...
    %     problemParams, nonDimParams, optimizationParams, refVals, ...
    %     targetState.delta_t / refVals.T_ref, verboseOutput);
end