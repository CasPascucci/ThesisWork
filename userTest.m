clear all; close all; clc; format short
%% Test function to see user experience for getParams() function calls
% All values are dimensional
PDIState = struct;
PDIState.altitude_km        = 15.24;
PDIState.lonInitDeg         = 41.85;
PDIState.latInitDeg         = -71.6;
PDIState.inertialVelocity   = 1698.3;
PDIState.flightPathAngleDeg = 0;
PDIState.azimuth            = pi;

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
targetState.landingLonDeg      = 41.85;
targetState.landingLatDeg      = -90.0;
targetState.rfLanding = [0;0;0];
targetState.vfLanding = [0;0;-1];
targetState.afLanding = [0;0;2*planetaryParams.gPlanet];
%targetState.afLanding = [0;0;vehicleParams.maxThrust/(6710)]; % This
                                                               %afStar seems to be too high for a solvable path
targetState.delta_t   = 5; % seconds dim, for btt

optimParams = struct;
optimParams.nodeCount = 997;
optimParams.glideSlopeFinalTheta = 45;




[gammaOpt, krOpt, tgoOpt, optFuelCost, simFuelCost, aTList] = getParams(PDIState, planetaryParams, targetState, vehicleParams, optimParams, 0.65, true);
% Print Values
gammaOpt
krOpt
tgoOpt
optFuelCost
simFuelCost