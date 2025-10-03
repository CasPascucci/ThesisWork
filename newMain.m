clear; clc; close all;
addpath([pwd, '/CoordinateFunctions']);
% Fractional-Polynomial Parameter IC guesses, currently at the optim for
% e-guidance
gamma_test = 1.0;
kr_test = 6.001;
tgo_test = 9.7025;

%% PDI Conditions
altitude_km        = 15.24;
lonInitDeg         = 41.85;
latInitDeg         = -71.6;
landingLonDeg      = 41.85;
landingLatDeg      = -90.0;
inertialVelocity   = 1698.3;
flightPathAngleDeg = 0;

%Convert PDI initial conditions into the MCMF frame of the problem
[r0Dim, v0Dim] = PDI2MCMF(altitude_km, lonInitDeg, latInitDeg, ...
                                   landingLonDeg, landingLatDeg, ...
                                   inertialVelocity, flightPathAngleDeg);

% Gives us the ENU vectors to base the local frame off of
[E0, N0, U0] = enuBasis(deg2rad(landingLatDeg),deg2rad(landingLonDeg));

%% Problem Params


rMoon = 1737.4 * 1000; % m
gMoon = 1.62; % m/s^2
g0 = 9.81; % m/s^2

rfDim = rMoon * U0;
vfTouch = -1;
vfDim = vfTouch * U0; % final velocity -1 m/s straight down at landing site
afTouch = 2*gMoon;
afDim = afTouch * U0; % final accel is 2g upwards at landing site

massInitDim = 15103.0; % kg
dryMassDim = massInitDim - 8248; % kg
ispDim = 311; % s
maxThrustDim = 45000; % N
minThrustDim = 4500; % N

% Struct to pass into optimizer
problemParams = struct;
problemParams.E0 = E0;
problemParams.N0 = N0;
problemParams.U0 = U0;
problemParams.r0Dim = r0Dim; % m
problemParams.v0Dim = v0Dim; % m
problemParams.rMoon = rMoon; % m
problemParams.gMoon = gMoon; % m/s^2
problemParams.g0 = g0; % m/s^2
problemParams.rfDim = rfDim; % m
problemParams.vfDim = vfDim; % m/s
problemParams.afDim = afDim; % m/s^2
problemParams.massInitDim = massInitDim; % kg
problemParams.dryMassDim = dryMassDim; % kg
problemParams.ispDim = ispDim; % s
problemParams.maxThrustDim = maxThrustDim; % N
problemParams.minThrustDim = minThrustDim; % N
problemParams.landingLatDeg = landingLatDeg;
problemParams.landingLonDeg = landingLonDeg;
%% NonDimensionalization

L_ref = 10000; % m
T_ref = sqrt(L_ref/gMoon); % m/s
A_ref = gMoon;
V_ref = L_ref / T_ref;
M_ref = 15103.0; % This should always be the original MInit, keeps scaling consistent as problem evolves


rMoonND  = rMoon / L_ref;
r0ND = r0Dim / L_ref;
v0ND = v0Dim / V_ref;
rfStarND = rfDim / L_ref;
vfStarND = vfDim / V_ref;
afStarND = afDim / A_ref;
gConst = -(rMoonND^2) * rfStarND / (norm(rfStarND)^3); % Constant landing site gravity
m0ND = massInitDim / M_ref;
mMinND = dryMassDim / M_ref;
ispND = ispDim * g0 / (A_ref * T_ref);

% Struct to pass into optimizer
refVals = struct;
refVals.L_ref = L_ref;
refVals.T_ref = T_ref;
refVals.A_ref = A_ref;
refVals.V_ref = V_ref;
refVals.M_ref = M_ref;

% Struct to pass into optimizer
nonDimParams = struct;
nonDimParams.rMoonND = rMoonND;
nonDimParams.r0ND = r0ND;
nonDimParams.v0ND = v0ND;
nonDimParams.rfStarND = rfStarND;
nonDimParams.vfStarND = vfStarND;
nonDimParams.afStarND = afStarND;
nonDimParams.gConst = gConst;
nonDimParams.m0ND = m0ND;
nonDimParams.mMinND = mMinND;
nonDimParams.ispND = ispND;

%% 
paramsX0 = [gamma_test,kr_test,tgo_test]; % Now tgo is nondim

%tStart = tic;

[optParams, optCost] = optimizationLoop(paramsX0, problemParams, nonDimParams, refVals);

%elapsed = toc(tStart);
%fprintf("Finished sim in %.1f s\n", elapsed);
fprintf("Gamma: %.4f, kr: %.4f, tgo: (%.4f) %.4f s \nOpt Estimated Cost: %.4f\n", optParams(1), optParams(2), optParams(3), optParams(3)*T_ref, optCost);

%% Closed Loop

gamma = optParams(1);
kr = optParams(2);
tgo0 = optParams(3);

[tTraj, stateTraj, aTList] = closedLoopSim(gamma, kr, tgo0, problemParams, nonDimParams, refVals);

simCost = M_ref*(stateTraj(1,7) - stateTraj(end,7));
fprintf("Closed Loop Sim Fuel Cost: %.4f kg\n",simCost);


%% 8) Plotting
newPlotting(tTraj, stateTraj, optParams, aTList, refVals, problemParams, nonDimParams)








% waitforbuttonpress;
% for i = 1:7
%     fig = figure(i);
%     waitforbuttonpress;
%     savefig(sprintf('C:/Users/casey/MATLAB Drive/Thesis Research/Media/Figure_%d',i));
% end