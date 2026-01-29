clear all; clc; close all;
addpath([pwd, '/CoordinateFunctions']);

fixedGamma1 = 0.07856;
fixedGamma2 = 0.7570;
tgoRange = [567,800];
% IF console prints message stating that too many constraints are active,
% adjust tgoRange to a safer region, constraints at beginning of
% activeTgoSweep mean the range should be moved to higher values, and vice
% versa for constraints at end of activeTgoSweep

PDIState = struct('altitude', 15240, ...
                  'lonInitDeg', 41.85, ...
                  'latInitDeg', -71.59, ...
                  'inertialVelocity', 1693.8, ...
                  'flightPathAngleDeg', 0, ...
                  'azimuth', 180);
planetaryParams = struct('rPlanet', 1736.01428 * 1000, ...
                         'gPlanet', 1.622, ...
                         'gEarth', 9.81);
vehicleParams = struct('massInit', 15103.0, ...
                       'dryMass', 6855, ...
                       'isp', 311, ...
                       'maxThrust', 45000, ...
                       'minThrust', 4500);
targetState = struct('landingLonDeg', 41.85, ...
                     'landingLatDeg', -90.00, ...
                     'rfLanding', [0;0;0], ...
                     'vfLanding', [0;0;-1], ...
                     'afLanding', [0;0;3.244], ...
                     'delta_t', 5, ...
                     'divertEnabled', false, ...
                     'altDivert', 1000, ...
                     'divertPoints', [0,0,0]);

optParams = struct;
optParams.nodeCount = 301; 

optParams.glideSlopeEnabled = false;
optParams.glideSlopeFinalTheta = 45;
optParams.glideSlopeHigh = 500;
optParams.glideSlopeLow = 250;
optParams.freeGlideNodes = 1;

optParams.pointingEnabled = false;
optParams.maxTiltAccel = 2;
optParams.minPointing = 10;

L_ref = 10000;
T_ref = sqrt(L_ref/planetaryParams.gPlanet);
refVals = struct('L_ref', L_ref, 'T_ref', T_ref, 'A_ref', planetaryParams.gPlanet, ...
    'V_ref', L_ref/T_ref, 'M_ref', vehicleParams.massInit);

[r0Dim, v0Dim] = PDI2MCMF(PDIState.altitude/1000, PDIState.lonInitDeg, PDIState.latInitDeg, ...
    targetState.landingLonDeg, targetState.landingLatDeg, PDIState.inertialVelocity, ...
    PDIState.flightPathAngleDeg, PDIState.azimuth*pi/180, planetaryParams.rPlanet);

rfDim = refVals.L_ref * ENU2MCMF(targetState.rfLanding/10000, targetState.landingLatDeg, targetState.landingLonDeg, true);
vfDim = ENU2MCMF(targetState.vfLanding, targetState.landingLatDeg, targetState.landingLonDeg, false);
afDim = ENU2MCMF(targetState.afLanding, targetState.landingLatDeg, targetState.landingLonDeg, false);


nonDim = struct();
nonDim.r0ND = r0Dim/L_ref;
nonDim.v0ND = v0Dim/refVals.V_ref;
nonDim.rfStarND = rfDim/L_ref;
nonDim.vfStarND = vfDim/refVals.V_ref;
nonDim.afStarND = afDim/refVals.A_ref;
nonDim.gConst = -(planetaryParams.rPlanet/L_ref)^2 * nonDim.rfStarND / (norm(nonDim.rfStarND)^3);
nonDim.m0ND = 1.0;
nonDim.ispND = vehicleParams.isp * 9.81 / refVals.V_ref;
nonDim.maxThrustND = vehicleParams.maxThrust / (vehicleParams.massInit * refVals.A_ref);
nonDim.minThrustND = vehicleParams.minThrust / (vehicleParams.massInit * refVals.A_ref);

fprintf('Tgo Sweep, Fixed Params: G1=%.2f, G2=%.2f\n', fixedGamma1, fixedGamma2);

numPoints = 100;

tgoVec = linspace(tgoRange(1),tgoRange(2),numPoints);
fuelTgoSweep = zeros(numPoints,1);
activeTgoSweep = zeros(numPoints,1);

for idx = 1: numPoints
    [fuelCost, nActive] = evaluateTraj(fixedGamma1,fixedGamma2, tgoVec(idx), nonDim, refVals, optParams);
    fuelTgoSweep(idx) = fuelCost;
    activeTgoSweep(idx) = nActive;
end

% Check if all tests are within constraints
maxConstraints = max(activeTgoSweep, [], 'all', 'omitnan');
if maxConstraints >= 2
    fprintf("Thrust Saturated at more than one point on at least one trajectory\n")
end

% Find max and minimum fuel, and slope approximation between them
[maxFuel, maxFuelIdx] = max(fuelTgoSweep);
[minFuel, minFuelIdx] = min(fuelTgoSweep);
maxFuelTime = tgoVec(maxFuelIdx);
minFuelTime = tgoVec(minFuelIdx);
fprintf("Max fuel cost of %.1f @ tgo = %.1f\n", maxFuel, maxFuelTime);
fprintf("Min fuel cost of %.1f @ tgo = %.1f\n", minFuel, minFuelTime);
fprintf("Approx slope: %.2f", (maxFuel-minFuel)/(maxFuelTime-minFuelTime));

%% Saving
folderName = sprintf('Tgo Sweep/G1_%.2f_G2_%.2f', fixedGamma1, fixedGamma2);
if ~exist(folderName,"dir")
    mkdir(folderName);
end

f1 = figure('Name','Fuel Consumption Sensitivty to Tgo Sweep'); hold on;
plot(tgoVec, fuelTgoSweep, 'b-', 'LineWidth',2);
xlabel('Time of Flight of Trajectory'); ylabel('Fuel Consumption from Optimized Trajectory');
title('Fuel Consumption Sensitivty to $t_{go}$ Sweep','Interpreter','latex');
subtitle(sprintf("$\\gamma_1$ = %.3f, $\\gamma_2$ = %.3f", fixedGamma1, fixedGamma2), 'Interpreter','latex');
f1axes = f1.findobj("Type", "axes");
f1axes.FontSize = 20;



function [fuelKg, activeCount] = evaluateTraj(gamma1, gamma2, tgoSec, nonDim, refVals, optParams)
    tgoND = tgoSec / refVals.T_ref;
    try
        [c1, c2] = calculateCoeffs(nonDim.r0ND, nonDim.v0ND, tgoND, gamma1, gamma2, nonDim.afStarND, nonDim.rfStarND, nonDim.vfStarND, nonDim.gConst);
    catch
        fuelKg = NaN; activeCount = NaN; return;
    end
    tgospan = linspace(0, tgoND, optParams.nodeCount);
    aT = nonDim.afStarND + c1.*(tgospan.^gamma1) + c2.*(tgospan.^gamma2);
    aTmag = vecnorm(aT, 2, 1);
    Q = cumtrapz(tgospan, aTmag ./ nonDim.ispND);
    Q = Q(end) - Q;
    
    mCurrent = nonDim.m0ND * exp(-Q);
    thrustMag = mCurrent .* aTmag; 
    fuelKg = refVals.M_ref * nonDim.m0ND * (1 - mCurrent(1)/nonDim.m0ND);
    tol = 1e-4;
    activeNodes = sum((thrustMag - nonDim.maxThrustND) > tol) + sum((nonDim.minThrustND - thrustMag) > tol);
    activeCount = activeNodes;

end