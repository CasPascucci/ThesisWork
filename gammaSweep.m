clear all; clc; close all;
addpath([pwd, '/CoordinateFunctions']);

fixedTgo = 594.69;
gammaRange = [0.01, 1.5];
gridResolution = 100;

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

fprintf('Gamma Sweep, Fixed Params: Tgo=%.0fs\n', fixedTgo);


g1Vec = linspace(gammaRange(1), gammaRange(2), gridResolution);
g2Vec = linspace(gammaRange(1), gammaRange(2), gridResolution);
[G1, G2] = meshgrid(g1Vec, g2Vec);

fuelGammaSweep = zeros(gridResolution,gridResolution);
activeGammaSweep = zeros(gridResolution,gridResolution);

for idx = 1:numel(G1)
    if abs(G1(idx) - G2(idx)) < 1e-4
        fuelGammaSweep(idx) = NaN;
        activeGammaSweep(idx) = NaN;
        continue;
    end
    [fuelCost, nActive] = evaluateTraj(G1(idx),G2(idx), fixedTgo, nonDim, refVals, optParams);
    fuelGammaSweep(idx) = fuelCost;
    activeGammaSweep(idx) = nActive;
end

maxConstraints = max(activeGammaSweep, [], 'all', 'omitnan');
if maxConstraints >= 2
    fprintf("Thrust Saturated at more than one point on at least one trajectory\n")
end

maxFuel = max(fuelGammaSweep, [], 'all', 'omitnan');
minFuel = min(fuelGammaSweep, [], 'all', 'omitnan');
validMask = ~isnan(fuelGammaSweep);
if any(validMask, 'all')
    [minRow, minCol] = find(fuelGammaSweep == minFuel, 1);
    [maxRow, maxCol] = find(fuelGammaSweep == maxFuel, 1);
    bestG1min = G1(minRow, minCol);
    bestG2min = G2(minRow, minCol);
    bestG1max = G1(maxRow, maxCol);
    bestG2max = G2(maxRow, maxCol);
    fprintf("Min fuel cost of %.1f @ G1=%.2f, G2=%.2f\n", minFuel, bestG1min, bestG2min);
    fprintf("Max fuel cost of %.1f @ G1=%.2f, G2=%.2f\n", maxFuel, bestG1max, bestG2max);
else
    fprintf("No valid solutions found.\n");
end

folderName = sprintf('Gamma Sweep/Tgo_%d', round(fixedTgo));
if ~exist(folderName,"dir")
    mkdir(folderName);
end

% Figure 1: Fuel Contour
f1 = figure('Name','Fuel Consumption Sensitivity to Gamma'); hold on;
contourf(G1, G2, fuelGammaSweep, 50); colorbar;
axis equal;
xlabel('$\gamma_1$','Interpreter','latex'); 
ylabel('$\gamma_2$','Interpreter','latex');
title('Fuel Consumption Sensitivity to Gamma','Interpreter','latex');
subtitle(sprintf("$t_{go}$ = %d s", round(fixedTgo)), 'Interpreter','latex');
f1axes = f1.findobj("Type", "axes");
f1axes.FontSize = 20;
saveas(f1, fullfile(folderName, 'Fuel_Contour.png'));
saveas(f1, fullfile(folderName, 'Fuel_Contour.fig'));

% Figure 2: Constraints Surface
f2 = figure('Name','Active Constraints Map'); hold on;
surf(G1, G2, activeGammaSweep); view(0, 90); colorbar; shading interp;
axis equal;
xlabel('$\gamma_1$','Interpreter','latex'); 
ylabel('$\gamma_2$','Interpreter','latex');
title('Active Constraints Map','Interpreter','latex');
subtitle(sprintf("$t_{go}$ = %d s", round(fixedTgo)), 'Interpreter','latex');
f2axes = f2.findobj("Type", "axes");
f2axes.FontSize = 20;
saveas(f2, fullfile(folderName, 'Constraints_Map.png'));
saveas(f2, fullfile(folderName, 'Constraints_Map.fig'));



T_Gamma = table(G1(:), G2(:), fuelGammaSweep(:), activeGammaSweep(:), ...
    'VariableNames', {'Gamma1', 'Gamma2', 'FuelCost_kg', 'ActiveConstraints'});
writetable(T_Gamma, fullfile(folderName, 'Gamma_Sweep_Data.csv'));

fprintf('Results saved to: %s\n', folderName);
%% Functions
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