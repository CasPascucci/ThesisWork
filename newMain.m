clear; clc; close all;
addpath([pwd, '/CoordinateFunctions']);
% Fractional-Polynomial Parameter IC guesses, currently at the optim for
% e-guidance
gamma_test = 1.0;
kr_test = 6.1;
tgo_test = 9.70249498737309;

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
problemParams.massInitDim = massInitDim; % kg
problemParams.dryMassDim = dryMassDim; % kg
problemParams.ispDim = ispDim; % s
problemParams.maxThrustDim = maxThrustDim; % N
problemParams.minThrustDim = minThrustDim; % N
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
%gConst = -(rMoonND^2) * r0ND / (norm(rfStarND)^3); % Constant PDI site gravity
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
paramsX0 = [gamma_test,kr_test,(tgo_test / T_ref)]; % Now tgo is nondim

tStart = tic;

[optParams, optCost] = optimizationLoop(paramsX0, problemParams, nonDimParams, refVals);

elapsed = toc(tStart);
fprintf("Finished sim in %.1f s\n", elapsed);
fprintf("Gamma: %.4f, kr: %.4f, tgo: (%.4f) %.4f s \nCost: %.4f\n", optParams(1), optParams(2), optParams(3), optParams(3)*T_ref, optCost);

%% Updated up to here only

% Build a summary table
idx_ok = find([Results.ok]);
n_ok = numel(idx_ok);

fprintf("%d tests were OK out of the total %d tests\n", n_ok, numTests);

%Pre allocate
gamma   = nan(numTests,1);
kr      = nan(numTests,1);
tgo     = nan(numTests,1);
tgoVirt = nan(numTests,1);
cost    = nan(numTests,1);
prop_kg = nan(numTests,1);
errMsg  = strings(numTests,1);

for i = 1:numTests
    if Results(i).ok %if ok, but so far all tests have been succesful, if it continues to do so then this if else can be removed
        S = Results(i).S;
        gamma(i)   = S.opt.gamma;
        kr(i)      = S.opt.kr;
        tgo(i)     = S.opt.tgo;
        tgoVirt(i) = S.opt.tgoVirtual;
        cost(i)    = S.opt.costEval;

        % propellant used in kg
        M_ref = S.refs.M_ref;            % equals massInitDim in current setup
        massInitDim = S.masses.massInitDim;
        mf_kg = S.stateTraj(end,7) * M_ref;
        prop_kg(i) = massInitDim - mf_kg;
    else
        errMsg(i) = Results(i).err;
    end
end

Summary = table(gamma, kr, tgo, tgoVirt, cost, prop_kg, errMsg); %Table helps for easy mass inspection
Summary.x0_gamma = X0_all(:,1);
Summary.x0_kr    = X0_all(:,2);
Summary.x0_tgo   = X0_all(:,3);

% Keep only successful rows for ranking, probably not necessary
SummaryOK = Summary( [Results.ok].', : );
SummaryOK = sortrows(SummaryOK, 'prop_kg');   % or 'cost'

% 7) Save everything
% ts = datestr(now, 'yyyymmdd_HHMMSS');
% md = matlabdrive;
% folder = "/Thesis Research/Batch Runs";
% fname = ['LEM_batch_' ts '.mat'];
% fullPath = fullfile(md, folder, fname);
% save(fullPath, 'Summary', 'SummaryOK', 'Results', 'X0_all', 'cfg', 'elapsed');
% fprintf('Saved batch to %s\n', fname);

%% 8) Plotting
% Run this section after loading desired results

% Find index of min propellant in the OK table
[~, minRow] = min(SummaryOK.prop_kg);
best = SummaryOK(minRow,:);

% Find this row in the original Summary, same idx as Results
% Probably unnecessary step if all tests are OK
bestIdx = find(Summary.gamma == best.gamma & ...
               Summary.kr    == best.kr & ...
               Summary.tgo   == best.tgo, 1, 'first');

if ~isempty(bestIdx) && Results(bestIdx).ok % Found a result and it actually ran
    plotLEMMassOpt(Results(bestIdx).S);
    fprintf('Plotted best run: gamma=%.4f, kr=%.4f, tgo=%.4f, prop=%.1f kg\n', ...
        best.gamma, best.kr, best.tgo, best.prop_kg);
else
    warning('Best run not found in Results.');
end


[minGamma, maxGamma] = bounds(SummaryOK.gamma);
rangeGamma = maxGamma - minGamma;
[minkr, maxkr] = bounds(SummaryOK.kr);
rangekr = maxkr - minkr;

% Figure 8: Gamma vs Kr spread for batch runs
figure(); hold on;
scatter(SummaryOK.gamma, SummaryOK.kr,10)
xscale('log');
title('Gamma vs Kr');
subtitle(sprintf('Gamma Spread: %.3e \n Kr Spread: %.3e', rangeGamma, rangekr));
xlabel('Gamma');
ylabel('Kr');

%Figure 9: Case Number vs Dimensional Tgo
figure();
plot(1:numTests, Summary.tgo*S.refs.T_ref,'.','MarkerSize',30);
title('Case Number vs Tgo (Dimensional)');
xlabel('Case Number');
ylabel('Tgo (seconds)');
grid;


% waitforbuttonpress;
% for i = 1:7
%     fig = figure(i);
%     waitforbuttonpress;
%     savefig(sprintf('C:/Users/casey/MATLAB Drive/Thesis Research/Media/8_15_25_Figures/Gamma 0.5/Figure_%d',i));
% end