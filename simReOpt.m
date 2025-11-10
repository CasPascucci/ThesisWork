function [tTraj, stateTraj, aTList, flag_thrustGotLimited, optHistory, exitFlags] = simReOpt(gamma0,kr0,tgo0, problemParams, nonDimParams, refVals, delta_t, optimizationParams, betaParam, optimParams, verboseOutput)

% Original IC's
r0 = nonDimParams.r0ND;
v0 = nonDimParams.v0ND;
m0 = nonDimParams.m0ND;
rfStar = nonDimParams.rfStarND;
vfStar = nonDimParams.vfStarND;
afStar = nonDimParams.afStarND;
isp = nonDimParams.ispND;
rMoonND = nonDimParams.rMoonND;
gConst = nonDimParams.gConst;

X0 = [r0; v0; m0];
t_curr = 0;
odeoptions = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

% Setup Results
tTraj = [];
stateTraj = [];
aTList = [];
flag_thrustGotLimited = false;
optHistory = [];
exitFlags = [];

%Start main loop
gamma = gamma0;
kr = kr0;
tgo = tgo0;












end