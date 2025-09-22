function S = runLEMMassOptClosedLoop(x0, cfg)
%RUNLEMMASSOPTCLOSEDLOOP Optimize FP guidance for the LEM and simulate
%
% Inputs
%   x0  [3x1] initial guess for [gamma; kr; tgo]
%   cfg struct with optional overrides for different test configurations from default, fields include:
%       .ub, .lb                               bounds on [gamma, kr, tgo]
%       .fminconOptions                        optimoptions('fmincon', ...)
%       .altitude_km, .lonInitDeg, .latInitDeg
%       .landingLonDeg, .landingLatDeg
%       .inertialVelocity_mps, .flightPathAngle_deg
%       .gMoon, .g0
%       .rfDim, .vfDim, .afDim
%       .tgoMinDim, .deltaTgoDim
%       .massInitDim, .massDryDim, .ispDim
%       .maxThrustDim, .minThrustDim
%       .L_ref
%
% Output
%   S struct with fields used for plotting and reporting:
%       tTraj, stateTraj, aTList, aT_norm, massList
%       refs:   L_ref, T_ref, A_ref, V_ref, M_ref
%       thrust: maxThrustDim, minThrustDim
%       masses: massInitDim, massDryDim
%       opt:    gamma, kr, tgo, tgoVirtual, costEval
%       targets: rfVirt, vfVirt, afVirt
%       params: echo of final dimensional and nondimensional params

    if nargin < 1 || isempty(x0)
        x0 = [1; 6.1; 10]; % Just the x0 I was running with before making this into a function
    end
    if nargin < 2
        cfg = struct(); % Need this to exist, but it can be blank if nothing special
    end
    % ========================================================================
    %% DEFAULTS
    %  ========================================================================
    eps_gamma2 = 1e-3; % This isn't being reached, bounds in the Aineq and bineq are saturated first, going below ~1e-4 there causes singularity
    kr_min = (x0(1) + 2) * (2+eps_gamma2);

    % Tgo for 762.3 seconds is 10.0438629990656 if g = 1.736, 
    % 9.70249498737309 if g=1.62

   % ======================= Bounds ====================
    lb = getfieldOr(cfg, 'lb', [1.0, kr_min, 9.70249498737309]); % Modified version of getfield that also allows a default option if N/A. At bottom of file
    ub = getfieldOr(cfg, 'ub', [1.0, 6.01, 9.70249498737309]);
    % fmincon options
    fminconOptions = getfieldOr(cfg, 'fminconOptions', ...
        optimoptions('fmincon', 'Display', 'none', 'MaxFunctionEvaluations', 1000, ...
        'FiniteDifferenceType', 'central', 'FiniteDifferenceStepSize', 1e-6, ...
        'Algorithm', 'sqp', 'HessianApproximation', 'lbfgs'));

    % ======================= Physical Parameters ====================

    gMoon = getfieldOr(cfg, 'gMoon', 1.62); % Sometimes 1.62 sometimes 1.736 in the paper
    g0    = getfieldOr(cfg, 'g0', 9.81);
    lunarRad_km = 1737.4; %km
    R_moon = lunarRad_km * 1000; %m

    altitude_km        = getfieldOr(cfg, 'altitude_km',        15.24); 
    lonInitDeg         = getfieldOr(cfg, 'lonInitDeg',         41.85); 
    latInitDeg         = getfieldOr(cfg, 'latInitDeg',        -71.6);
    landingLonDeg      = getfieldOr(cfg, 'landingLonDeg',      41.85);
    landingLatDeg      = getfieldOr(cfg, 'landingLatDeg',     -90.0);
    inertialVelocity   = getfieldOr(cfg, 'inertialVelocity_mps', 1698.3);
    flightPathAngleDeg = getfieldOr(cfg, 'flightPathAngle_deg',  0);

    [rDim, vDim] = PDI2MCMF(altitude_km, lonInitDeg, latInitDeg, ...
                                       landingLonDeg, landingLatDeg, ...
                                       inertialVelocity, flightPathAngleDeg);

    [E0, N0, U0] = enuBasis(deg2rad(landingLatDeg),deg2rad(landingLonDeg));

    rfDim = getfieldOr(cfg, 'rfDim', R_moon*U0); % South Pole
    vf_touch = getfieldOr(cfg, 'vf_touch_mps', -1.0);
    vfDim = getfieldOr(cfg, 'vfDim', vf_touch*U0);
    gDim  = [0; 0; -gMoon];
    af_touch = getfieldOr(cfg,'af_touch_mps2', 2*gMoon);
    afDim = getfieldOr(cfg, 'afDim', af_touch * U0);

    
    deltaTgoDim = getfieldOr(cfg, 'deltaTgoDim', 0); % I don't like how much results change with this, but maybe thats fair

    massInitDim = getfieldOr(cfg, 'massInitDim', 15103.0);
    massDryDim  = getfieldOr(cfg, 'massDryDim',  massInitDim - 8248);
    ispDim      = getfieldOr(cfg, 'ispDim',      311);

    maxThrustDim = getfieldOr(cfg, 'maxThrustDim', 45000);
    minThrustDim = getfieldOr(cfg, 'minThrustDim',  4500);

    

    % =========================================================================
    %% NON-DIMENSIONALIZE
    %  ========================================================================
    
    L_ref = getfieldOr(cfg, 'L_ref', 10000);
    T_ref = sqrt(L_ref / gMoon);
    A_ref = gMoon;
    V_ref = L_ref / T_ref;
    M_ref = massInitDim; % This should be changed to the true starting mass,
                         % as a fixed constant once this code is run in a
                         % loop, otherwise M_ref changes each iter
    RMoon_nd = R_moon / L_ref;

    r0      = rDim / L_ref;
    rfStar = rfDim / L_ref;
    afStar = afDim / A_ref;
    v0      = vDim / V_ref;
    vfStar = vfDim / V_ref;

    dtgo   = deltaTgoDim / T_ref;

    %gConst = -(RMoon_nd^2) * r0 / (norm(r0)^3); % Gravity vector at IC, used right now for open loop
    %gConst = -(RMoon_nd^2) * rfStar / (norm(rfStar)^3); % Gravity vector at Landing, used right now for open loop
    %gConst = -(gMoon/A_ref) * U0;
    

    m0     = massInitDim / M_ref;
    mfMin  = massDryDim / M_ref;
    isp    = ispDim * g0 / (A_ref * T_ref);
    %maxThrust = max_thrust_dim / (M_ref * A_ref);
    %minThrust = min_thrust_dim / (M_ref * A_ref);

    % Inequality constraints for fmincon
    Aineq = [-1  0  0;      % gamma >= 0
              2 -1  0;      % kr >= 2*(gamma + 2)
              0  0 -1];     % tgo >= tgoMin
    bineq = [-1e-6; -4-1e-4; -0.001];

    % =========================================================================
    %% OPTIMIZATION
    %  ========================================================================

    obj = @(params) objectiveFunction(params, afStar, rfStar, r0, vfStar, v0, dtgo, RMoon_nd, isp, m0);

    [paramOpt, costEval] = fmincon(obj, x0, Aineq, bineq, [], [], lb, ub, [], fminconOptions);

    gammaOpt = paramOpt(1);
    krOpt    = paramOpt(2);
    tgoOpt   = paramOpt(3);
    tgoOptVirtual = tgoOpt + dtgo;

    gamma1Opt = gammaOpt;
    gamma2Opt = krOpt/(gammaOpt + 2) - 2;

% =========================================================================
%%  TRAJECTORY SIMULATION ODE
%  ========================================================================

    % [rfVirt, vfVirt, afVirt] = computeBeyondTerminationTargeting(r0, v0, gammaOpt, krOpt, rfStar, vfStar, afStar, dtgo, tgoOpt, gConst);
    % [c1, c2] = computeCoefficients(r0, v0, tgoOptVirtual, gamma1Opt, gamma2Opt, afVirt, rfVirt, vfVirt, gConst);


    X0 = [r0; v0; m0];
    odeoptions = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
    
    
    rhs = @(t,x) trajectory(t, x, tgoOpt, dtgo, gammaOpt, gamma1Opt, gamma2Opt, krOpt, ...
                                   afStar, rfStar, vfStar, RMoon_nd, isp);
    N = 1001;
    [tTraj,stateTraj] = ode45(rhs, linspace(0,tgoOpt,N), X0, odeoptions);
    
    % [tTraj, stateTraj] = ode45(@(t, X) trajectory(t, X, gamma1Opt, gamma2Opt, ...
    %     tgoOptVirtual, afVirt, isp, c1, c2, gConst, RMoon_nd, rfStar), tspan, X0, odeoptions);

% =========================================================================
%%  POST-PROCESSING AND ANALYSIS
%  ========================================================================
    aTList = zeros(3, numel(tTraj));
    massList = stateTraj(:, 7)';
    gList = zeros(3,N);
    c1List = zeros(3,N);
    c2List = zeros(3,N);
    rfVirtList = zeros(3,N);
    vfVirtList = zeros(3,N);
    afVirtList = zeros(3,N);

    % for i = 1:length(tTraj)
    %     tgo = tgoOptVirtual - tTraj(i);
    %     aT = afVirt + c1*tgo^gamma1Opt + c2*tgo^gamma2Opt;
    %     aTList(i, :) = aT;
    % end
    % aT_norm = vecnorm(aTList, 2, 2);

    for i = 1:N
        r = stateTraj(i,1:3)';
        v = stateTraj(i,4:6)';
        t = tTraj(i);
        [aTi, gi, coeffsi, rfVirti, vfVirti, afVirti] = trajectoryGetVars(t, r, v, tgoOpt, dtgo, gammaOpt, gamma1Opt, gamma2Opt, krOpt, afStar, rfStar, vfStar, RMoon_nd, isp);
        aTList(:,i) = aTi;
        gList(:,i) = gi';
        c1List(:,i) = coeffsi(1,:);
        c2List(:,i) = coeffsi(2,:);
        rfVirtList(:,i) = rfVirti;
        vfVirtList(:,i) = vfVirti;
        afVirtList(:,i) = afVirti;
    end
    aT_norm = vecnorm(aTList, 2, 1);
% =========================================================================
%%  ORGANIZE OUTPUTS
%  ========================================================================
    
    %Organize outputs to be more usable when run in bulk
    S = struct();

    S.tTraj     = tTraj;
    S.stateTraj = stateTraj;
    S.aTList    = aTList;
    S.aT_norm   = aT_norm;
    S.massList  = massList;

    S.refs = struct('L_ref', L_ref, 'T_ref', T_ref, 'A_ref', A_ref, ...
                    'V_ref', V_ref, 'M_ref', M_ref, 'R_moon', R_moon, 'R_nd', RMoon_nd);

    S.thrust = struct('maxThrustDim', maxThrustDim, 'minThrustDim', minThrustDim);
    S.masses = struct('massInitDim', massInitDim, 'massDryDim', massDryDim);

    S.opt = struct('gamma', gammaOpt, 'kr', krOpt, 'tgo', tgoOpt, ...
                   'tgoVirtual', tgoOptVirtual, 'costEval', costEval, ...
                   'gamma1', gamma1Opt, 'gamma2', gamma2Opt);

    S.targets = struct('rfVirtList', rfVirtList, 'vfVirtList', vfVirtList, 'afVirtList', afVirtList);

    % Echo key dimensional and nondimensional params for reference
    S.params = struct( ...
        'rfDim', rfDim, 'vfDim', vfDim, 'afDim', afDim, 'gDim', gDim, ...
        'rDim', rDim, 'vDim', vDim, 'gMoon', gMoon, 'g0', g0, ...
        'deltaTgoDim', deltaTgoDim, 'ispDim', ispDim, 'mfMin', mfMin,...
        'm0', m0, 'dtgo', dtgo, 'rfStar', rfStar, 'vfStar', vfStar,...
        'afStar', afStar,'lonInitDeg', lonInitDeg, 'latInitDeg', latInitDeg, ...
        'landingLonDeg', landingLonDeg, 'landingLatDeg', landingLatDeg, ...
        'lunarRad_km', lunarRad_km, 'E0', E0, 'N0', N0, 'U0', U0);
    S.optCost = S.opt.costEval * S.refs.M_ref;
    S.odeCost = S.refs.M_ref - (S.massList(end) * S.refs.M_ref);

end


% =========================================================================
%%  SUPPORTING FUNCTIONS
%  ========================================================================

function cost = objectiveFunction(params, afStar, rfStar, r0, vfStar, v0, dtgo, RMoon_nd, isp, m0)
    gamma  = params(1);
    kr     = params(2);
    tgo0   = params(3);
    tgoVirt = tgo0 + dtgo;

    

    gamma1 = gamma;
    gamma2 = kr/(gamma+2) - 2;
    if gamma2 < 0
        cost = 1e10;
        return;
    end
    x0 = [r0(:); v0(:); m0];
    opts = odeset('RelTol',1e-8, 'AbsTol',1e-9);

    rhs = @(t,x) trajectory(t, x, tgo0, dtgo, gamma, gamma1, gamma2, kr, ...
                                   afStar, rfStar, vfStar, RMoon_nd, isp);
    N = 1001;
    [tspan,X] = ode45(rhs, linspace(0,tgo0,N), x0, opts);
    n = numel(tspan);
    aTVec = zeros(3,n);
    aT2 = zeros(1,n);
    % gVec = zeros(3,n);
    % coeffsVec = zeros(2,n);
    % rfVirtVec = zeros(3,n);
    % vfVirtVec = zeros(3,n);
    % afVirtVec = zeros(3,n);

    for i = 1:n
        r = X(i,1:3)';
        v = X(i,4:6)';
        m = X(i,7);
        t = tspan(i);
        [aTi, ~, ~, ~, ~, ~] = trajectoryGetVars(t, r, v, tgo0, dtgo, gamma, gamma1, gamma2, kr, afStar, rfStar, vfStar, RMoon_nd, isp);
        aTVec(:,i) = aTi;
        % gVec(:,i) = gi;
        aT2(i) = dot(aTi,aTi);
    end

    %aT_norm = vecnorm(aTVec, 2, 1);


    % integrate up to real tgo (exclude the beyond segment)
    [~, idx] = min(abs(tspan - tgo0));
    if mod(idx-1, 2) ~= 0
        idx = idx - 1; % Simpson needs even number of intervals
    end
    %J1 = simpsonIntegral(tspan(1:idx), aT_norm(1:idx)); %Use this for fuel optimal
                                                
    J2 = simpsonIntegral(tspan(1:idx), aT2(1:idx)); % This is at^2 optimizing

    cost = J2;
    %For J1 cost, that is total impulse, mf=m0* exp(-I/isp), but can't
    %extract consumed mass from J2
end

function xdot = trajectory(t, x, tgo0, dtgo, gamma, gamma1, gamma2, kr, afStar, rfStar, vfStar, RMoon_nd, isp)
    r = x(1:3);
    v = x(4:6);
    mass = x(7);
    g = -(RMoon_nd^2) * r / (norm(r)^3);
    
    tgo_true = tgo0-t;
    tgoVirt = tgo_true + dtgo;
    
    [rfVirt, vfVirt, afVirt] = computeBeyondTerminationTargeting(r, v, gamma, kr, rfStar, vfStar, afStar, dtgo, tgo0, g);
    [c1, c2] = computeCoefficients(r, v, tgoVirt, gamma1, gamma2, afVirt, rfVirt, vfVirt, g);
    aT = afVirt + c1*tgoVirt^gamma1 + c2*tgoVirt^gamma2;
    
    rdot = v;
    vdot = aT + g;
    F_mag = norm(aT) * mass;

    mdot = -F_mag / (isp);
    
    xdot = [rdot; vdot; mdot];
end

function [aT, g, coeffs, rfVirt, vfVirt, afVirt] = trajectoryGetVars(t, r, v, tgo0, dtgo, gamma, gamma1, gamma2, kr, afStar, rfStar, vfStar, RMoon_nd, isp)
    g = -(RMoon_nd^2) * r / (norm(r)^3);
    
    tgo_true = tgo0-t;
    tgoVirt = tgo_true + dtgo;
    
    [rfVirt, vfVirt, afVirt] = computeBeyondTerminationTargeting(r, v, gamma, kr, rfStar, vfStar, afStar, dtgo, tgo0, g);
    [c1, c2] = computeCoefficients(r, v, tgoVirt, gamma1, gamma2, afVirt, rfVirt, vfVirt, g);
    aT = afVirt + c1*tgoVirt^gamma1 + c2*tgoVirt^gamma2;
    coeffs = [c1;c2];
end

function [c1, c2] = computeCoefficients(r, v, tgo, gamma1, gamma2, af_star, rfStar, vfStar, g)
    phi1_bar = -1/(gamma1 + 1) * tgo^(gamma1 + 1);
    phi2_bar = -1/(gamma2 + 1) * tgo^(gamma2 + 1);

    phi1_hat = (tgo^(gamma1+2))/((gamma1+1)*(gamma1+2));
    phi2_hat = (tgo^(gamma2+2))/((gamma2+1)*(gamma2+2));

    delta = phi1_hat*phi2_bar - phi2_hat*phi1_bar;

    if abs(delta) < 1e-10
        c1 = 0;
        c2 = 0;
        return;
    end

    r_err = r - rfStar + vfStar*tgo - 0.5*(g + af_star)*tgo^2;
    v_err = v - vfStar + (g + af_star)*tgo;

    c1 = ( -phi2_hat * v_err +  phi2_bar * r_err) / delta;
    c2 = (  phi1_hat * v_err -  phi1_bar * r_err) / delta;
end

% function dXdt = trajectory(t, X, gamma1, gamma2, tgoVirt, afVirt, isp, c1, c2, g, R_nd, rfStar)
%     r    = X(1:3);
%     V    = X(4:6); 
%     mass = X(7);
% 
%     %g = -(R_nd^2) * r / (norm(r)^3);
%     g = -(R_nd^2) * rfStar / (norm(rfStar)^3);
% 
%     tgo  = tgoVirt - t;
%     aT = afVirt + c1*tgo^gamma1 + c2*tgo^gamma2;
%     F_mag = norm(aT) * mass;
% 
%     dm_dt = -F_mag / (isp);
%     dXdt = [V; aT + g; dm_dt];
% end

function I = simpsonIntegral(t, y)
    N = length(t) - 1;
    if mod(N, 2) ~= 0
        error('Simpson''s rule requires an even number of intervals');
    end
    h = (t(end) - t(1)) / N;
    I = y(1) + y(end) + 4*sum(y(2:2:end-1)) + 2*sum(y(3:2:end-2));
    I = I * (h/3);
end

function [rfVirtual, vfVirtual, afVirtual, tgoVirtual] = computeBeyondTerminationTargeting(r, v, gamma, kr, rfStar, vfStar, afStar, delta_t, tgo_true, g)
    r = r(:); v = v(:);
    rfStar = rfStar(:); vfStar = vfStar(:); afStar = afStar(:);

    if gamma < 0, error('gamma must be >= 0'); end
    if kr < 2*(gamma + 2), error('kr must be >= 2*(gamma+2)'); end

    tgoVirtual = tgo_true + delta_t;
    if delta_t < 1e-15
        rfVirtual = rfStar; vfVirtual = vfStar; afVirtual = afStar;
        tgoVirtual = tgo_true;
        return;
    end

    gamma1 = gamma;
    gamma2 = kr/(gamma + 2) - 2;

    if abs(gamma1 - gamma2) < 1e-15
        rfVirtual = rfStar; vfVirtual = vfStar; afVirtual = afStar;
        tgoVirtual = tgo_true;
        return;
    end

    tgo = tgoVirtual;

    phi1_bar = -(1/(gamma1 + 1)) * tgo^(gamma1 + 1);
    phi2_bar = -(1/(gamma2 + 1)) * tgo^(gamma2 + 1);
    phi1_hat = tgo^(gamma1 + 2) / ((gamma1 + 1)*(gamma1 + 2));
    phi2_hat = tgo^(gamma2 + 2) / ((gamma2 + 1)*(gamma2 + 2));
    delta = phi1_hat * phi2_bar - phi2_hat * phi1_bar;
    
    
    %fprintf("Value of tgo^(gamma1+1): %.5f, Value of tgo^(gamma2+1): %.5f\n",tgo^(gamma1+1),tgo^(gamma2+1));

    k1r = -phi2_bar / delta;
    k1v = (phi2_bar * tgo + phi2_hat) / delta;
    k1a = -(0.5 * tgo * phi2_bar + phi2_hat) * tgo / delta;

    k2r =  phi1_bar / delta;
    k2v = -(phi1_bar * tgo + phi1_hat) / delta;
    k2a =  (0.5 * tgo * phi1_bar + phi1_hat) * tgo / delta;

    d1 = ((r - 0.5*g*tgo^2)*phi2_bar - (v + g*tgo)*phi2_hat)/delta;
    d2 = -((r - 0.5*g*tgo^2)*phi1_bar - (v + g*tgo)*phi1_hat)/delta;

    phi1_delta = delta_t^gamma1;
    phi2_delta = delta_t^gamma2;
    phi1_bar_delta = -(1/(gamma1 + 1)) * delta_t^(gamma1 + 1);
    phi2_bar_delta = -(1/(gamma2 + 1)) * delta_t^(gamma2 + 1);
    phi1_hat_delta = (1/((gamma1 + 1)*(gamma1 + 2))) * delta_t^(gamma1 + 2);
    phi2_hat_delta = (1/((gamma2 + 1)*(gamma2 + 2))) * delta_t^(gamma2 + 2);

    m11 = k1r * phi1_hat_delta + k2r * phi2_hat_delta + 1;
    m12 = k1v * phi1_hat_delta + k2v * phi2_hat_delta - delta_t;
    m13 = k1a * phi1_hat_delta + k2a * phi2_hat_delta + 0.5 * delta_t^2;

    m21 = k1r * phi1_bar_delta + k2r * phi2_bar_delta;
    m22 = k1v * phi1_bar_delta + k2v * phi2_bar_delta + 1;
    m23 = k1a * phi1_bar_delta + k2a * phi2_bar_delta - delta_t;

    m31 = k1r * phi1_delta + k2r * phi2_delta;
    m32 = k1v * phi1_delta + k2v * phi2_delta;
    m33 = k1a * phi1_delta + k2a * phi2_delta + 1;

    b1 = phi1_hat_delta * d1 + phi2_hat_delta * d2 + 0.5 * g * delta_t^2;
    b2 = phi1_bar_delta * d1 + phi2_bar_delta * d2 - g * delta_t;
    b3 = phi1_delta * d1 + phi2_delta * d2;

    M = zeros(9, 9);
    M(1:3,1:3) = m11*eye(3); M(1:3,4:6) = m12*eye(3); M(1:3,7:9) = m13*eye(3);
    M(4:6,1:3) = m21*eye(3); M(4:6,4:6) = m22*eye(3); M(4:6,7:9) = m23*eye(3);
    M(7:9,1:3) = m31*eye(3); M(7:9,4:6) = m32*eye(3); M(7:9,7:9) = m33*eye(3);

    b_vec = [rfStar - b1; vfStar - b2; afStar - b3];
    %cond(M)
    solution = M \ b_vec;

    rfVirtual = solution(1:3);
    vfVirtual = solution(4:6);
    afVirtual = solution(7:9);
end

function val = getfieldOr(s, name, default)
    if isfield(s, name)
        val = s.(name);
    else
        val = default;
    end
end
function [E, N, U] = enuBasis(lat, lon)
    clat = cos(lat); slat = sin(lat);
    clon = cos(lon); slon = sin(lon);
    U = [clat*clon; clat*slon; slat];
    E = [-slon; clon; 0];
    N = [-slat*clon; -slat*slon; clat];
end