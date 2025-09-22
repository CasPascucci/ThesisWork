% Apollo Lunar Lander, Optimized Fractional-Polynomial Guidance 
clear; clc; close all;

%% ========================================================================
%  OPTIMIZATION SETUP
%  ========================================================================

% Initial guess for optimization variables: [gamma, kr, tgo]
x0 = [1; 6.1; 10]; % Initial around E-Guidance
%x0 = [1; 12.1; 10]; % Initial around Apollo Guidance
%x0 = [2; 8.1; 10]; % New initial guess folllowing gamma kr relation

% Variable bounds: [gamma, kr, tgo]
lb = [0.001, 0.001, 0.001];    % Lower bounds - avoid zero/negative values
ub = [5, 25, 20];        % Upper bounds

% Optimization options
fminconOptions = optimoptions('fmincon', ...
    'Display', 'iter-detailed', ...
    'MaxFunctionEvaluations', 1000,...
    'FiniteDifferenceType','central','FiniteDifferenceStepSize',1e-6, ...
    'Algorithm','sqp','HessianApproximation','lbfgs');

%% ========================================================================
%  DIMENSIONAL PHYSICAL PARAMETERS
%  ========================================================================

% Gravitational and environmental constants
%gMars = 3.73;                    % Gravitational acceleration for Mars (m/s²)
gMoon = 1.736;                    % Gravitational acceleration for Moon (m/s²)
g0 = 9.81;                         % Constant for earth, used for isp

%Get PDI IC from Lunar Paper
altitude = 15.24; %km
lonInitDeg = 41.85;
latInitDeg = -71.6;
landingLonDeg = 41.85;
landingLatDeg = -90;
inertialVelocity = 1698.3; %m/s
flightPathAngle = 0; %deg
%lunarRadius = ; %Optional value, but a default is in PDIToLocalCartesian()
[rDim, vDim] = PDIToLocalCartesian(altitude, lonInitDeg, latInitDeg, landingLonDeg, landingLatDeg, inertialVelocity, flightPathAngle);
% ^ Values returned in m and m/s


% Initial and final conditions (dimensional)

rfDim = [0; 0; 0];                % Target position (m)
vfDim = [0; 0; -1.0];             % Target velocity (m/s)
gDim = [0; 0; -gMoon];          % Gravity vector (m/s²)
afDim = [0; 0; -2*gMoon];        % Final acceleration (m/s²)

% Time constraints
tgoMinDim = 1;                   % Minimum time-to-go (s)
deltaTgoDim = 10;                   % deltaT for Beyond-Terimnation Targeting

% Spacecraft mass properties
massInitDim = 15103.0;             % Initial mass (kg)
massDryDim = massInitDim - 8248;   % Dry mass (kg)
ispDim = 311;                      % Specific impulse of Apollo Lunar Module, found, not from paper

% Thrust constraints
maxThrustDim = 45000;           % Maximum thrust (N)
minThrustDim = 4500; % Minimum thrust (N)

%% ========================================================================
%  NON-DIMENSIONALIZATION
%  ========================================================================

% Reference values for non-dimensionalization
L_ref = 10000;                    % Reference length (m)
T_ref = sqrt(L_ref/gMoon);      % Reference time (s)
A_ref = gMoon;          % Reference acceleration (m/s²)
V_ref = L_ref / T_ref;            % Reference velocity (m/s)
M_ref = 15103.0;            % Reference mass (kg)

% Convert to non-dimensional values
r = rDim / L_ref;
rfStar = rfDim / L_ref;
afStar = afDim / A_ref;
v = vDim / V_ref;
vfStar = vfDim / V_ref;
tgoMin = tgoMinDim / T_ref;
dtgo = deltaTgoDim / T_ref;
g = gDim / A_ref;
m0 = massInitDim / M_ref;
mfMin = massDryDim / M_ref;
isp = ispDim * g0/ (A_ref * T_ref);
%maxThrust = max_thrust_dim / (M_ref * A_ref);
%minThrust = min_thrust_dim / (M_ref * A_ref);

%% ========================================================================
%  CONSTRAINT SETUP
%  ========================================================================

% Linear inequality constraints: Ax <= b
linearIneqMatrix = [-1  0  0;    % gamma >= 0
                     2 -1  0;    % 2*gamma - kr <= -4 (kr >= 2*gamma + 4)
                     0  0 -1];   % tgo >= tgo_min
linearIneqVec = [-1e-2; -4-1e-2; -tgoMin];

%% ========================================================================
%  OPTIMIZATION EXECUTION
%  ========================================================================


% Solve the optimization problem
[paramOpt, costEval] = fmincon(...
    @(params) objectiveFunction(params, afStar, g, rfStar, r, vfStar, v, m0, mfMin, isp, dtgo), ...
    x0, linearIneqMatrix, linearIneqVec, [], [], lb, ub, ...
    [], fminconOptions);

% Extract optimal parameters
gammaOpt = paramOpt(1);
krOpt = paramOpt(2);
tgoOpt = paramOpt(3); 
tgoOptVirtual = tgoOpt + dtgo;
% Initial time-to-go (non-dimensional)
tgoOptDim = tgoOpt * T_ref;     % Convert to dimensional time

fprintf('\nOptimal Parameters:\n');
fprintf('  gamma = %.4f\n', gammaOpt);
fprintf('  kr = %.4f\n', krOpt);
fprintf('  tgo = %.4f (%.2f s)\n\n', tgoOpt, tgoOptDim);
fprintf('  Minimum cost = %.3f kg of original %.3f kg\n\n', costEval*M_ref, M_ref); % mass

gamma1Opt = gammaOpt;
gamma2Opt = krOpt/(gammaOpt+2)-2;
fprintf('   gamma1 = %.4f\n', gamma1Opt);
fprintf('   gamma2 = %.4f\n', gamma2Opt);

%% ========================================================================
%  TRAJECTORY SIMULATION
%  ========================================================================

fprintf('\nSimulating optimal trajectory...\n');

X0 = [r; v; m0];
tspan = linspace(0, tgoOpt, 1000); % Ensures trajectory isn't sparse, especially towards beginning


odeoptions = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);%, 'Events',@(t,X) stopWhenTgoEqualsDeltaTgo(t, X, tgoOptVirtual, dtgo));

[rfVirt, vfVirt, afVirt] = computeBeyondTerminationTargeting(r, v, gammaOpt, krOpt, g, rfStar, vfStar, afStar, dtgo, tgoOpt);
[c1,c2] = computeCoefficients(r, v, tgoOptVirtual, gamma1Opt, gamma2Opt, afVirt, g, rfVirt, vfVirt);

[tTraj, stateTraj] = ode45(@(t, X) trajectory(t, X, gamma1Opt, gamma2Opt, ...
    tgoOptVirtual, afVirt, g, isp, c1, c2), tspan, X0, odeoptions);

%% ========================================================================
%  POST-PROCESSING AND ANALYSIS
%  ========================================================================

% Calculate thrust acceleration profile
aTList = zeros(length(tTraj), 3);
massList = stateTraj(:, 7);

%initialRad = stateTraj(1, 1:3).';
%initialVel = stateTraj(1, 4:6).';
%[c1, c2] = computeCoefficients(initialRad, initialVel, tgoOptVirtual, gamma1Opt, gamma2Opt, afVirt, g, rfVirt, vfVirt);

for i = 1:length(tTraj)
    
    tgo = tgoOptVirtual - tTraj(i);
    aT = afVirt + c1*tgo^gamma1Opt + c2*tgo^gamma2Opt;
    aTList(i, :) = aT;
end

aT_norm = vecnorm(aTList, 2, 2);

% Final trajectory analysis
vf_real_norm_dim = sqrt(stateTraj(end,4)^2 + stateTraj(end,5)^2 + ...
                       stateTraj(end,6)^2) * V_ref*sign(stateTraj(end,6));
xf_real_dim = stateTraj(end, 1) * L_ref;
yf_real_dim = stateTraj(end, 2) * L_ref;
zf_real_dim = stateTraj(end, 3) * L_ref;
rf_real_norm_dim = sqrt(stateTraj(end,1)^2 + stateTraj(end,2)^2 + ...
                       stateTraj(end,3)^2) * L_ref;

fprintf('\nFinal Trajectory Results:\n');
fprintf('  Final velocity magnitude: %.2f m/s\n', vf_real_norm_dim);
fprintf('  Final position: [%.2f, %.2f, %.2f] m\n', xf_real_dim, yf_real_dim, zf_real_dim);
fprintf('  Final position magnitude: %.2f m\n', rf_real_norm_dim);
fprintf('  Consumed mass from sim: %.3f kg of original %.3f kg\n', M_ref - M_ref*massList(end), M_ref);

%% ========================================================================
%  VISUALIZATION
%  ========================================================================

% Figure 1: 3D trajectory with thrust acceleration coloring
figure(1); hold on;
scatter3(stateTraj(:,1), stateTraj(:,2), stateTraj(:,3), 20, aT_norm, 'filled');
xlabel('East (Non-dimensional)');
ylabel('North (Non-dimensional)');
zlabel('Up (Non-dimensional)');
title('3D Trajectory with Thrust Acceleration Magnitude');
grid on;
view(45, 45);
%axis equal;
xGrid = -1:0.1:1;
yGrid = -10:5:60;
zGrid = 0;
for i = 1:length(yGrid)
    line(xGrid, yGrid(i)*ones(size(xGrid)), zGrid*ones(size(xGrid)), 'Color',[0.5, 0.5, 0.5], 'LineStyle', '-');
end
for i = 1:length(xGrid)
    line(xGrid(i)*ones(size(yGrid)), yGrid, zGrid*ones(size(yGrid)),'Color',[0.5, 0.5, 0.5], 'LineStyle', '-');
end
colormap(brewermap([], '-RdYlBu'));
c = colorbar;
c.Label.String = 'Thrust Acceleration Magnitude';

% Figure 2: 2D trajectory
figure(2); hold on;
plot(stateTraj(:,2)*L_ref/1000,stateTraj(:,3)*L_ref/1000);
xlabel("North km");
ylabel("Up km");
title("2D Trajectory Plot");
grid on;
yline(0,'Color',[1,0.2,0.2]);
%ylim([0, 16]);
%xlim([0, 600]);
%yticks([0 5 10 15 20]);


% Figure 3: Thrust acceleration components (non-dimensional)
figure(3);
plot(tTraj, aTList(:,1), 'r-', 'LineWidth', 1.5);
hold on;
plot(tTraj, aTList(:,2), 'g-', 'LineWidth', 1.5);
plot(tTraj, aTList(:,3), 'c-', 'LineWidth', 1.5);
plot(tTraj, vecnorm(aTList, 2, 2), '--', 'LineWidth', 2,'Color',[1 0.65 0]);
legend('East Thrust Accel', 'North Thrust Accel', 'Up Thrust Accel', ...
       'Magnitude', 'Location', 'best');
xlabel('Non-dimensional Time');
ylabel('Non-dimensional Acceleration');
title('Thrust Acceleration Profile (Non-dimensional)');
grid on;

% Figure 4: Thrust acceleration components (dimensional)
figure(4);
plot(tTraj * T_ref, aTList(:,1) * A_ref, 'r-', 'LineWidth', 1.5);
hold on;
plot(tTraj * T_ref, aTList(:,2) * A_ref, 'g-', 'LineWidth', 1.5);
plot(tTraj * T_ref, aTList(:,3) * A_ref, 'c-', 'LineWidth', 1.5);
plot(tTraj * T_ref, vecnorm(aTList, 2, 2) * A_ref, '--', 'LineWidth', 2,'Color',[1 0.65 0]);
legend('East Thrust Accel', 'North Thrust Accel', 'Up Thrust Accel', ...
       'Magnitude', 'Location', 'best');
xlabel('Time (s)');
ylabel('Acceleration (m/s²)');
title('Thrust Acceleration Profile (Dimensional)');
grid on;

% Figure 5, Velocity Components
figure(5);
plot(tTraj * T_ref, stateTraj(:,4) * V_ref, 'r-', 'LineWidth', 1.5);
hold on;
plot(tTraj * T_ref, stateTraj(:,5) * V_ref, 'g-', 'LineWidth', 1.5);
plot(tTraj * T_ref, stateTraj(:,6) * V_ref, 'c-', 'LineWidth', 1.5);
plot(tTraj * T_ref, vecnorm([stateTraj(:,4),stateTraj(:,5),stateTraj(:,6)], 2, 2) * V_ref, '--', 'LineWidth', 2,'Color',[1 0.65 0]);
legend('East Vel Component', 'North Vel Component', 'Up Vel Component', ...
       'Magnitude', 'Location', 'best');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Velocity Profile (Dimensional)');
grid on;

% Figure 6, Throttle Profile
figure(6); hold on
plot(tTraj*T_ref, (vecnorm(aTList.*massList,2,2)*A_ref*M_ref)/maxThrustDim);
yline(maxThrustDim/maxThrustDim, 'r--', 'LineWidth', 1, 'DisplayName', 'Max Thrust');
yline(minThrustDim/maxThrustDim, 'r--', 'LineWidth', 1, 'DisplayName', 'Min Thrust');
xlabel('Time');
ylabel('Throttle Fraction');
title('Throttle Profile (Limits just for Display)')
ylim([-0.1,1.1]);


% Figure 7: Mass depletion over time
figure(7);
plot(tTraj, massList, 'c-', 'LineWidth', 2);
hold on;
yline(massInitDim / M_ref, 'r--', 'LineWidth', 1, 'DisplayName', 'Initial Mass');
yline(massDryDim / M_ref, 'g--', 'LineWidth', 1, 'DisplayName', 'Dry Mass');
xlabel('Non-dimensional Time');
ylabel('Non-dimensional Mass');
title('Mass Depletion Profile');
legend('Current Mass', 'Initial Mass', 'Dry Mass', 'Location','northeast');
grid on;

%% ========================================================================
%  SUPPORTING FUNCTIONS
%  ========================================================================

function cost = objectiveFunction(params, afStar, g, rfStar, r, vfStar, v, m, mf_min, isp, dtgo)
    gamma  = params(1);
    kr     = params(2);
    tgo0   = params(3);
    tgoVirt = tgo0 + dtgo;

    % compute gamma1, gamma2
    gamma1 = gamma;
    gamma2 = kr/(gamma+2) - 2;
    if gamma2 < 0
        cost = 1e10;
        return;
    end

    % --- SOLVE c1,c2 ONCE AT THE START ---

    [rfVirt, vfVirt, afVirt] = computeBeyondTerminationTargeting(r, v, gamma, kr, norm(g), rfStar, vfStar, afStar, dtgo, tgo0);
    [c1, c2] = computeCoefficients(...
        r, v, tgoVirt, gamma1, gamma2, afVirt, g, rfVirt, vfVirt);

    N = 1001;
    tspan = linspace(0, tgoVirt, N);
    tgoTrueIdx = find(diff(sign(tspan-tgo0)));
    aTList = zeros(3,length(tspan));
    for i = 1:length(tspan)
        tgo = tgoVirt - tspan(i);
        aT = afVirt + c1*tgo^gamma1 + c2*tgo^gamma2;
        aTList(:, i) = aT;
    end
    
    

 
    aT_norm = vecnorm(aTList,2,1);

    
    gn = -g(3); %Positive value
    mdot = aT_norm / (isp*gn);

    if mod(tgoTrueIdx, 2) ~= 0
        I = simpsonIntegral(tspan(1:tgoTrueIdx), mdot(1:tgoTrueIdx));
    else
        I = simpsonIntegral(tspan(1:tgoTrueIdx-1), mdot(1:tgoTrueIdx-1));
    end

    m_final = m * exp(-I);
    if m_final < mf_min
        deficit = mf_min - m_final;
        cost = 1000*deficit^2 + (m-m_final);
    else
        cost    = m - m_final;
    end
end

%----------------------------------------------------------------------------------------------
function [c1, c2] = computeCoefficients(r, v, tgo, gamma1, gamma2, af_star, g, rf_star, vf_star)
    %phi1 = tgo^gamma1;
    %phi2 = tgo^gamma2;

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

    r_err = r - rf_star + vf_star*tgo - 0.5*(g + af_star)*tgo^2;
    v_err = v - vf_star + (g + af_star)*tgo;

    c1 = -(phi2_hat * v_err) + (phi2_bar * r_err);
    c1 = c1/delta;
    c2 = (phi1_hat * v_err) - (phi1_bar * r_err);
    c2 = c2/delta;

end
%----------------------------------------------------------------------------------------------
function dXdt = trajectory(t, X, gamma1, gamma2, tgoVirt, afVirt, g, isp, c1, c2)
    %Now with BTT
    r    = X(1:3);
    V    = X(4:6);
    mass    = X(7);
    
    tgo  = tgoVirt - t;
    
    % no computeCoefficients here!
    aT = afVirt + c1*tgo^gamma1 + c2*tgo^gamma2;
    norm_aT = norm(aT);
    F_mag = norm_aT * mass;
    
    dm_dt = -F_mag / (isp * (-g(3)));
    dXdt = [V; aT+g; dm_dt];
end

%----------------------------------------------------------------------------------------------
function I = simpsonIntegral(t, y)
    % Simpson's 1/3 rule for numerical integration
    
    N = length(t) - 1;
    if mod(N, 2) ~= 0
        error("Simpson's rule requires an even number of intervals (odd number of points).");
    end
    
    h = (t(end) - t(1)) / N;
    I = y(1) + y(end) + 4*sum(y(2:2:end-1)) + 2*sum(y(3:2:end-2));
    I = I * (h/3);
end

%----------------------------------------------------------------------------------------------
function [rfVirtual, vfVirtual, afVirtual, tgoVirtual] = computeBeyondTerminationTargeting(r, v, gamma, kr, g0, rfStar, vfStar, afStar, delta_t, tgo_true)


    % Ensure column vectors
    r = r(:);
    v = v(:);
    rfStar = rfStar(:);
    vfStar = vfStar(:);
    afStar = afStar(:);
    
    if isscalar(g0)
        g = [0; 0; -g0];
    else
        g = g0(:);
    end
    
    % Check input parameters
    if gamma < 0
        error('Incorrect value for gamma: must be >= 0');
    elseif kr < 2*(gamma + 2)
        error('Incorrect value for kr: must be >= 2*(gamma+2)');
    end
    

    tgoVirtual = tgo_true + delta_t;
    tgo = tgoVirtual;
    %tgo = tgo_true;
    if delta_t < 1e-15
        tgoVirtual = tgo_true;
        rfVirtual = rfStar;
        vfVirtual = vfStar;
        afVirtual = afStar;
        return;
    end
    

    gamma1 = gamma;
    gamma2 = kr/(gamma + 2) - 2;
    

    if abs(gamma1 - gamma2) < 1e-15
        tgoVirtual = tgo_true;
        rfVirtual = rfStar;
        vfVirtual = vfStar;
        afVirtual = afStar;
        return;
    end
    
    % Compute basis functions and their integrals for tgo
    %phi1 = tgo^gamma1;
    %phi2 = tgo^gamma2;
    
    
    % First integrals
    phi1_bar = -(1/(gamma1 + 1)) * tgo^(gamma1 + 1);
    phi2_bar = -(1/(gamma2 + 1)) * tgo^(gamma2 + 1);
    
    % Second integrals
    phi1_hat = tgo^(gamma1 + 2) / ((gamma1 + 1)*(gamma1 + 2));
    phi2_hat = tgo^(gamma2 + 2) / ((gamma2 + 1)*(gamma2 + 2));
    
    delta = phi1_hat * phi2_bar - phi2_hat * phi1_bar;
    
    % Coefficients
    k1r = -phi2_bar / delta;
    k1v = (phi2_bar * tgo + phi2_hat) / delta;
    k1a = -(0.5 * tgo * phi2_bar + phi2_hat) * tgo / delta;
    
    k2r = phi1_bar / delta;
    k2v = -(phi1_bar * tgo + phi1_hat) / delta;
    k2a = (0.5 * tgo * phi1_bar + phi1_hat) * tgo / delta; % Unsure if this should be +-, paper and Fortran do not agree

    d1 = ((r - 0.5*g*tgo^2)*phi2_bar - (v + g*tgo)*phi2_hat)/delta;
    d2 = -((r - 0.5*g*tgo^2)*phi1_bar - (v + g*tgo)*phi1_hat)/delta; % This negative sign is also only in the Fortran, but combining with a + k2a, we reach target
    
    % Now for delta_t
    phi1_delta = delta_t^gamma1;
    phi2_delta = delta_t^gamma2;
    
    
    phi1_bar_delta = -(1/(gamma1 + 1)) * delta_t^(gamma1 + 1);
    phi2_bar_delta = -(1/(gamma2 + 1)) * delta_t^(gamma2 + 1);

    phi1_hat_delta = (1/((gamma1 + 1)*(gamma1 + 2))) * delta_t^(gamma1 + 2);
    phi2_hat_delta = (1/((gamma2 + 1)*(gamma2 + 2))) * delta_t^(gamma2 + 2);
    
    % Build the 3x3 coefficient matrix A
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
    
    % Build the 9x9 matrix M for the linear system
    M = zeros(9, 9);
    
    % Fill the matrix blockwise (3x3 blocks)
    M(1:3,1:3) = m11*eye(3);
    M(1:3,4:6) = m12*eye(3);
    M(1:3,7:9) = m13*eye(3);
    M(4:6,1:3) = m21*eye(3);
    M(4:6,4:6) = m22*eye(3);
    M(4:6,7:9) = m23*eye(3);
    M(7:9,1:3) = m31*eye(3);
    M(7:9,4:6) = m32*eye(3);
    M(7:9,7:9) = m33*eye(3);
  
    % Build the right-hand side vector
    b_vec = [rfStar - b1;
             vfStar - b2;
             afStar - b3];
    
    
    % Solve the full 9x9 system
    solution = M \ b_vec;
    rfVirtual = solution(1:3);
    vfVirtual = solution(4:6);
    afVirtual = solution(7:9);
   
end

function [value, isterminal, direction] = stopWhenTgoEqualsDeltaTgo(t, X, tgoTotal, deltaTgo)
    tgo = tgoTotal - t;
    value = tgo - deltaTgo;  % Trigger when this crosses zero
    isterminal = 1;          % Stop the integration
    direction = -1;          % Only trigger when tgo is decreasing through deltaTgo
end
