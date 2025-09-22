% Setup using Eq43 formulation from the FP-PDG paper, this is a backup
% without the fixes described in the 7/22 Meeting
clear; clc; close all;

%% ========================================================================
%  OPTIMIZATION SETUP
%% ========================================================================

% Initial guess for optimization variables: [gamma, kr, tgo]
x0 = [1; 6; 1]; % Initial around E-Guidance
%x0 = [1; 12; 1.5]; % Initial around Apollo Guidance
%x0 = [2; 8; 1]; % New initial guess folllowing gamma kr relation

% Variable bounds: [gamma, kr, tgo]
lb = [0.1, 4, 0.5];    % Lower bounds - avoid zero/negative values
ub = [5, 25, 5];     % Upper bounds

% Optimization options
options = optimoptions('fmincon', ...
    'Display', 'iter-detailed', ...
    'MaxFunctionEvaluations', 1000,...
    'FiniteDifferenceType','central','FiniteDifferenceStepSize',1e-6);

%% ========================================================================
%  DIMENSIONAL PHYSICAL PARAMETERS
%% ========================================================================

% Gravitational and environmental constants
g_const = 3.73;                    % Gravitational acceleration (m/s²)

% Initial and final conditions (dimensional)
r_dim = [20000; 3000; 5000];       % Initial position (m)
rf_dim = [0; 0; 0];                % Target position (m)
v_dim = [-500; 0; -200];           % Initial velocity (m/s)
vf_dim = [0; 0; -0.1];             % Target velocity (m/s)
g_dim = [0; 0; -g_const];          % Gravity vector (m/s²)
af_dim = [0; 0; 2*g_const];        % Final acceleration (m/s²)

% Time constraints
tgo_min_dim = 1;                   % Minimum time-to-go (s)

% Spacecraft mass properties
mass_init_dim = 62000;             % Initial mass (kg)
mass_dry_dim = 18000;              % Dry mass (kg)
isp_dim = 330;                     % Specific impulse (s)

% Thrust constraints
max_thrust_dim = 800000;           % Maximum thrust (N)
min_thrust_dim = 0.3 * max_thrust_dim; % Minimum thrust (N)

%% ========================================================================
%  NON-DIMENSIONALIZATION
%% ========================================================================

% Reference values for non-dimensionalization
L_ref = 1000;                      % Reference length (m)
T_ref = 100;                       % Reference time (s)
A_ref = L_ref / T_ref^2;          % Reference acceleration (m/s²)
V_ref = L_ref / T_ref;            % Reference velocity (m/s)
M_ref = mass_init_dim;            % Reference mass (kg)

% Convert to non-dimensional values
r = r_dim / L_ref;
rfStar = rf_dim / L_ref;
afStar = af_dim / A_ref;
v = v_dim / V_ref;
vfStar = vf_dim / V_ref;
tgo_min = tgo_min_dim / T_ref;
g = g_dim / A_ref;
m0 = mass_init_dim / M_ref;
mf_min = mass_dry_dim / M_ref;
isp = isp_dim / T_ref;
maxThrust = max_thrust_dim / (M_ref * A_ref);
minThrust = min_thrust_dim / (M_ref * A_ref);

%% ========================================================================
%  CONSTRAINT SETUP
%% ========================================================================

% Nonlinear constraints function handle
nonLinearCons = @(x) nonlinearConstraints(x, afStar, g, rfStar, r, v, ...
    vfStar, m0, maxThrust, minThrust, isp);

% Linear inequality constraints: Ax <= b
linearIneqMatrix = [-1  0  0;    % gamma >= 0
                       2 -1  0;    % 2*gamma - kr <= -4 (kr >= 2*gamma + 4)
                       0  0 -1];   % tgo >= tgo_min
linearIneqVec = [0; -4; -tgo_min];

%% ========================================================================
%  OPTIMIZATION EXECUTION
%% ========================================================================

fprintf('Starting trajectory optimization...\n');

% Solve the optimization problem
[x_opt, fval] = fmincon(...
    @(params) objective(params, afStar, g, rfStar, r, vfStar, v, m0, ...
                       maxThrust, minThrust, isp), ...
    x0, linearIneqMatrix, linearIneqVec, [], [], lb, ub, ...
    nonLinearCons, options);

% Extract optimal parameters
gammaOpt = x_opt(1);
krOpt = x_opt(2);
tgo_opt = x_opt(3);                % Initial time-to-go (non-dimensional)
tgo_opt_dim = tgo_opt * T_ref;     % Convert to dimensional time

fprintf('\nOptimal Parameters:\n');
fprintf('  gamma = %.4f\n', gammaOpt);
fprintf('  kr = %.4f\n', krOpt);
fprintf('  tgo = %.4f (%.2f s)\n', tgo_opt, tgo_opt_dim);
fprintf('  Minimum cost = %.6f\n', fval);

%% ========================================================================
%  TRAJECTORY SIMULATION
%% ========================================================================

fprintf('\nSimulating optimal trajectory...\n');

% Initial state vector: [position; velocity; mass]
X0 = [r; v; m0];
tspan = [0, tgo_opt];

% ODE solver options for high accuracy
odeoptions = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'MaxStep', 0.001);

% Simulate the trajectory
[t_traj, state_traj] = ode45(@(t, X) trajectory(t, X, gammaOpt, krOpt, ...
    tgo_opt, afStar, g, rfStar, vfStar, maxThrust, minThrust, isp), ...
    tspan, X0, odeoptions);

%% ========================================================================
%  POST-PROCESSING AND ANALYSIS
%% ========================================================================

% Calculate thrust acceleration profile
aT_list = zeros(length(t_traj), 3);
mass_list = state_traj(:, 7);

for i = 1:length(t_traj)
    current_r = state_traj(i, 1:3).';
    current_v = state_traj(i, 4:6).';
    tgo = max((tgo_opt - t_traj(i)), 0.001);
    
    aT = compute_aT(gammaOpt, krOpt, tgo, afStar, g, rfStar, ...
                   current_r, vfStar, current_v).';
    norm_aT = norm(aT);
    F_mag = norm_aT * mass_list(i);
    
    % Apply thrust limits
    if F_mag > maxThrust
        F_mag = maxThrust;
        aT = (aT / norm_aT) * (maxThrust / mass_list(i));
    elseif F_mag < minThrust
        F_mag = minThrust;
        aT = (aT / norm_aT) * (minThrust / mass_list(i));
    end
    
    aT_list(i, :) = aT;
end

aT_norm = vecnorm(aT_list, 2, 2);

% Final trajectory analysis
vf_real_norm_dim = sqrt(state_traj(end,4)^2 + state_traj(end,5)^2 + ...
                       state_traj(end,6)^2) * V_ref;
zf_real_dim = state_traj(end, 3) * L_ref;
rf_real_norm_dim = sqrt(state_traj(end,1)^2 + state_traj(end,2)^2 + ...
                       state_traj(end,3)^2) * L_ref;

fprintf('\nFinal Trajectory Results:\n');
fprintf('  Final velocity magnitude: %.2f m/s\n', vf_real_norm_dim);
fprintf('  Final altitude: %.2f m\n', zf_real_dim);
fprintf('  Final position magnitude: %.2f m\n', rf_real_norm_dim);

%% ========================================================================
%  VISUALIZATION
%% ========================================================================

% Figure 1: 3D trajectory with thrust acceleration coloring
figure(1);
scatter3(state_traj(:,1), state_traj(:,2), state_traj(:,3), 20, aT_norm, 'filled');
xlabel('X (Non-dimensional)');
ylabel('Y (Non-dimensional)');
zlabel('Z (Non-dimensional)');
title('3D Trajectory with Thrust Acceleration Magnitude');
grid on;
view(45, 45);
axis equal;
colormap(brewermap([], '-RdYlBu'));
c = colorbar;
c.Label.String = 'Thrust Acceleration Magnitude';

% Figure 2: Thrust acceleration components (non-dimensional)
figure(2);
plot(t_traj, aT_list(:,1), 'r-', 'LineWidth', 1.5);
hold on;
plot(t_traj, aT_list(:,2), 'g-', 'LineWidth', 1.5);
plot(t_traj, aT_list(:,3), 'c-', 'LineWidth', 1.5);
plot(t_traj, vecnorm(aT_list, 2, 2), '--', 'LineWidth', 2,'Color',[1 0.65 0]);
plot(t_traj, maxThrust ./ mass_list, 'm--', 'LineWidth', 1);
plot(t_traj, minThrust ./ mass_list, 'm--', 'LineWidth', 1);
legend('X Thrust Accel', 'Y Thrust Accel', 'Z Thrust Accel', ...
       'Magnitude', 'Max Limit', 'Min Limit', 'Location', 'best');
xlabel('Non-dimensional Time');
ylabel('Non-dimensional Acceleration');
title('Thrust Acceleration Profile (Non-dimensional)');
grid on;

% Figure 3: Thrust acceleration components (dimensional)
figure(3);
plot(t_traj * T_ref, aT_list(:,1) * A_ref, 'r-', 'LineWidth', 1.5);
hold on;
plot(t_traj * T_ref, aT_list(:,2) * A_ref, 'g-', 'LineWidth', 1.5);
plot(t_traj * T_ref, aT_list(:,3) * A_ref, 'c-', 'LineWidth', 1.5);
plot(t_traj * T_ref, vecnorm(aT_list, 2, 2) * A_ref, '--', 'LineWidth', 2,'Color',[1 0.65 0]);
plot(t_traj * T_ref, maxThrust ./ mass_list * A_ref, 'm--', 'LineWidth', 1);
plot(t_traj * T_ref, minThrust ./ mass_list * A_ref, 'm--', 'LineWidth', 1);
legend('X Thrust Accel', 'Y Thrust Accel', 'Z Thrust Accel', ...
       'Magnitude', 'Max Limit', 'Min Limit', 'Location', 'best');
xlabel('Time (s)');
ylabel('Acceleration (m/s²)');
title('Thrust Acceleration Profile (Dimensional)');
grid on;

% Figure 4: Mass depletion over time
figure(4);
plot(t_traj, mass_list, 'c-', 'LineWidth', 2);
hold on;
yline(mass_init_dim / M_ref, 'r--', 'LineWidth', 1, 'DisplayName', 'Initial Mass');
yline(mass_dry_dim / M_ref, 'g--', 'LineWidth', 1, 'DisplayName', 'Dry Mass');
xlabel('Non-dimensional Time');
ylabel('Non-dimensional Mass');
title('Mass Depletion Profile');
legend('Current Mass', 'Initial Mass', 'Dry Mass', 'Location','northeast');
grid on;

%% ========================================================================
%  SUPPORTING FUNCTIONS
%% ========================================================================

function cost = objective(params, af_star, g, rf_star, r, vf_star, v, m, ...
                         max_thrust, min_thrust, isp)
    % Objective function: minimize fuel consumption (integrated thrust)
    
    gamma = params(1);
    kr = params(2);
    tgo = params(3);
    
    N = 300;  % Number of integration points
    tspan = linspace(0, tgo, N);
    odeoptions = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
    X0 = [r; v; m];

    % Simulate current trajectory
    [T, X] = ode45(@(t_nd, X_nd) trajectory(t_nd, X_nd, gamma, kr, tgo, ...
                                           af_star, g, rf_star, vf_star, ...
                                           max_thrust, min_thrust, isp), ...
                   tspan, X0, odeoptions);
    
    % Calculate thrust magnitude profile
    thrustProfile = zeros(1, N);
    for i = 1:N
        r_i = X(i, 1:3);
        v_i = X(i, 4:6);
        tau = max(tgo - T(i), 0.001);
        aT_i = compute_aT(gamma, kr, tau, af_star, g, rf_star, r_i, vf_star, v_i);
        thrustProfile(i) = norm(aT_i) * X(i, 7);  % Force = acceleration × mass
    end
    
    % Ensure even number of intervals for Simpson's rule
    if mod(length(T) - 1, 2) == 1
        T = T(1:end-1);
        thrustProfile = thrustProfile(1:end-1);
    end
    
    % Integrate thrust over time to get fuel consumption
    cost = simpsonIntegral(T, thrustProfile);
end

function aT = compute_aT(gamma, kr, tgo, af_star, g, rf_star, r, V_star, V)
    % Compute commanded thrust acceleration using optimal control law
    
    term1 = gamma * (((kr / (2*(gamma + 2))) - 1) * af_star);
    term2 = ((gamma * kr / (2*(gamma + 2))) - gamma - 1) * g;
    term3 = (gamma + 1) / tgo * (1 - (kr / (gamma + 2))) * (V_star - V);
    term4 = (kr / tgo^2) * (rf_star - r - V * tgo);
    
    aT = term1 + term2 + term3 + term4;
end

function dXdt = trajectory(t, X, gamma, kr, tgo0, af_star, g, rf_star, ...
                          V_star, max_thrust, min_thrust, isp)
    % Trajectory dynamics for ODE integration
    
    r = X(1:3);      % Position
    V = X(4:6);      % Velocity
    mass = X(7);     % Mass
    
    g_nondim = g(3);
    tgo = max((tgo0 - t), 0.001);
    
    % Compute thrust acceleration
    aT = compute_aT(gamma, kr, tgo, af_star, g, rf_star, r, V_star, V);
    norm_aT = norm(aT);
    F_mag = norm_aT * mass;
    
    % Mass flow rate (negative because mass decreases)
    dm_dt = -F_mag / (isp * (-g_nondim));
    
    % State derivatives: [velocity; acceleration; mass_rate]
    dXdt = [V; aT; dm_dt];
end

function [c, ceq] = nonlinearConstraints(params, af_star, g, rf_star, r0, ...
                                        vf_star, v0, m0, max_thrust, ...
                                        min_thrust, isp)
    % Nonlinear constraints for optimization
    
    gamma = params(1);
    kr = params(2);
    tgo = params(3);
    
    X0 = [r0; v0; m0];
    tspan = linspace(0, tgo, 300);
    
    odeoptions = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
    [T, state_traj_nd] = ode45(@(t_nd, X_nd) trajectory(t_nd, X_nd, gamma, ...
                                                        kr, tgo, af_star, g, ...
                                                        rf_star, vf_star, ...
                                                        max_thrust, min_thrust, isp), ...
                               tspan, X0, odeoptions);
    
    N = size(state_traj_nd, 1);
    thrustProfile = zeros(N, 1);
    
    % Calculate thrust profile for constraint checking
    for i = 1:N
        radius = state_traj_nd(i, 1:3);
        velocity = state_traj_nd(i, 4:6);
        currentMass = state_traj_nd(i, 7);
        currentTgo = max(tgo - T(i), 0.001);
        aT = compute_aT(gamma, kr, currentTgo, af_star, g, rf_star, ...
                       radius, vf_star, velocity);
        thrustProfile(i) = norm(aT) * currentMass;
    end
    
    % Inequality constraints (c <= 0)
    c1 = -(state_traj_nd(end, 3));                    % Final altitude >= 0
    c2 = state_traj_nd(end, 7) - (18000/62000) * m0; % Final mass >= dry mass
    c3 = max(thrustProfile - max_thrust);             % Thrust <= max_thrust
    c4 = max(min_thrust - thrustProfile);             % Thrust >= min_thrust
    
    c = [c1; c2; c3; c4];
    
    % Equality constraints (none in this problem)
    ceq = [];
end

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