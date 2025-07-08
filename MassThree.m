clear; clc; close all;

%-----------------------
% Initial Guess & Bounds
%-----------------------
x0 = [1; 6; 1]; % Initial Guess for [gamma; kr; tgo]

% Variable bounds for [gamma, kr, tgo]
lb = [0.01; 0.01; 0.01];    % Lower bounds
ub = [5;    15;   3];       % Upper bounds

%-----------------------
% fmincon options
%-----------------------
options = optimoptions('fmincon', ...
    'Display','iter-detailed', ...
    'MaxFunctionEvaluations',1e4, ...
    'FiniteDifferenceType','central');

%-----------------------
% Dimensional values
%-----------------------
g_const = 3.73;                         % Lunar gravity (m/s^2)
r_dim      = [20000; 5000; 10000];      % Initial position (m)
rf_dim     = [0; 0; 0];                  % Target position
v_dim      = [-500; 0; -200];            % Initial velocity (m/s)
vf_dim     = [0; 0; -0.1];               % Target velocity (m/s)
g_dim      = [0; 0; -g_const];          % Grav accel vector
af_dim     = [0; 0; 2*g_const];          % Final thrust accel guess
mass_init  = 62000;                      % Initial mass (kg)
mass_dry   = 18000;                      % Dry mass (kg)
isp_dim    = 330;                        % Isp (s)
max_thrust_dim = 8e5;                   % Max thrust (N)
min_thrust_dim = 0.3*max_thrust_dim;    % Min thrust (N)

%-----------------------
% Reference scales
%-----------------------
L_ref = 1e3;            % Length scale (m)
T_ref = 1e2;            % Time scale (s)
V_ref = L_ref/T_ref;    % Velocity scale (m/s)
A_ref = L_ref/T_ref^2;  % Acceleration scale (m/s^2)
M_ref = mass_init;      % Mass scale (kg)


%-----------------------
% Non-dimensionalize
%-----------------------
r        = r_dim/L_ref;
rf_star = rf_dim/L_ref;
af_star = af_dim/A_ref;
v        = v_dim/V_ref;
vf_star = vf_dim/V_ref;
g        = g_dim/A_ref;

m0      = mass_init/M_ref;
mf_min  = mass_dry/M_ref;
isp     = isp_dim/T_ref;
max_thrust = max_thrust_dim/(M_ref*A_ref);
min_thrust = min_thrust_dim/(M_ref*A_ref);

%-----------------------
% Nonlinear constraints handle
%-----------------------
nonLinearCons = @(x) nonlinearConstraints(x, af_star, g, rf_star, r, vf_star, v, m0, mf_min, max_thrust, min_thrust, isp);

%-----------------------
% Linear inequality on params (example)
%-----------------------
Aineq = [-1  0  0;
          2 -1  0;
          0  0 -1];
bineq = [0; -4; -1];  % ensure tgo >= 1

%-----------------------
% Run optimization
%-----------------------
[x_opt, fval] = fmincon(@(p) objective(p, af_star, g, rf_star, r, vf_star, v), ...
                        x0, Aineq, bineq, [], [], lb, ub, nonLinearCons, options);

gamma_opt = x_opt(1);
kr_opt    = x_opt(2);
tgo_opt   = x_opt(3);

%-----------------------
% Simulate optimal trajectory
%-----------------------
X0 = [r; v; m0];
tspan = [0, tgo_opt];
ode = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',1e-3);
[t_traj, state_traj] = ode45(@(t,X) trajectory(t, X, gamma_opt, kr_opt, tgo_opt, ...
                                             af_star, g, rf_star, vf_star, isp), ...
                                              tspan, X0, node);

%-----------------------
% Plotting
%-----------------------
aT_list = zeros(length(t_traj),3);
mass_list = state_traj(:,7);
for i=1:length(t_traj)
    current_r = state_traj(i,1:3).';
    current_v = state_traj(i,4:6).';
    tgo = max((tgo_opt - t_traj(i)),0.001);
    aT = compute_aT(gamma_opt, kr_opt, tgo_opt, af_star,g,rf_star,current_r,vf_star,current_v).';
    norm_aT = norm(aT);
    F_mag = norm_aT * mass_list(i);
    
    if F_mag > max_thrust
        F_mag = max_thrust;
        aT = (aT/norm_aT)*(max_thrust/mass_list(i));
        norm_aT = norm(aT);
    elseif F_mag < min_thrust
        F_mag = min_thrust;
        aT = (aT/norm_aT)*(min_thrust/mass_list(i));
        norm_aT = norm(aT);
    end
    aT_list(i,:) = aT;
end
aT_norm = vecnorm(aT_list,2,2);

%Check Values at lowest height before landing
% lowestIdx = find(diff(sign(state_traj(:,3))),1,"first");
% lowestHeight = state_traj(lowestIdx,3);
% lowestHeight_Vel = sqrt(state_traj(lowestIdx,4)^2 + state_traj(lowestIdx,5)^2 + state_traj(lowestIdx,6)^2);
% lowestHeightDim = lowestHeight*L_ref
% lowestHeightVelDim = lowestHeight_Vel*V_ref

vf_real_norm_dim = sqrt(state_traj(end,4)^2+state_traj(end,5)^2+state_traj(end,6)^2)*V_ref
zf_real_dim = state_traj(end,3)*L_ref
rf_real_norm_dim = sqrt(state_traj(end,1)^2+state_traj(end,2)^2+state_traj(end,3)^2)*L_ref

figure(1); hold on;
scatter3(state_traj(:,1),state_traj(:,2),state_traj(:,3), 20, aT_norm, 'filled');
xlabel('X');
ylabel('Y');
zlabel('Z');
grid();
view(45,45);
axis('equal');
colormap(brewermap([],'-RdYlBu'));
c=colorbar;
c.Label.String = 'Thrust Acceleration Magnitude';

%plot3(state_traj(:,1),state_traj(:,2),state_traj(:,3));


figure(2); hold on;
plot(t_traj,aT_list(:,1),'r-');
plot(t_traj,aT_list(:,2),'g-');
plot(t_traj,aT_list(:,3),'c-');
plot(t_traj,vecnorm(aT_list,2,2),'w-');
% Below two lines need to not plot min/max thrust, but plot the
% acceleration limits, but that changes with the mass through flight
plot(t_traj, max_thrust./mass_list,'m--');
plot(t_traj, min_thrust./mass_list,'m--');
legend('X Thrust Accel','Y Thrust Accel','Z Thrust Accel','Norm Thrust Accel');
xlabel('Non Dimensional Time');
ylabel('Non Dimensional Acceleration');

figure(3); hold on;
plot(t_traj * T_ref, aT_list(:,1) * A_ref,'r-');
plot(t_traj * T_ref, aT_list(:,2) * A_ref,'g-');
plot(t_traj * T_ref, aT_list(:,3) * A_ref,'c-');
plot(t_traj * T_ref, vecnorm(aT_list,2,2) * A_ref,'w--');
plot(t_traj * T_ref,  max_thrust ./ mass_list*  A_ref,'m--');
plot(t_traj * T_ref,  min_thrust ./ mass_list*  A_ref,'m--');
legend('X Thrust Accel','Y Thrust Accel','Z Thrust Accel','Norm Thrust Accel');
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');

figure(4); hold on;
plot(t_traj,mass_list)
yline(mass_init_dim/M_ref,'r--');
yline(mass_dry_dim/M_ref,'g--');
%%----------------------
%% FUNCTIONS
%%----------------------

function cost = objective(params, af_star, g, rf_star, r0, vf_star, v0)
    gamma = params(1);
    kr    = params(2);
    tgo   = params(3);
    aT = compute_aT(gamma, kr, tgo, af_star, g, rf_star, r0, vf_star, v0);
    cost = norm(aT);
end

function aT = compute_aT(gamma, kr, tgo, af_star, g, rf_star, r, V_star, V)
    term1 = gamma * ((kr / (2*(gamma + 2)) - 1) * af_star);
    term2 = (gamma * kr / (2*(gamma + 2)) - gamma - 1) * g;
    term3 = (gamma + 1)/tgo * (1 - kr / (gamma + 2)) * (V_star - V);
    term4 = kr / tgo^2 * (rf_star - r - V * tgo);
    aT    = term1 + term2 + term3 + term4;
end

% Propagate dynamics WITHOUT clamping thrust inside
function dXdt = trajectory(t, X, gamma, kr, tgo0, af_star, g, rf_star, V_star, isp)
    r    = X(1:3);
    V    = X(4:6);
    mass = X(7);
    g_nd = g(3);
    tau  = max(tgo0 - t, 1e-3);
    aT   = compute_aT(gamma, kr, tau, af_star, g, rf_star, r, V_star, V);
    dm_dt = -norm(aT)*mass/(isp * -g_nd);
    dXdt  = [V; aT; dm_dt];
end

% Enforce final altitude, final mass and thrust path limits as constraints
function [c, ceq] = nonlinearConstraints(params, af_star, g, rf_star, r0, vf_star, v0, m0, mf_min, max_thrust, min_thrust, isp)
    gamma = params(1);
    kr    = params(2);
    tgo0  = params(3);
    X0    = [r0; v0; m0];
    tspan = [0, tgo0];
    opts  = odeset('RelTol',1e-6,'AbsTol',1e-6);

    [t_nd, state_nd] = ode45(@(t,X) trajectory(t, X, gamma, kr, tgo0, af_star, g, rf_star, vf_star, isp), tspan, X0, opts);

    % Final altitude and mass constraints
    r_end = state_nd(end,3);
    m_end = state_nd(end,7);
    c1    = -r_end;             % r_end >= 0  ->  -r_end <= 0
    c2    = mf_min - m_end;     % m_end >= mf_min -> mf_min - m_end <= 0

    % Thrust path constraints at every time node
    N = length(t_nd);
    c_thrust = zeros(2*N,1);
    for i=1:N
        ri     = state_nd(i,1:3)';
        Vi     = state_nd(i,4:6)';
        mi     = state_nd(i,7);
        tau_i  = max(tgo0 - t_nd(i), 1e-3);
        aTi    = compute_aT(gamma, kr, tau_i, af_star, g, rf_star, ri, vf_star, Vi);
        Fmag   = norm(aTi)*mi;
        c_thrust(i)     = Fmag - max_thrust;   % <= 0
        c_thrust(N+i)   = min_thrust - Fmag;   % <= 0
    end

    % Combine all inequalities
    c   = [c1; c2; c_thrust];
    ceq = [];
end