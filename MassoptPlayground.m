clear; clc; close all;


x0 = [1; 6; 100]; %Initial Guess for gamma, kr, tgo

% Variable bounds: [gamma, kr, tgo]
lb = [0.01, 0.01, 0.01];  % Avoid zero values or negative time to go
ub = [10, 25, 250];

% fmincon options
options = optimoptions('fmincon','Display','iter-detailed','MaxFunctionEvaluations',1000); % Options for fmincon
%'Algorithm','sqp',
%% Dimensional Values
g_const = 9.81;
r_dim = [20000; 2500; 5000]; % Initial position, can be varied
rf_dim = [0; 0; 0]; % Target position
v_dim = [-350; 0; -200]; % Initial velocity
vf_dim = [0; 0; 0]; % Target velocity
g_dim = [0; 0; -g_const]; % grav accel in -z direction
af_dim = [0; 0; 2*g_const];
tgo_min_dim = 1;
mass_init_dim = 62000;
mass_dry_dim = 18000;
isp_dim = 300;
max_thrust_dim = 800000;
min_thrust_dim = 0.3*max_thrust_dim;


%% Reference Values
L_ref = 1000; % Arbitrarily Set, can be redefined if there is a reason to use another set of values
T_ref = 1;
A_ref = L_ref/T_ref^2;
V_ref = L_ref/T_ref;
M_ref = mass_init_dim;

%% Non Dim Values
r = r_dim/L_ref;
rf_star = rf_dim/L_ref;
af_star = af_dim/A_ref;
v = v_dim/V_ref;
vf_star = vf_dim/V_ref;
tgo_min = tgo_min_dim/T_ref;
g = g_dim/A_ref;
m_init = mass_init_dim/M_ref;
m_dry = mass_dry_dim/M_ref;
isp = isp_dim/T_ref;
max_thrust = max_thrust_dim/(M_ref*A_ref);
min_thrust = min_thrust_dim/(M_ref*A_ref);


nonLinearCons = @(x) nonlinearConstraints(x, af_star, g_const, rf_star, r, v, vf_star, m_init, max_thrust, min_thrust, isp); % Pass all needed parameters
%Linear inequalities to handle the constraints on our three parameters
linear_ineq_matrix = [-1 0 0;
                      2 -1 0;
                      0 0 -1];
linear_ineq_vec = [0; -4; -tgo_min];

% Call fmincon to optimize 
[x_opt, fval] = fmincon(@(params) objective(params, af_star, g, rf_star, r, vf_star, v), ...
                          x0, linear_ineq_matrix, linear_ineq_vec, [], [], lb, ub,nonLinearCons, options);

gamma_opt = x_opt(1)
kr_opt = x_opt(2)
tgo_opt = x_opt(3) %This is initial tgo

X0 = [r;v;m_init];
tspan = [0, tgo_opt];
odeoptions = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t_traj, state_traj] = ode45(@(t, X) trajectory(t, X, gamma_opt, kr_opt, tgo_opt,...
                        af_star, g_const, rf_star, vf_star,max_thrust, min_thrust, isp), tspan, X0,odeoptions);

aT_list = zeros(length(t_traj),3);
mass_list = state_traj(:,7)*M_ref;
for i=1:length(t_traj)
    current_r = state_traj(i,1:3).';
    current_v = state_traj(i,4:6).';
    tgo = max((tgo_opt - t_traj(i)),0.01);
    aT_list(i,:) = compute_aT(gamma_opt, kr_opt, tgo_opt, af_star,g,rf_star,current_r,vf_star,current_v).';
end

figure(1); hold on;

plot3(state_traj(:,1),state_traj(:,2),state_traj(:,3));
xlabel('X');
ylabel('Y');
zlabel('Z');
grid();
view(45,45);
axis('equal');

figure(2); hold on;
plot(t_traj,aT_list(:,1),'r-');
plot(t_traj,aT_list(:,2),'g-');
plot(t_traj,aT_list(:,3),'c-');
plot(t_traj,vecnorm(aT_list,2,2),'k--');
legend('X Thrust Accel','Y Thrust Accel','Z Thrust Accel','Norm Thrust Accel');

figure(3); hold on;
plot(t_traj,aT_list(:,1)*A_ref,'r-');
plot(t_traj,aT_list(:,2)*A_ref,'g-');
plot(t_traj,aT_list(:,3)*A_ref,'c-');
plot(t_traj,vecnorm(aT_list,2,2)*A_ref,'k--');
legend('X Thrust Accel','Y Thrust Accel','Z Thrust Accel','Norm Thrust Accel');
%% Functions
%Cost function to optimize
function cost = objective(params, af_star, g, rf_star, r, V_star, V)
    gamma = params(1);
    kr = params(2);
    tgo = params(3);
    
    aT = compute_aT(gamma, kr, tgo, af_star, g, rf_star, r, V_star, V);
    cost = norm(aT);
end

%helper function to calculate current commanded thrust
function aT = compute_aT(gamma, kr, tgo, af_star, g, rf_star, r, V_star, V)
    term1 = gamma * ((kr / (2*(gamma + 2)) - 1) * af_star);
    term2 = (gamma * kr / (2*(gamma + 2)) - gamma - 1) * g;
    term3 = (gamma + 1)/tgo * (1 - kr / (gamma + 2)) * (V_star - V);
    term4 = kr / tgo^2 * (rf_star - r - V * tgo);
    aT = term1 + term2 + term3 + term4;
end

%Trajectory function for ode propagation
function dXdt = trajectory(t, X, gamma, kr, tgo0, af_star, g, rf_star, V_star, max_thrust, min_thrust, isp)
    r= X(1:3);
    V = X(4:6);
    mass = X(7);
    
    tgo = max((tgo0 - t),0.01);
    
    aT = compute_aT(gamma, kr, tgo, af_star, g, rf_star, r, V_star, V);
    norm_aT = norm(aT);
    F_mag = norm_aT * mass;
    
    if F_mag > max_thrust
        F_mag = max_thrust;
        aT = (aT/norm_aT)*(max_thrust/mass);
    elseif F_mag < min_thrust
        F_mag = min_thrust;
        aT = (aT/norm_aT)*(min_thrust/mass);
    end
    dm_dt = -F_mag/ (isp*g);



dXdt = [V;aT;dm_dt];
end

function [c, ceq] = nonlinearConstraints(params, af_star, g, rf_star, r0, v0, vf_star, m0, max_thrust, min_thrust, isp)
    gamma = params(1);
    kr = params(2);
    tgo = params(3);


    X0 = [r0; v0; m0];

    % T_ref is fixed, as per your previous successful change.
    % tgo_opt (initial tgo) is in non-dimensional units
    % The ode45 tspan must be in non-dimensional units, scaled by T_ref.
    tspan = [0, tgo]; % Non-dimensional time span for ODE

    odeoptions = odeset('RelTol',1e-9,'AbsTol',1e-9);
    [~, state_traj_nd] = ode45(@(t_nd, X_nd) trajectory(t_nd, X_nd, gamma, kr, tgo,...
                                    af_star, g, rf_star, vf_star, max_thrust, min_thrust, isp), tspan, X0, odeoptions);


    c = min(state_traj_nd(:,3)); % This will be a vector of constraints, one for each time step

    ceq = [];
end