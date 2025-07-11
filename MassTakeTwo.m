clear; clc; close all;


x0 = [1; 6; 1]; %Initial Guess for gamma, kr, tgo

% Variable bounds: [gamma, kr, tgo]
lb = [0.1, 4.2, 0.5];  % Avoid zero values or negative time to go
ub = [5, 25, 3];

% fmincon options
options = optimoptions('fmincon','Display','iter-detailed','MaxFunctionEvaluations',1000);%,'EnableFeasibilityMode',true,'SubproblemAlgorithm','cg'); % Options for fmincon
%'Algorithm','sqp',
%% Dimensional Values
g_const = 3.73;
r_dim = [20000; 3000; 5000]; % Initial position, can be varied
rf_dim = [0; 0; 0]; % Target position
v_dim = [-500; 0; -200]; % Initial velocity
vf_dim = [0; 0; -0.1]; % Target velocity
g_dim = [0; 0; -g_const]; % grav accel in -z direction
af_dim = [0; 0; 2*g_const];
tgo_min_dim = 1;

mass_init_dim = 62000;
mass_dry_dim = 18000;
isp_dim = 330;
max_thrust_dim = 800000;
min_thrust_dim = 0.3*max_thrust_dim;

%% Reference Values
L_ref = 1000; % Arbitrarily Set, can be redefined if there is a reason to use another set of values
T_ref = 100;
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
m0 = mass_init_dim/M_ref;
mf_min = mass_dry_dim/M_ref;
isp = isp_dim/T_ref;
max_thrust = max_thrust_dim/(M_ref*A_ref);
min_thrust = min_thrust_dim/(M_ref*A_ref);

nonLinearCons = @(x) nonlinearConstraints(x, af_star, g, rf_star, r, v, vf_star, m0, max_thrust, min_thrust, isp); % Pass all needed parameters
%Linear inequalities to handle the constraints on our three parameters
linear_ineq_matrix = [-1 0 0; % gamma >= 0
                      2 -1 0; % 2gamma - kr <= -4  --> 2(gamma+2) < = kr
                      0 0 -1]; % - tgo <= -tgo_min --> tgo >= tgo_min
linear_ineq_vec = [0; -4; -tgo_min];

% Call fmincon to optimize 
[x_opt, fval] = fmincon(@(params) objective(params, af_star, g, rf_star, r, vf_star, v, m0, max_thrust, min_thrust, isp), ...
                          x0, linear_ineq_matrix, linear_ineq_vec, [], [], lb, ub,nonLinearCons, options);

gamma_opt = x_opt(1)
kr_opt = x_opt(2)
tgo_opt = x_opt(3) %This is initial tgo
tgo_opt_dim = tgo_opt*T_ref

X0 = [r;v;m0];
tspan = [0, tgo_opt];
odeoptions = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',0.001);
[t_traj, state_traj] = ode45(@(t, X) trajectory(t, X, gamma_opt, kr_opt, tgo_opt,...
                        af_star, g, rf_star, vf_star, max_thrust, min_thrust, isp), tspan, X0,odeoptions);

aT_list = zeros(length(t_traj),3);
mass_list = state_traj(:,7);
for i=1:length(t_traj)
    current_r = state_traj(i,1:3).';
    current_v = state_traj(i,4:6).';
    tgo = max((tgo_opt - t_traj(i)),0.001);
    aT = compute_aT(gamma_opt, kr_opt, tgo, af_star,g,rf_star,current_r,vf_star,current_v).';
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
%% Functions
%Cost function to optimize
function cost = objective(params, af_star, g, rf_star, r, vf_star, v, m, max_thrust, min_thrust, isp)
    gamma = params(1);
    kr = params(2);
    tgo = params(3);
    N = 300;
    tspan = linspace(0,tgo,N);
    odeoptions = odeset('RelTol',1e-6,'AbsTol',1e-6);
    X0 = [r; v; m];

    [T, X] = ode45(@(t_nd, X_nd) trajectory(t_nd, X_nd, gamma, kr, tgo,...
                                    af_star, g, rf_star, vf_star, max_thrust, min_thrust, isp), tspan, X0, odeoptions);
    thrustProfile = zeros(1,N);
    for i=1:N
        r_i = X(i,1:3);
        v_i = X(i,4:6);
        tau = max(tgo - T(i),0.001);
      aT_i   = compute_aT(gamma, kr, tau, af_star, g, rf_star, r_i, vf_star, v_i);
      thrustProfile(i) = norm(aT_i)*X(i,7);  % aT*m → force
    end
    if mod(length(T)-1,2)==1
        T = T(1:end-1);
        thrustProfile = thrustProfile(1:end-1);
    end
    cost = simpsonIntegral(T, thrustProfile);

    %aT = compute_aT(gamma, kr, tgo, af_star, g, rf_star, r, V_star, V);
    %cost = norm(aT);
end

%helper function to calculate current commanded thrust
function aT = compute_aT(gamma, kr, tgo, af_star, g, rf_star, r, V_star, V)
    term1 = gamma * (((kr / (2*(gamma + 2))) - 1) * af_star);
    term2 = ((gamma * kr / (2*(gamma + 2))) - gamma - 1) * g;
    term3 = (gamma + 1)/tgo * (1 - (kr / (gamma + 2))) * (V_star - V);
    term4 = (kr / tgo^2) * (rf_star - r - V * tgo);
    aT = term1 + term2 + term3 + term4;
end

%Trajectory function for ode propagation
function dXdt = trajectory(t, X, gamma, kr, tgo0, af_star, g, rf_star, V_star, max_thrust, min_thrust, isp)
    r= X(1:3);
    V = X(4:6);
    mass = X(7);
    
    g_nondim = g(3);
    tgo = max((tgo0 - t),0.001);
    
    aT = compute_aT(gamma, kr, tgo, af_star, g, rf_star, r, V_star, V);
    norm_aT = norm(aT);
    F_mag = norm_aT * mass;
    
    % if F_mag > max_thrust
    %     F_mag = max_thrust;
    %     aT = (aT/norm_aT)*(max_thrust/mass);
    %     norm_aT = norm(aT);
    % elseif F_mag < min_thrust
    %     F_mag = min_thrust;
    %     aT = (aT/norm_aT)*(min_thrust/mass);
    %     norm_aT = norm(aT);
    % end
    dm_dt = -F_mag/ (isp*-g_nondim);

    dXdt = [V;aT;dm_dt];
end

% Nonlinear Constraint Equations, keeps Z positive, and Mass within range
function [c, ceq] = nonlinearConstraints(params, af_star, g, rf_star, r0, vf_star, v0, m0, max_thrust, min_thrust, isp)
    gamma = params(1);
    kr = params(2);
    tgo = params(3);


    X0 = [r0; v0; m0];

    tspan = linspace(0,tgo,300); % Non-dimensional time span for ODE

    odeoptions = odeset('RelTol',1e-6,'AbsTol',1e-6);
    [T, state_traj_nd] = ode45(@(t_nd, X_nd) trajectory(t_nd, X_nd, gamma, kr, tgo,...
                                    af_star, g, rf_star, vf_star, max_thrust, min_thrust, isp), tspan, X0, odeoptions);

    N = size(state_traj_nd,1);
   thrustProfile = zeros(N,1);
   for i = 1:N
       radius = state_traj_nd(i,1:3);
       velocity = state_traj_nd(i,4:6);
       currentMass = state_traj_nd(i,7);
       currentTgo = max(tgo - T(i),0.001);
       aT = compute_aT(gamma, kr, currentTgo, af_star, g, rf_star, radius, vf_star, velocity);
       thrustProfile(i) = norm(aT)*currentMass;
   end
    
    
    % Ineq > 0
    %c1 = -min(state_traj_nd(:,3)); % This will be a vector of constraints, one for each time step


    c1 = -(state_traj_nd(end,3)); % Z constraint
    c2 = state_traj_nd(end, 7) - 18000/62000 * m0; % Final mass constraint
    c3 = max(thrustProfile - max_thrust);
    c4 = max(min_thrust - thrustProfile);
    c = [c1 c2 c3 c4];
    
    %ceq1 = min(state_traj_nd(:,3));
    ceq = [];
end
% Simspons Rule Integral: Composite 1/3
function I = simpsonIntegral(t, y)
    N = length(t)-1;
    if mod(N,2) ~=0
        error("Simpson's needs an even number of intervals, odd number of points.");
    end

    h = (t(end) - t(1))/N;
    I = y(1) + y(end) + 4*sum(y(2:2:end-1)) + 2*sum(y(3:2:end-2));
    I = I * (h/3);
end