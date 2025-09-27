function plotLEMMassOpt(S)
%PLOTLEMMASSOPT Produce the seven figures from the simulation results
%   plotLEMMassOpt(S) where S is from runLEMMassOptOpenLoop


%Arranging Variables
    tTraj     = S.tTraj;
    X         = S.stateTraj;
    aT_nd    = S.aTList;
    aT_norm_nd   = S.aT_norm;
    mass_nd  = S.massList;

    L_ref = S.refs.L_ref;  % m
    T_ref = S.refs.T_ref;  % s
    A_ref = S.refs.A_ref;  % m/s^2
    V_ref = S.refs.V_ref;  % m/s
    M_ref = S.refs.M_ref;  % kg
    RMoon = S.refs.R_moon; % m

    Rdim = X(:,1:3) * L_ref; % Trajectory radius in meters MCMF
    Vdim = X(:,4:6) * V_ref; % Trajectory velcoity in meters MCMF
    aTdim = aT_nd * A_ref;   % Commanded thrust in m/s^2
    aT_norm_dim = aT_norm_nd * A_ref; % norm of commnded thrust in m/s^2
    Mdim = mass_nd * M_ref; % Dimensional in kg mass list

    maxThrustDim = S.thrust.maxThrustDim;
    minThrustDim = S.thrust.minThrustDim;

    massInitDim = S.masses.massInitDim;
    massDryDim  = S.masses.massDryDim;


    lat0_deg = S.params.landingLatDeg;
    lon0_deg = S.params.landingLonDeg;
    E0 = S.params.E0;
    N0 = S.params.N0;
    U0 = S.params.U0;
    rLanding = RMoon * U0;


    deltaR = Rdim - rLanding';
    East  = deltaR * E0;
    North = deltaR * N0;
    Up    = deltaR * U0;

    alt_m = vecnorm(Rdim,2,2) - RMoon;

    vE = Vdim * E0;
    vN = Vdim * N0;
    vU = Vdim * U0;

    aE = aTdim * E0;
    aN = aTdim * N0;
    aU = aTdim * U0;
    aT_ENU = [aE, aN, aU];
    aT_norm_ENU = vecnorm(aT_ENU,2,2);

% Figure 1: 3D trajectory colored by thrust accel magnitude
    figure(); hold on;
    scatter3(East/1000, North/1000, Up/1000, 20, aT_norm_ENU, 'filled');
    xlabel('East (km)');
    ylabel('North (km)');
    zlabel('Up (km)');
    title('3D Trajectory with Thrust Accel Magnitude');
    grid on; view(90, 0);
    subtitle(sprintf('gamma = %.4f\n kr = %.4f\n tgo = %.4f (s)',S.opt.gamma, S.opt.kr, S.opt.tgo*T_ref));
    colorbar; colormap(parula);
    axis equal;
    thetaCircle = linspace(0,2*pi,1000);
    x = (RMoon * cos(thetaCircle))/1000;
    y = (RMoon * sin(thetaCircle) - RMoon)/1000;
    plot3(zeros(1000),x,y,'w-');
    xlim([min(East)/1000, max(East)/1000]);
    ylim([min(North)/1000, max(North)/1000]);
    zlim([min(Up)/1000, max(Up)/1000]);

% Figure 2: Range vs Altitude
    rhat = Rdim ./ vecnorm(Rdim,2,2);
    central = acos(rhat*U0);
    arcLength = RMoon * central;

    figure(); hold on;
    plot(arcLength/1000, alt_m/1000);
    %fprintf("Final X,Y: [%.3f, %.3f]\n",X(end,2),X(end,3))
    xlabel('North km'); ylabel('Up km'); title('Range vs Altitude');
    grid on; yline(0, 'Color', [1 0.2 0.2]);
    subtitle(sprintf("East Error: %.2f m\n North Error: %.2f\n Up Error: %.2f", East(end), North(end), Up(end)));
    

    % Figure 3: Thrust acceleration components (dimensional)
    figure(); hold on;
    plot(tTraj*T_ref, aE, 'LineWidth', 1.5);
    plot(tTraj*T_ref, aN, 'LineWidth', 1.5);
    plot(tTraj*T_ref, aU, 'LineWidth', 1.5);
    plot(tTraj*T_ref, vecnorm(aT_nd, 2, 2)*A_ref, '-', 'LineWidth', 2);
    legend('East', 'North', 'Up', 'Magnitude', 'Location', 'best');
    xlabel('Time s'); ylabel('Accel m/s^2'); title('Thrust Accel Profile (Dim)');
    grid on;

    % Figure 4: Velocity components (dimensional)
    figure(); hold on;
    plot(tTraj*T_ref, vE, 'LineWidth', 1.5);
    plot(tTraj*T_ref, vN, 'LineWidth', 1.5);
    plot(tTraj*T_ref, vU, 'LineWidth', 1.5);
    plot(tTraj*T_ref, vecnorm([X(:,4), X(:,5), X(:,6)], 2, 2)*V_ref, '-', 'LineWidth', 2);
    legend('East', 'North', 'Up', 'Magnitude', 'Location', 'best');
    xlabel('Time s'); ylabel('Velocity m/s'); title('Velocity Profile (Dim)');
    grid on;

    % Figure 5: Throttle profile (display limits only)
    figure(); hold on;
    % Thrust magnitude = ||aT|| * m, dimensional thrust = a * m * A_ref * M_ref
    thrustDim = vecnorm(aT_nd .* mass_nd, 2, 2) * A_ref * M_ref;
    plot(tTraj*T_ref, thrustDim / maxThrustDim,'DisplayName','Throttle Profile');
    yline(1.0, 'r--', 'LineWidth', 1, 'DisplayName', 'Max Thrust');
    yline(minThrustDim/maxThrustDim, 'r--', 'LineWidth', 1, 'DisplayName', 'Min Thrust');
    xlabel('Time s'); ylabel('Throttle Fraction'); title('Time vs Throttle'); subtitle('Limits only for show, not applied to trajectory')
    legend()

    % Figure 6: Mass depletion over ND time
    figure(); hold on;
    plot(tTraj*T_ref, mass_nd*M_ref, 'LineWidth', 2);
    yline(massInitDim, 'r--', 'LineWidth', 1, 'DisplayName', 'Initial');
    yline(massDryDim, 'g--', 'LineWidth', 1, 'DisplayName', 'Dry');
    xlabel('Time (s)'); ylabel('Mass (kg)'); title('Mass Depletion');
    legend('Current', 'Initial', 'Dry', 'Location', 'northeast');
    grid on;
    subtitle(sprintf("Consumed Fuel: %.2f kg", massInitDim - mass_nd(end)*M_ref));

    % Figure 7: Time vs Altituide
    figure(); hold on;
    plot(tTraj*T_ref, alt_m/1000, 'LineWidth', 1.5);
    xlabel('Time s'); ylabel('Altitude m'); title('Time vs Altituide (Dimensional)');
    grid on;

    

end
function [E, N, U] = enuBasis(lat, lon)
    clat = cos(lat); slat = sin(lat);
    clon = cos(lon); slon = sin(lon);
    U = [clat*clon; clat*slon; slat];
    E = [-slon; clon; 0];
    N = [-slat*clon; -slat*slon; clat];
end