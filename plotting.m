% Plotting function
function plotting(tTraj, stateTraj, optParams, aTList, refVals, problemParams, nonDimParams, flag_thrustGotLimited)
    
    gamma = optParams(1);
    kr = optParams(2);
    tgo0 = optParams(3);
    L_ref = refVals.L_ref;  % m
    T_ref = refVals.T_ref;  % s
    A_ref = refVals.A_ref;  % m/s^2
    V_ref = refVals.V_ref;  % m/s
    M_ref = refVals.M_ref;  % kg
    rMoon = problemParams.rMoon; % m
    E0 = problemParams.E0;
    N0 = problemParams.N0;
    U0 = problemParams.U0;


    rState = stateTraj(:,1:3);
    rDim = rState * L_ref;
    vState = stateTraj(:,4:6);
    vDim = vState * V_ref;
    mState = stateTraj(:,7);
    mDim = mState * M_ref;
    tgoState = tgo0 - tTraj;
    tgoDim = tgoState * T_ref;
    tgo0Dim = tgoDim(1);

    aTDim = aTList * A_ref;
    aTDimNorm = vecnorm(aTDim,2,1);

    %rLanding = rMoon*U0;
    deltaR = MCMF2ENU(rDim',problemParams.landingLatDeg,problemParams.landingLonDeg,true,true);
    East  = deltaR(1,:)';
    North = deltaR(2,:)';
    Up    = deltaR(3,:)';

    alt_m = vecnorm(rDim,2,2) - rMoon;

    vDimENU = MCMF2ENU(vDim',problemParams.landingLatDeg,problemParams.landingLonDeg,false,true);
    vE = vDimENU(1,:)';
    vN = vDimENU(2,:)';
    vU = vDimENU(3,:)';

    aTDimENU = MCMF2ENU(aTDim,problemParams.landingLatDeg,problemParams.landingLonDeg,false,true);
    aE = aTDimENU(1,:)';
    aN = aTDimENU(2,:)';
    aU = aTDimENU(3,:)';
    aT_ENU = [aE, aN, aU];
    aT_norm_ENU = vecnorm(aT_ENU,2,2);

    maxThrustDim = problemParams.maxThrustDim;
    minThrustDim = problemParams.minThrustDim;
    massInitDim = problemParams.massInitDim;
    massDryDim = problemParams.dryMassDim;



% Figure 1: 3D trajectory colored by thrust accel magnitude
    figure(); hold on;
    scatter3(East/1000, North/1000, Up/1000, 20, aT_norm_ENU, 'filled');
    xlabel('East (km)');
    ylabel('North (km)');
    zlabel('Up (km)');
    title('3D Trajectory with Thrust Accel Magnitude');
    grid on; view(90, 0);
    subtitle(sprintf('gamma = %.4f\n kr = %.4f\n tgo = %.4f (s)',gamma, kr, tgo0*T_ref));
    colorbar; colormap(parula);
    axis equal;
    thetaCircle = linspace(0,2*pi,1000);
    x = (rMoon * cos(thetaCircle))/1000;
    y = (rMoon * sin(thetaCircle) - rMoon)/1000;
    plot3(zeros(1000),x,y,'w-');
    xlim([min(East)/1000, max(East)/1000]);
    ylim([min(North)/1000, max(North)/1000]);
    zlim([min(Up)/1000, max(Up)/1000]);

% Figure 2: Range vs Altitude
    rhat = rDim ./ vecnorm(rDim,2,2);
    central = acos(rhat*U0);
    arcLength = rMoon * central;

    figure(); hold on;
    plot(arcLength/1000, alt_m/1000);
    %fprintf("Final X,Y: [%.3f, %.3f]\n",X(end,2),X(end,3))
    xlabel('North km'); ylabel('Up km'); title('Range vs Altitude');
    grid on; yline(0, 'Color', [1 0.2 0.2]);
    ylim([0,50]);
    subtitle(sprintf("East Error: %.2f m\n North Error: %.2f\n Up Error: %.2f", East(end), North(end), Up(end)));

% Figure 3: Thrust acceleration components (dimensional) ENU
    figure(); hold on;
    plot(tTraj*T_ref, aE, 'LineWidth', 1.5);
    plot(tTraj*T_ref, aN, 'LineWidth', 1.5);
    plot(tTraj*T_ref, aU, 'LineWidth', 1.5);
    plot(tTraj*T_ref, aT_norm_ENU, '-', 'LineWidth', 2);
    plot(tTraj(end)*T_ref,norm(nonDimParams.afStarND)*A_ref,'.','MarkerSize',10); % Plot afStar
    plot(tTraj(end)*T_ref,A_ref,'x','MarkerSize',10); % Plot 1g
    plot(tTraj(end)*T_ref,3*A_ref,'x','MarkerSize',10); % Plot 3g
    %xline(tgo0Dim-3, 'r--','LineWidth',1); % Line 3 seconds before, when BTT kicks in ( if enabled and set to 3 seconds)
    legend('East', 'North', 'Up', 'Magnitude','afStar','1g', '3g', 'Location', 'best');
    xlabel('Time s'); ylabel('Accel m/s^2'); title('Thrust Accel Profile (Dim) in ENU frame');
    grid on;

% Figure 4: Thrust acceleration components (dimensional) MCMF
    figure(); hold on;
    plot(tTraj*T_ref, aTDim(1,:), 'LineWidth', 1.5);
    plot(tTraj*T_ref, aTDim(2,:), 'LineWidth', 1.5);
    plot(tTraj*T_ref, aTDim(3,:), 'LineWidth', 1.5);
    plot(tTraj*T_ref, aTDimNorm, '-', 'LineWidth', 2);
    plot(tTraj(end)*T_ref,norm(nonDimParams.afStarND)*A_ref,'.','MarkerSize',10); % Plot afStar
    plot(tTraj(end)*T_ref,A_ref,'x','MarkerSize',10); % Plot 1g
    plot(tTraj(end)*T_ref,3*A_ref,'x','MarkerSize',10); % Plot 3g
    %xline(tgo0Dim-3, 'r--','LineWidth',1); % Line 3 seconds before, when BTT kicks in ( if enabled and set to 3 seconds)
    legend('X', 'Y', 'Z', 'Magnitude','afStar','1g', '3g', 'Location', 'best');
    xlabel('Time s'); ylabel('Accel m/s^2'); title('Thrust Accel Profile (Dim) in MCMF frame');
    if flag_thrustGotLimited
        subtitle("Thrust is being Throttled");
        fprintf("Thrust is being Throttled");
    end
    grid on;

% Figure 5: Velocity components (dimensional)
    figure(); hold on;
    plot(tTraj*T_ref, vE, 'LineWidth', 1.5);
    plot(tTraj*T_ref, vN, 'LineWidth', 1.5);
    plot(tTraj*T_ref, vU, 'LineWidth', 1.5);
    plot(tTraj*T_ref, vecnorm([stateTraj(:,4), stateTraj(:,5), stateTraj(:,6)], 2, 2)*V_ref, '-', 'LineWidth', 2);
    legend('East', 'North', 'Up', 'Magnitude', 'Location', 'best');
    xlabel('Time s'); ylabel('Velocity m/s'); title('Velocity Profile (Dim)');
    grid on;

% Figure 6: Throttle profile (display limits only)
    figure(); hold on;
    % Thrust magnitude = ||aT|| * m, dimensional thrust = a * m * A_ref * M_ref
    thrustDim = aTDimNorm' .* mDim;%/maxThrustDim;
    plot(tTraj*T_ref, thrustDim/maxThrustDim,'DisplayName','Throttle Profile');
    yline(1.0, 'r--', 'LineWidth', 1, 'DisplayName', 'Max Thrust');
    yline(minThrustDim/maxThrustDim, 'r--', 'LineWidth', 1, 'DisplayName', 'Min Thrust');
    xlabel('Time s'); ylabel('Throttle Fraction'); title('Time vs Throttle'); subtitle('Limits only for show, not applied to trajectory')
    legend()

% Figure 7: Mass depletion over ND time
    figure(); hold on;
    plot(tTraj*T_ref, mState*M_ref, 'LineWidth', 2);
    yline(massInitDim, 'r--', 'LineWidth', 1, 'DisplayName', 'Initial');
    yline(massDryDim, 'g--', 'LineWidth', 1, 'DisplayName', 'Dry');
    xlabel('Time (s)'); ylabel('Mass (kg)'); title('Mass Depletion');
    legend('Current', 'Initial', 'Dry', 'Location', 'northeast');
    grid on;
    subtitle(sprintf("Consumed Fuel: %.2f kg", massInitDim - mState(end)*M_ref));

% Figure 8: Time vs Altituide
    figure(); hold on;
    plot(tTraj*T_ref, alt_m/1000, 'LineWidth', 1.5);
    xlabel('Time s'); ylabel('Altitude m'); title('Time vs Altitude (Dimensional)');
    grid on;
end