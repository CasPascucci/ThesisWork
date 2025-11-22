% Plotting function for Single Runs
function plotting(tTraj, stateTraj, optParams, optCost, aTOptim, mOptim, rdOptim, vdOptim, aTList, refVals, problemParams, nonDimParams, optimParams, flag_thrustGotLimited, unconstrained)
    % Check if simulation data is available
    hasSimData = ~isempty(tTraj) && ~isempty(stateTraj) && ~isempty(aTList);
    
    gamma = optParams(1);
    kr = optParams(2);
    tgo0 = optParams(3);
    L_ref = refVals.L_ref;  % m
    T_ref = refVals.T_ref;  % s
    A_ref = refVals.A_ref;  % m/s^2
    V_ref = refVals.V_ref;  % m/s
    M_ref = refVals.M_ref;  % kg
    rMoon = problemParams.rMoon; % m
    [E0, N0, U0] = enuBasis(deg2rad(problemParams.landingLatDeg),deg2rad(problemParams.landingLonDeg));
    isp = nonDimParams.ispND;

    % Process simulation data if available
    if hasSimData
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
    end
    
    
    hasUC = exist('unconstrained','var') && ~isempty(unconstrained);
    if hasUC
        tTrajUC = unconstrained.tTraj;
        rStateUC = unconstrained.stateTraj(:,1:3);
        vStateUC = unconstrained.stateTraj(:,4:6);
        mStateUC = unconstrained.stateTraj(:,7);
        rDimUC = rStateUC * L_ref;
        vDimUC = vStateUC * V_ref;
        mDimUC = mStateUC * M_ref;
        deltaR_UC   = MCMF2ENU(rDimUC', problemParams.landingLatDeg, problemParams.landingLonDeg, true,  true);
        vDimENU_UC  = MCMF2ENU(vDimUC', problemParams.landingLatDeg, problemParams.landingLonDeg, false, true);
        EastUC = deltaR_UC(1,:)';
        NorthUC = deltaR_UC(2,:)';
        UpUC = deltaR_UC(3,:)';

        alt_m_UC = vecnorm(rDimUC,2,2) - rMoon;
        aTDimUC = unconstrained.aTList * A_ref;
        aTDimENU_UC = MCMF2ENU(aTDimUC, problemParams.landingLatDeg, problemParams.landingLonDeg, false, true);
        aE_UC = aTDimENU_UC(1,:)';
        aN_UC = aTDimENU_UC(2,:)';
        aU_UC = aTDimENU_UC(3,:)';
        aT_norm_ENU_UC = vecnorm([aE_UC, aN_UC, aU_UC], 2, 2);

        aTDimNormUC = vecnorm(aTDimUC,2,1)';
    end



    maxThrustDim = problemParams.maxThrustDim;
    minThrustDim = problemParams.minThrustDim;
    massInitDim = problemParams.massInitDim;
    massDryDim = problemParams.dryMassDim;
%% Optim Figures
    nodeCount = optimParams.nodeCount;
    aTOptim = aTOptim';
    rdOptim = rdOptim';
    vdOptim = vdOptim';
    mOptim = mOptim';
    rdOptimUC = unconstrained.rdOptim';

    rdOptimTOPO = MCMF2ENU(rdOptim',problemParams.landingLatDeg,problemParams.landingLonDeg,true,false);
    rdOptimTOPODim = rdOptimTOPO*L_ref;

    tgospanOpt = linspace(0,tgo0,nodeCount);
    tspanOpt = tgo0 - tgospanOpt;
    aTNormOpt = vecnorm(aTOptim,2,2);

    EOpt = rdOptimTOPODim(1, :);
    NOpt = rdOptimTOPODim(2, :);
    UOpt = rdOptimTOPODim(3, :);
% Optim Throttle
    figure('Name',"Optim Throttle"); hold on;
    thrustDim = aTNormOpt .* mOptim *(refVals.M_ref*refVals.A_ref);
    plot(tspanOpt(:)*T_ref, thrustDim(:)/problemParams.maxThrustDim,'DisplayName','Throttle Profile');
    if hasUC && isfield(unconstrained,'aTOptim') && ~isempty(unconstrained.aTOptim) ...
         && isfield(unconstrained,'mOptim')  && ~isempty(unconstrained.mOptim)

        aTOptimUC = unconstrained.aTOptim';
        mOptimUC  = unconstrained.mOptim';
        aTNormOptUC = vecnorm(aTOptimUC, 2, 2);
        thrustDimUC_optim = aTNormOptUC .* mOptimUC * (refVals.M_ref*refVals.A_ref);
        tgo0UC = unconstrained.optParams(3);

        nodeCountUC = size(aTOptimUC,1);
        tgospanOptUC = linspace(0, tgo0UC, nodeCountUC);
        tspanOptUC   = tgo0UC - tgospanOptUC;
    
        plot(tspanOptUC(:)*T_ref, thrustDimUC_optim(:)/problemParams.maxThrustDim, ...
            'r-', 'LineWidth', 1.3, 'DisplayName','Throttle Profile (Unconstrained)');
        
    end
    yline(1.0, 'r--', 'LineWidth', 1, 'DisplayName', 'Max Thrust');
    if optimParams.pointingEnabled
        yline(0.95, 'k:', 'LineWidth', 1, 'DisplayName', '95% Max (If Pointing)');
    end
    yline(problemParams.minThrustDim/problemParams.maxThrustDim, 'r--', 'LineWidth', 1, 'DisplayName', 'Min Thrust');
    xlabel('Time s'); ylabel('Throttle Fraction'); title('Time vs Throttle');
    if hasUC
        subtitle(sprintf("Optim - gamma: %.2f, kr: %.2f, tgo: %.2f s\n UC - gamma: %.2f, kr: %.2f, tgo: %.2f s", optParams(1),optParams(2),optParams(3)*T_ref,unconstrained.optParams(1),unconstrained.optParams(2),unconstrained.optParams(3)*T_ref));
    end
    optFuelCost = (mOptim(end)-mOptim(1))*M_ref;
    subtitle(sprintf("Consumed Fuel: %.2f kg", optFuelCost));
    legend('Location','best');

% Optim Range vs Altitude
    rdOptimDim = rdOptim*L_ref;
    rhatOpt = rdOptimDim ./ vecnorm(rdOptimDim,2,2);
    centralOpt = acos(rhatOpt*U0);
    arcLengthOpt = rMoon * centralOpt;
    alt_opt = vecnorm(rdOptimDim,2,2) - rMoon;

    figure('Name','Opt Range vs Altitude'); hold on;

    if ~isempty(rdOptimUC)
        rdOptimUCDim = rdOptimUC * L_ref;
        rhatOptUC = rdOptimUCDim ./ vecnorm(rdOptimUCDim,2,2);
        centralOptUC = acos(rhatOptUC*U0);
        arcLengthOptUC = rMoon * centralOptUC;
        alt_optUC = vecnorm(rdOptimUCDim,2,2) - rMoon;
        plot(arcLengthOptUC/1000, alt_optUC/1000, 'r-', 'DisplayName',"Unconstrained Trajectory");
    end

    plot(arcLengthOpt/1000, alt_opt/1000, 'b-', 'DisplayName',"Constrained Trajectory");
    %fprintf("Final X,Y: [%.3f, %.3f]\n",X(end,2),X(end,3))
    xlabel('North km'); ylabel('Up km'); title('Range vs Altitude');
    grid on;
    altMax = 1.5 * max(max(alt_opt),max(alt_optUC)) / 1000;
    ylim([0,altMax]);
    subtitle(sprintf("East Error: %.2f m\n North Error: %.2f\n Up Error: %.2f", EOpt(1), NOpt(1), UOpt(1)));
    legend();

% Optim 3D Plot
    Ntotal = size(rdOptimTOPODim, 2);% total node count
    idx_hi = floor(0.15 * Ntotal); % match the high index from constraints
    idx_lo = floor(0.02 * Ntotal) + 1; % match the low index from constraints
    idx = idx_lo:idx_hi; % range
    
    EOpt = EOpt(1:idx_hi);
    NOpt = NOpt(1:idx_hi);
    UOpt = UOpt(1:idx_hi);
    
    ETraj = rdOptimTOPODim(1, 1:idx_hi);
    NTraj = rdOptimTOPODim(2, 1:idx_hi);
    UTraj = rdOptimTOPODim(3, 1:idx_hi);
    
    theta_fun = @(u) min(89.9, 45 + 45 * max(0, min(1, (u - 250) ./ 250))); % enforces theta piecewise nature
    Umin = min(UOpt); % lowest height value in the idx range
    Umax = min(500, max(UOpt, [], 'omitnan')); % highest height value in idx, capped at 500
    U_samp = linspace(Umin, Umax, 360);
    theta_deg = theta_fun(U_samp);
    
    cosBound = cosd(theta_deg);
    cosBound = max(cosBound, 1e-3);
    R_samp = U_samp .* sqrt(1 - cosBound.^2) ./ cosBound; % R=U*tan(theta), and then sin = sqrt(1-cos^2)
    Rcap = 3000;
    R_samp = min(R_samp, Rcap);
    
    azimuth = linspace(0, 2*pi, 360);
    [azGrid, altGrid] = meshgrid(azimuth, U_samp);
    [~, rGrid] = meshgrid(azimuth, R_samp);
    eastGrid  = rGrid .* cos(azGrid);
    northGrid = rGrid .* sin(azGrid);
    
    % Plateau ring at top
    UPlat = Umax;
    RPlat = Rcap;
    azPlat = linspace(0, 2*pi, 360);
    eastPlat = RPlat * cos(azPlat);
    northPlat = RPlat * sin(azPlat);
    UPlat = UPlat * ones(size(azPlat));
    
    figure('Name','Opt 3D'); hold on; grid on; axis equal;
    plot3(ETraj/1000, NTraj/1000, UTraj/1000, 'b-', 'LineWidth', 2, 'MarkerSize',15);
    surf(eastGrid/1000, northGrid/1000, altGrid/1000,altGrid/1000,'EdgeAlpha',0.15,'FaceAlpha',0.4,'MeshStyle','row','LineWidth',0.8);
    plot3(eastPlat/1000, northPlat/1000, UPlat/1000, 'k--', 'LineWidth', 0.8);
    if ~isempty(rdOptimUC)
        rdOptimUC_TOPO = MCMF2ENU(rdOptimUC',problemParams.landingLatDeg,problemParams.landingLonDeg,true,false);
        rdOptimUC_TOPO_DIM = rdOptimUC_TOPO * L_ref;
        ETrajUC = rdOptimUC_TOPO_DIM(1, 1:idx_hi);
        NTrajUC = rdOptimUC_TOPO_DIM(2, 1:idx_hi);
        UTrajUC = rdOptimUC_TOPO_DIM(3, 1:idx_hi);
        plot3(ETrajUC/1000, NTrajUC/1000, UTrajUC/1000, 'r-', 'LineWidth', 2);
    end
    xlabel('East (km)');
    ylabel('North (km)');
    zlabel('Up (km)');
    title('3D Trajectory with Cosine-Based Glide-Slope Envelope (2%–20% nodes)');
    if hasUC
        subtitle(sprintf("Constrained - gamma: %.2f, kr: %.2f, tgo: %.2f \n Thrust Only - gamma: %.2f, kr: %.2f, tgo: %.2f s", optParams(1),optParams(2),optParams(3)*T_ref,unconstrained.optParams(1),unconstrained.optParams(2),unconstrained.optParams(3)*T_ref));
    end
    view(80, 15); camproj orthographic;
    zlim([0,2]);
    xlim([-3,3]);
    ylim([-3,3]);
    axis square;
    legend("Constrained Trajectory","","","Thrust Constraint Only Trajectory","Location","bestoutside");



% Optim Velocity Components 
    vdOptimDim = vdOptim * V_ref;
    vdOptimTOPO = MCMF2ENU(vdOptimDim',problemParams.landingLatDeg,problemParams.landingLonDeg,false,true);
    figure('Name','Opt Vel TOPO'); hold on;
    plot(tspanOpt*T_ref, vdOptimTOPO(1,:), 'LineWidth', 1.5);
    plot(tspanOpt*T_ref, vdOptimTOPO(2,:), 'LineWidth', 1.5);
    plot(tspanOpt*T_ref, vdOptimTOPO(3,:), 'LineWidth', 1.5);
    plot(tspanOpt*T_ref, vecnorm(vdOptimTOPO, 2, 1), '-', 'LineWidth', 2);
    legend('East', 'North', 'Up', 'Magnitude', 'Location', 'best');
    xlabel('Time s'); ylabel('Velocity m/s'); title('Velocity Profile (Dim)');
    grid on;
    xlim([0, ceil(max(tspanOpt*T_ref)/100)*100]);

% Optim Acceleration Components Topo
    aTOptimDim = aTOptim * A_ref;
    aTOptimTOPO = MCMF2ENU(aTOptimDim',problemParams.landingLatDeg,problemParams.landingLonDeg,false,true);
    figure('Name','Opt TOPO Accel'); hold on;
    plot(tspanOpt(:)*T_ref, aTOptimTOPO(1,:), 'LineWidth', 1.5);
    plot(tspanOpt(:)*T_ref, aTOptimTOPO(2,:), 'LineWidth', 1.5);
    plot(tspanOpt(:)*T_ref, aTOptimTOPO(3,:), 'LineWidth', 1.5);
    plot(tspanOpt(:)*T_ref, vecnorm(aTOptimTOPO(:,:),2,1), '-', 'LineWidth', 2);
    plot(tspanOpt(1)*T_ref,norm(nonDimParams.afStarND)*A_ref,'.','MarkerSize',10); % Plot afStar
    plot(tspanOpt(1)*T_ref,A_ref,'x','MarkerSize',10); % Plot 1g
    plot(tspanOpt(1)*T_ref,3*A_ref,'x','MarkerSize',10); % Plot 3g
    legend('East', 'North', 'Up', 'Magnitude','afStar','1g', '3g', 'Location', 'best');
    xlabel('Time s'); ylabel('Accel m/s^2'); title('Commanded Accel Profile (Dim) in TOPO frame');
    grid on;

% Optim Acceleration Components MCMF
    aTOptimDim = aTOptimDim';
    figure('Name','Opt MCMF Accel'); hold on;
    plot(tspanOpt(:)*T_ref, aTOptimDim(1,:), 'LineWidth', 1.5);
    plot(tspanOpt(:)*T_ref, aTOptimDim(2,:), 'LineWidth', 1.5);
    plot(tspanOpt(:)*T_ref, aTOptimDim(3,:), 'LineWidth', 1.5);
    plot(tspanOpt(:)*T_ref, vecnorm(aTOptimDim(:,:),2,1), '-', 'LineWidth', 2);
    plot(tspanOpt(1)*T_ref,norm(nonDimParams.afStarND)*A_ref,'.','MarkerSize',10); % Plot afStar
    plot(tspanOpt(1)*T_ref,A_ref,'x','MarkerSize',10); % Plot 1g
    plot(tspanOpt(1)*T_ref,3*A_ref,'x','MarkerSize',10); % Plot 3g
    %xline(tgo0Dim-3, 'r--','LineWidth',1); % Line 3 seconds before, when BTT kicks in ( if enabled and set to 3 seconds)
    legend('X', 'Y', 'Z', 'Magnitude','afStar','1g', '3g', 'Location', 'best');
    xlabel('Time s'); ylabel('Accel m/s^2'); title('Commanded Accel Profile (Dim) in MCMF frame');
    grid on;
    
%% Sim Figures
if hasSimData
    % Only generate simulation figures if simulation data is available
    % Figure 1: 3D trajectory with glideslope Cone
    idx2km = find(alt_m <= 2000, 1, 'first');
    if isempty(idx2km)
        % Fallback: if never reaches 2 km, just use the last 20%
        NtotalSim = numel(alt_m);
        idx2km = max(1, NtotalSim - floor(0.20*NtotalSim) + 1);
    end
    idx_sim = idx2km:numel(alt_m);

    ETrajSim = East(idx_sim);
    NTrajSim = North(idx_sim);
    UTrajSim = Up(idx_sim);

    if hasUC
        idx2kmUC = find(alt_m_UC <= 2000, 1, 'first');
        if isempty(idx2kmUC)
            NtotalSimUC = numel(alt_m_UC);
            idx2kmUC = max(1, NtotalSimUC - floor(0.20*NtotalSimUC) + 1);
        end
        idx_simUC = idx2kmUC:numel(alt_m_UC);

        ETrajSimUC = EastUC(idx_simUC);
        NTrajSimUC = NorthUC(idx_simUC);
        UTrajSimUC = UpUC(idx_simUC);
    end

    theta_fun = @(u) min(89.9, 45 + 45*max(0, min(1, (u-250)./250)));
    UminEnv = min(UTrajSim);                       % use sim heights
    UmaxEnv = min(500, max(UTrajSim,[],'omitnan'));
    U_samp  = linspace(UminEnv, UmaxEnv, 360);
    theta_deg = theta_fun(U_samp);
    cosBound  = max(cosd(theta_deg), 1e-3);
    R_samp    = U_samp .* sqrt(1 - cosBound.^2) ./ cosBound;
    Rcap      = 3000;
    R_samp    = min(R_samp, Rcap);
    azimuth         = linspace(0, 2*pi, 360);
    [azGrid, altGrid] = meshgrid(azimuth, U_samp);
    [~, rGrid]        = meshgrid(azimuth, R_samp);
    eastGrid  = rGrid .* cos(azGrid);
    northGrid = rGrid .* sin(azGrid);

    % Plateau ring at top
    UPlat = UmaxEnv; RPlat = Rcap;
    azPlat = linspace(0, 2*pi, 360);
    eastPlat  = RPlat * cos(azPlat);
    northPlat = RPlat * sin(azPlat);
    UPlat     = UPlat * ones(size(azPlat));

    figure('Name','Sim 3D (Constrained vs Unconstrained)'); hold on; grid on; axis equal;
    surf(eastGrid/1000, northGrid/1000, altGrid/1000,altGrid/1000,'EdgeAlpha',0.15,'FaceAlpha',0.4,'MeshStyle','row','LineWidth',0.8);
    plot3(eastPlat/1000, northPlat/1000, UPlat/1000, 'k--', 'LineWidth', 0.8);
    plot3(ETrajSim/1000, NTrajSim/1000, UTrajSim/1000, 'b-', 'LineWidth', 2);
    if hasUC
        plot3(ETrajSimUC/1000, NTrajSimUC/1000, UTrajSimUC/1000, 'r-', 'LineWidth', 2);
    end

    xlabel('East (km)'); ylabel('North (km)'); zlabel('Up (km)');
    title('3D Trajectory (Sim) with Cosine-Based Glide-Slope Envelope, Final 2 KM');
    subtitle(sprintf("Constrained - gamma: %.2f, kr: %.2f, tgo: %.2f \n Thrust Only - gamma: %.2f, kr: %.2f, tgo: %.2f s", optParams(1),optParams(2),optParams(3)*T_ref,unconstrained.optParams(1),unconstrained.optParams(2),unconstrained.optParams(3)*T_ref));
    view(80, 15); camproj orthographic;
    zlim([0, 2]); xlim([-3,3]); ylim([-3,3]); axis square;
    legend("","","Constrained Trajectory","Thrust Constraint Only Trajectory","Location","bestoutside");

% Figure 2: Range vs Altitude
    rhat      = rDim ./ vecnorm(rDim,2,2);
    central   = acos(rhat*U0);
    arcLength = rMoon * central;   % meters
    
    figure('Name','Sim Range vs Altitude'); hold on;
    plot(arcLength/1000, alt_m/1000, 'b-', 'DisplayName','Constrained (Sim)');
    
    if hasUC
        rhatUC      = rDimUC ./ vecnorm(rDimUC,2,2);
        centralUC   = acos(rhatUC*U0);
        arcLengthUC = rMoon * centralUC;
        plot(arcLengthUC/1000, alt_m_UC/1000, 'r-', 'DisplayName','Unconstrained (Sim)');
    end

    xlabel('Range (km)'); ylabel('Up (km)');
    title('Range vs Altitude (Sim)'); grid on;
    legend('Location','best');
    if hasUC
        subtitle(sprintf(['Final positions (East, North, Up)\n' ...
            'Constrained: (%.3f, %.3f, %.3f) m\n' ...
            'Unconstrained: (%.3f, %.3f, %.3f) m'], ...
            East(end), North(end), Up(end), ...
            EastUC(end), NorthUC(end), UpUC(end)));
    else
        subtitle(sprintf('Final position (East, North, Up): (%.1f, %.1f, %.1f) m', ...
            East(end), North(end), Up(end)));
    end

% Figure 3: Thrust acceleration components (dimensional) TOPO
    figure('Name','Sim TOPO Accel'); hold on;
    plot(tTraj*T_ref, aE, 'LineWidth', 1.5);
    plot(tTraj*T_ref, aN, 'LineWidth', 1.5);
    plot(tTraj*T_ref, aU, 'LineWidth', 1.5);
    plot(tTraj*T_ref, aT_norm_ENU, '-', 'LineWidth', 2);
    plot(tTraj(end)*T_ref,norm(nonDimParams.afStarND)*A_ref,'.','MarkerSize',10); % Plot afStar
    plot(tTraj(end)*T_ref,A_ref,'x','MarkerSize',10); % Plot 1g
    plot(tTraj(end)*T_ref,3*A_ref,'x','MarkerSize',10); % Plot 3g
    %xline(tgo0Dim-3, 'r--','LineWidth',1); % Line 3 seconds before, when BTT kicks in ( if enabled and set to 3 seconds)
    legend('East', 'North', 'Up', 'Magnitude','afStar','1g', '3g', 'Location', 'best');
    xlabel('Time s'); ylabel('Accel m/s^2'); title('Commanded Accel Profile (Dim) in TOPO frame');
    grid on;

% Figure 4: Thrust acceleration components (dimensional) MCMF
    figure('Name','Sim MCMF Accel'); hold on;
    plot(tTraj*T_ref, aTDim(1,:), 'LineWidth', 1.5);
    plot(tTraj*T_ref, aTDim(2,:), 'LineWidth', 1.5);
    plot(tTraj*T_ref, aTDim(3,:), 'LineWidth', 1.5);
    plot(tTraj*T_ref, aTDimNorm, '-', 'LineWidth', 2);
    plot(tTraj(end)*T_ref,norm(nonDimParams.afStarND)*A_ref,'.','MarkerSize',10); % Plot afStar
    plot(tTraj(end)*T_ref,A_ref,'x','MarkerSize',10); % Plot 1g
    plot(tTraj(end)*T_ref,3*A_ref,'x','MarkerSize',10); % Plot 3g
    %xline(tgo0Dim-3, 'r--','LineWidth',1); % Line 3 seconds before, when BTT kicks in ( if enabled and set to 3 seconds)
    legend('X', 'Y', 'Z', 'Magnitude','afStar','1g', '3g', 'Location', 'best');
    xlabel('Time s'); ylabel('Accel m/s^2'); title('Commanded Accel Profile (Dim) in MCMF frame');
    if flag_thrustGotLimited
        subtitle("Thrust is being Throttled");
        fprintf("Thrust is being Throttled");
    end
    grid on;

% Figure 5: Velocity components (dimensional)
    figure('Name','Sim Vel'); hold on;
    plot(tTraj*T_ref, vE, 'LineWidth', 1.5);
    plot(tTraj*T_ref, vN, 'LineWidth', 1.5);
    plot(tTraj*T_ref, vU, 'LineWidth', 1.5);
    plot(tTraj*T_ref, vecnorm([stateTraj(:,4), stateTraj(:,5), stateTraj(:,6)], 2, 2)*V_ref, '-', 'LineWidth', 2);
    legend('East', 'North', 'Up', 'Magnitude', 'Location', 'best');
    xlabel('Time s'); ylabel('Velocity m/s'); title('Velocity Profile (Dim)');
    grid on;

% Figure 6: Throttle profile (display limits only)
    figure('Name','Sim Throttle'); hold on;
    thrustDim = aTDimNorm' .* mDim;
    plot(tTraj*T_ref, thrustDim/maxThrustDim, 'b-', 'DisplayName','Constrained (Sim)');
    if hasUC
        thrustDimUC = aTDimNormUC .* mDimUC;
        plot(tTrajUC*T_ref, thrustDimUC/maxThrustDim, 'r-', 'DisplayName','Unconstrained (Sim)');
    end
    yline(1.0,'r--','LineWidth',1,'DisplayName','Max Throttle'); yline(minThrustDim/maxThrustDim,'r--','LineWidth',1,'DisplayName','Min Throttle');
    xlabel('Time s'); ylabel('Throttle Fraction'); title('Time vs Throttle (Sim)'); legend('Location','best'); grid on;

% Figure 7: Mass depletion over time
    figure('Name','Sim Mass Depletion'); hold on; grid on;
    
    % Constrained
    plot(tTraj*T_ref, mDim, 'b-', 'LineWidth', 2, 'DisplayName','Constrained (Sim)');
    
    % Unconstrained
    if hasUC
        plot(tTrajUC*T_ref, mDimUC, 'r-', 'LineWidth', 1.5, 'DisplayName','Unconstrained (Sim)');
    end
    
    xlabel('Time s');
    ylabel('Mass kg');
    title('Vehicle Mass vs Time');
    legend('Location','best');
    
    % --- Summaries for the subtitle ---
    propUsedC  = max(0, mDim(1) - mDim(end));
    if hasUC
        propUsedU = max(0, mDimUC(1) - mDimUC(end));
        subtitle(sprintf(['Final mass and propellant used\n' ...
            'Constrained: m_f = %.1f kg, used = %.1f kg\n' ...
            'Unconstrained: m_f = %.1f kg, used = %.1f kg'], ...
            mDim(end), propUsedC, mDimUC(end), propUsedU));
    else
        subtitle(sprintf('Final mass = %.1f kg, propellant used = %.1f kg', ...
            mDim(end), propUsedC));
    end

% Figure 8: Time vs Altituide
    figure('Name','Sim Alt'); hold on;
    plot(tTraj*T_ref, alt_m/1000, 'LineWidth', 1.5);
    xlabel('Time s'); ylabel('Altitude m'); title('Time vs Altitude (Dimensional)');
    grid on;

    
% Figure 9: Pointing angle vs limit (Sim)
% Compute pointing angle phi from thrust direction in ENU

epsMag = 1e-12;
aMag   = vecnorm([aE aN aU], 2, 2);
thrustU_ENU = [aE aN aU] ./ aMag;          % unit thrust direction
dotUp       = max(-1, min(1, thrustU_ENU(:,3)));   % dot with [0 0 1]
phiSim      = acosd(dotUp);                         % deg

phi0_deg   = optimParams.minPointing;
phiA_deg   = 0.5 * (optimParams.maxTiltAccel) .* (tgoDim.^2);
ThetaSim   = min(180, phi0_deg + phiA_deg);

if hasUC
    % Unconstrained pointing
    aMagUC   = max(aT_norm_ENU_UC, epsMag);
    thrustU_ENU_UC = [aE_UC aN_UC aU_UC] ./ aMagUC;
    dotUpUC  = max(-1, min(1, thrustU_ENU_UC(:,3)));
    phiSimUC = acosd(dotUpUC);
end

figure('Name','Pointing vs Limit (Sim)'); tiledlayout(2,1);

% Top: angles
nexttile; hold on; grid on;
plot(tTraj*T_ref, phiSim, 'b-', 'LineWidth', 1.6, 'DisplayName','\phi (Sim)');
plot(tTraj*T_ref, ThetaSim, 'b--', 'LineWidth', 1.6, 'DisplayName','\Theta limit');
if hasUC
    plot(tTrajUC*T_ref, phiSimUC, 'r-', 'LineWidth', 1.2, 'DisplayName','\phi (Sim UC)');
end
xlabel('Time s'); ylabel('Angle deg');
title('Pointing angle versus limit');
legend('Location','best');

% Bottom: margin
nexttile; hold on; grid on;
plot(tTraj*T_ref, ThetaSim - phiSim, 'b-', 'LineWidth', 1.6, 'DisplayName','Margin = \Theta - \phi');

xlabel('Time s'); ylabel('deg');
title('Pointing margin (positive means within limit)');
legend('Location','best');

else
    % Simulation data not available - only optimization plots were generated
    fprintf('\n=== NOTE: Simulation plots skipped (runSimulation = false) ===\n');
    fprintf('Only optimization-based plots have been generated.\n\n');
end

end