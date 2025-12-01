function plotting(tTraj, stateTraj, optParams, optCost, aTOptim, mOptim, rdOptim, vdOptim, aTList, refVals, problemParams, nonDimParams, optimParams, flag_thrustGotLimited, secondaryData, optHistory, ICstates)
% PLOTTING Visualizes trajectory optimization and simulation results.
%   Generates figures comparing the primary run (Initial Plan + Simulation)
%   against a secondary baseline (Comparison Plan + Simulation).

    %% 1. Configuration & Data Processing
    
    % Determine plotting mode
    isReOpt = optimParams.updateOpt;
    if isReOpt
        % Re-Optimization Mode
        labelMainSim = 'Re-Optimized Flight';
        labelSecSim  = 'Static Baseline Flight';
        
        labelMainOpt = 'Initial Plan'; 
        % In Re-Opt, the "Secondary" plan is identical to the "Main" plan at t=0.
        % We disable plotting the secondary plan in the Optimization section to avoid redundancy.
        plotSecondaryOptim = false; 
        
        titleSuffixSim  = ' (Re-Opt vs Static)';
        titleSuffixOpt  = ' (Initial Plan)';
    else
        % Constrained vs Unconstrained Mode
        labelMainSim = 'Constrained Flight';
        labelSecSim  = 'Unconstrained Flight';
        
        labelMainOpt = 'Constrained Plan';
        labelSecOpt  = 'Unconstrained Plan';
        
        plotSecondaryOptim = true; % We want to compare the plans here
        
        titleSuffixSim  = ' (Constrained vs Unconstrained)';
        titleSuffixOpt  = ' (Constrained vs Unconstrained)';
    end

    % Check for simulation data availability
    hasSimData = ~isempty(tTraj) && ~isempty(stateTraj) && ~isempty(aTList);
    
    % Extract Parameters
    tgo0 = optParams(3);
    L_ref = refVals.L_ref;  
    T_ref = refVals.T_ref;  
    A_ref = refVals.A_ref;  
    V_ref = refVals.V_ref;  
    M_ref = refVals.M_ref;  
    rMoon = problemParams.rMoon; 
    [~, ~, U0] = enuBasis(deg2rad(problemParams.landingLatDeg),deg2rad(problemParams.landingLonDeg));
    isp = nonDimParams.ispND;

    % Process Primary Simulation Data
    if hasSimData
        rState = stateTraj(:,1:3);
        rDim = rState * L_ref;
        vState = stateTraj(:,4:6);
        vDim = vState * V_ref;
        mState = stateTraj(:,7);
        mDim = mState * M_ref;
        tgoState = tgo0 - tTraj;
        tgoDim = tgoState * T_ref;

        aTDim = aTList * A_ref;
        aTDimNorm = vecnorm(aTDim,2,1);

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
    
    % Process Secondary Data (if available)
    hasSec = exist('secondaryData','var') && ~isempty(secondaryData);
    if hasSec
        tTrajSec = secondaryData.tTraj;
        rStateSec = secondaryData.stateTraj(:,1:3);
        vStateSec = secondaryData.stateTraj(:,4:6);
        mStateSec = secondaryData.stateTraj(:,7);
        rDimSec = rStateSec * L_ref;
        vDimSec = vStateSec * V_ref;
        mDimSec = mStateSec * M_ref;
        deltaR_Sec   = MCMF2ENU(rDimSec', problemParams.landingLatDeg, problemParams.landingLonDeg, true,  true);
        vDimENU_Sec  = MCMF2ENU(vDimSec', problemParams.landingLatDeg, problemParams.landingLonDeg, false, true);
        EastSec = deltaR_Sec(1,:)';
        NorthSec = deltaR_Sec(2,:)';
        UpSec = deltaR_Sec(3,:)';

        alt_m_Sec = vecnorm(rDimSec,2,2) - rMoon;
        aTDimSec = secondaryData.aTList * A_ref;
        aTDimENU_Sec = MCMF2ENU(aTDimSec, problemParams.landingLatDeg, problemParams.landingLonDeg, false, true);
        aE_Sec = aTDimENU_Sec(1,:)';
        aN_Sec = aTDimENU_Sec(2,:)';
        aU_Sec = aTDimENU_Sec(3,:)';
        aT_norm_ENU_Sec = vecnorm([aE_Sec, aN_Sec, aU_Sec], 2, 2);

        aTDimNormSec = vecnorm(aTDimSec,2,1)';
    end

    maxThrustDim = problemParams.maxThrustDim;
    minThrustDim = problemParams.minThrustDim;

    %% 2. Optimization Results Figures (Planned Trajectories)
    nodeCount = optimParams.nodeCount;
    aTOptim = aTOptim';
    rdOptim = rdOptim';
    vdOptim = vdOptim';
    mOptim = mOptim';
    
    rdOptimTOPO = MCMF2ENU(rdOptim',problemParams.landingLatDeg,problemParams.landingLonDeg,true,false);
    rdOptimTOPODim = rdOptimTOPO*L_ref;

    tgospanOpt = linspace(0,tgo0,nodeCount);
    tspanOpt = tgo0 - tgospanOpt;
    aTNormOpt = vecnorm(aTOptim,2,2);

    EOpt = rdOptimTOPODim(1, :);
    NOpt = rdOptimTOPODim(2, :);
    UOpt = rdOptimTOPODim(3, :);
    
    % Figure: Optimization Throttle Profile
    figure('Name',"Optim Throttle"); hold on;
    thrustDim = aTNormOpt .* mOptim *(refVals.M_ref*refVals.A_ref);
    
    % Plot Main (Blue)
    plot(tspanOpt(2:end)*T_ref, thrustDim(2:end)/problemParams.maxThrustDim, 'b--', 'LineWidth', 1.5, 'DisplayName', labelMainOpt);
    
    % Plot Secondary (Red) - Only if enabled (Constrained vs Unconstrained)
    if plotSecondaryOptim && hasSec && isfield(secondaryData,'aTOptim') && ~isempty(secondaryData.aTOptim) ...
         && isfield(secondaryData,'mOptim')  && ~isempty(secondaryData.mOptim)

        aTOptimSec = secondaryData.aTOptim';
        mOptimSec  = secondaryData.mOptim';
        aTNormOptSec = vecnorm(aTOptimSec, 2, 2);
        thrustDimSec_optim = aTNormOptSec .* mOptimSec * (refVals.M_ref*refVals.A_ref);
        tgo0Sec = secondaryData.optParams(3);

        nodeCountSec = size(aTOptimSec,1);
        tgospanOptSec = linspace(0, tgo0Sec, nodeCountSec);
        tspanOptSec   = tgo0Sec - tgospanOptSec;
    
        plot(tspanOptSec(2:end)*T_ref, thrustDimSec_optim(2:end)/problemParams.maxThrustDim, ...
            'r--', 'LineWidth', 1.3, 'DisplayName', labelSecOpt);
    end
    
    yline(1.0, 'k-', 'LineWidth', 1, 'DisplayName', 'Max Thrust');
    if optimParams.pointingEnabled
        yline(0.95, 'k:', 'LineWidth', 1, 'DisplayName', '95% Max');
    end
    yline(problemParams.minThrustDim/problemParams.maxThrustDim, 'k-', 'LineWidth', 1, 'DisplayName', 'Min Thrust');
    xlabel('Time s'); ylabel('Throttle Fraction'); title(['Planned Throttle Profile' titleSuffixOpt]);
    
    if plotSecondaryOptim && hasSec
        subtitle(sprintf("Main - gamma: %.2f, kr: %.2f, tgo: %.2f s\n Sec - gamma: %.2f, kr: %.2f, tgo: %.2f s", optParams(1),optParams(2),optParams(3)*T_ref,secondaryData.optParams(1),secondaryData.optParams(2),secondaryData.optParams(3)*T_ref));
    else
        subtitle(sprintf("Plan - gamma: %.2f, kr: %.2f, tgo: %.2f s", optParams(1),optParams(2),optParams(3)*T_ref));
    end
    legend('Location','best');

    % Figure: Optimization Range vs Altitude
    rdOptimDim = rdOptim*L_ref;
    rhatOpt = rdOptimDim ./ vecnorm(rdOptimDim,2,2);
    centralOpt = acos(rhatOpt*U0);
    arcLengthOpt = rMoon * centralOpt;
    alt_opt = vecnorm(rdOptimDim,2,2) - rMoon;

    figure('Name','Opt Range vs Altitude'); hold on;

    if plotSecondaryOptim && hasSec && isfield(secondaryData,'rdOptim')
        rdOptimSec = secondaryData.rdOptim';
        rdOptimSecDim = rdOptimSec * L_ref;
        rhatOptSec = rdOptimSecDim ./ vecnorm(rdOptimSecDim,2,2);
        centralOptSec = acos(rhatOptSec*U0);
        arcLengthOptSec = rMoon * centralOptSec;
        alt_optSec = vecnorm(rdOptimSecDim,2,2) - rMoon;
        plot(arcLengthOptSec/1000, alt_optSec/1000, 'r--', 'LineWidth', 1.5, 'DisplayName', labelSecOpt);
    end

    plot(arcLengthOpt/1000, alt_opt/1000, 'b--', 'LineWidth', 1.5, 'DisplayName', labelMainOpt);
    xlabel('Range km'); ylabel('Up km'); title(['Planned Range vs Altitude' titleSuffixOpt]);
    grid on;
    legend('Location','best');
    
    % Figure: Optimization 3D Trajectory
    ETraj = rdOptimTOPODim(1, :);
    NTraj = rdOptimTOPODim(2, :);
    UTraj = rdOptimTOPODim(3, :);
    
    % Generate Glideslope Cone Visualization
    theta_fun = @(u) min(89.9, 45 + 45 * max(0, min(1, (u - 250) ./ 250))); 
    Umin = 0; 
    Umax = min(500, max(UOpt, [], 'omitnan')); 
    if isempty(Umax) || isnan(Umax); Umax = 500; end
    if isempty(Umin) || isnan(Umin); Umin = 0; end
    U_samp = linspace(Umin, Umax, 360);
    theta_deg = theta_fun(U_samp);
    
    cosBound = cosd(theta_deg);
    cosBound = max(cosBound, 1e-3);
    R_samp = U_samp .* sqrt(1 - cosBound.^2) ./ cosBound; 
    Rcap = 3000;
    R_samp = min(R_samp, Rcap);
    
    azimuth = linspace(0, 2*pi, 360);
    [azGrid, altGrid] = meshgrid(azimuth, U_samp);
    [~, rGrid] = meshgrid(azimuth, R_samp);
    eastGrid  = rGrid .* cos(azGrid);
    northGrid = rGrid .* sin(azGrid);
    
    % Generate Plateau Ring
    UPlat = Umax;
    RPlat = Rcap;
    azPlat = linspace(0, 2*pi, 360);
    eastPlat = RPlat * cos(azPlat);
    northPlat = RPlat * sin(azPlat);
    UPlat = UPlat * ones(size(azPlat));
    
    figure('Name','Opt 3D'); hold on; grid on; axis equal;
    surf(eastGrid/1000, northGrid/1000, altGrid/1000,altGrid/1000,'EdgeAlpha',0.15,'FaceAlpha',0.4,'MeshStyle','row','LineWidth',0.8);
    plot3(ETraj/1000, NTraj/1000, UTraj/1000, 'b-', 'LineWidth', 2, 'MarkerSize',10, 'DisplayName', labelMainOpt);
    plot3(eastPlat/1000, northPlat/1000, UPlat/1000, 'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
   

    if plotSecondaryOptim && hasSec && isfield(secondaryData,'rdOptim')
        rdOptimSec = secondaryData.rdOptim';
        rdOptimSecTOPO = MCMF2ENU(rdOptimSec',problemParams.landingLatDeg,problemParams.landingLonDeg,true,false);
        rdOptimSecTOPODim = rdOptimSecTOPO * L_ref;
        ETrajSec = rdOptimSecTOPODim(1, :);
        NTrajSec = rdOptimSecTOPODim(2, :);
        UTrajSec = rdOptimSecTOPODim(3, :);
        plot3(ETrajSec/1000, NTrajSec/1000, UTrajSec/1000, 'r--', 'LineWidth', 1.5, 'DisplayName', labelSecOpt);
    end
    
    xlabel('East (km)');
    ylabel('North (km)');
    zlabel('Up (km)');
    title(['Planned 3D Trajectory' titleSuffixOpt]);
    view(80, 15); camproj orthographic;
    zlim([0,2]);
    xlim([-3,3]);
    ylim([-3,3]);
    axis square;
    legend('Location','bestoutside');

    % Figure: Optimization Velocity Profile
    vdOptimDim = vdOptim * V_ref;
    vdOptimTOPO = MCMF2ENU(vdOptimDim',problemParams.landingLatDeg,problemParams.landingLonDeg,false,true);
    figure('Name','Opt Vel TOPO'); hold on;
    
    % Main Trace (Blue)
    plot(tspanOpt*T_ref, vdOptimTOPO(1,:), 'b-', 'LineWidth', 1.0, 'DisplayName', [labelMainOpt ' East']);
    plot(tspanOpt*T_ref, vdOptimTOPO(2,:), 'b--', 'LineWidth', 1.0, 'DisplayName', [labelMainOpt ' North']);
    plot(tspanOpt*T_ref, vdOptimTOPO(3,:), 'b:', 'LineWidth', 1.0, 'DisplayName', [labelMainOpt ' Up']);
    plot(tspanOpt*T_ref, vecnorm(vdOptimTOPO, 2, 1), 'b-', 'LineWidth', 2, 'DisplayName', [labelMainOpt ' Mag']);
    
    % Secondary Trace (Red - Magnitude Only)
    if plotSecondaryOptim && hasSec && isfield(secondaryData,'vdOptim')
        vdOptimSec = secondaryData.vdOptim';
        vdOptimSecDim = vdOptimSec * V_ref;
        vdOptimSecTOPO = MCMF2ENU(vdOptimSecDim',problemParams.landingLatDeg,problemParams.landingLonDeg,false,true);
        tspanOptSec = tgo0Sec - linspace(0, tgo0Sec, size(vdOptimSec,1));
        plot(tspanOptSec*T_ref, vecnorm(vdOptimSecTOPO, 2, 1), 'r--', 'LineWidth', 2, 'DisplayName', [labelSecOpt ' Mag']);
    end
    
    legend('Location', 'best');
    xlabel('Time s'); ylabel('Velocity m/s'); title(['Planned Velocity Profile' titleSuffixOpt]);
    grid on;
    xlim([0, ceil(max(tspanOpt*T_ref)/100)*100]);

    % Figure: Optimization TOPO Acceleration
    aTOptimDim = aTOptim * A_ref;
    aTOptimTOPO = MCMF2ENU(aTOptimDim',problemParams.landingLatDeg,problemParams.landingLonDeg,false,true);
    figure('Name','Opt TOPO Accel'); hold on;
    
    % Main Trace (Blue)
    plot(tspanOpt(2:end)*T_ref, aTOptimTOPO(1,2:end), 'b-', 'LineWidth', 1.0, 'DisplayName', [labelMainOpt ' East']);
    plot(tspanOpt(2:end)*T_ref, aTOptimTOPO(2,2:end), 'b--', 'LineWidth', 1.0, 'DisplayName', [labelMainOpt ' North']);
    plot(tspanOpt(2:end)*T_ref, aTOptimTOPO(3,2:end), 'b:', 'LineWidth', 1.0, 'DisplayName', [labelMainOpt ' Up']);
    plot(tspanOpt(2:end)*T_ref, vecnorm(aTOptimTOPO(:,2:end),2,1), 'b-', 'LineWidth', 2, 'DisplayName', [labelMainOpt ' Mag']);
    
    % Secondary Trace (Red - Magnitude Only)
    if plotSecondaryOptim && hasSec && isfield(secondaryData,'aTOptim')
        aTOptimSec = secondaryData.aTOptim';
        aTOptimSecDim = aTOptimSec * A_ref;
        aTOptimSecTOPO = MCMF2ENU(aTOptimSecDim',problemParams.landingLatDeg,problemParams.landingLonDeg,false,true);
        tspanOptSec = tgo0Sec - linspace(0, tgo0Sec, size(aTOptimSec,1));
        plot(tspanOptSec(2:end)*T_ref, vecnorm(aTOptimSecTOPO(:,2:end),2,1), 'r--', 'LineWidth', 2, 'DisplayName', [labelSecOpt ' Mag']);
    end

    plot(tspanOpt(1)*T_ref,norm(nonDimParams.afStarND)*A_ref,'.','MarkerSize',10, 'HandleVisibility', 'off'); 
    plot(tspanOpt(1)*T_ref,A_ref,'x','MarkerSize',10, 'HandleVisibility', 'off'); 
    
    legend('Location', 'best');
    xlabel('Time s'); ylabel('Accel m/s^2'); title(['Planned Accel Profile (TOPO)' titleSuffixOpt]);
    grid on;

    % Figure: Optimization MCMF Acceleration
    aTOptimDim = aTOptimDim';
    figure('Name','Opt MCMF Accel'); hold on;
    
    % Main Trace (Blue)
    plot(tspanOpt(2:end)*T_ref, aTOptimDim(1,2:end), 'b-', 'LineWidth', 1.0, 'DisplayName', [labelMainOpt ' X']);
    plot(tspanOpt(2:end)*T_ref, aTOptimDim(2,2:end), 'b--', 'LineWidth', 1.0, 'DisplayName', [labelMainOpt ' Y']);
    plot(tspanOpt(2:end)*T_ref, aTOptimDim(3,2:end), 'b:', 'LineWidth', 1.0, 'DisplayName', [labelMainOpt ' Z']);
    plot(tspanOpt(2:end)*T_ref, vecnorm(aTOptimDim(:,2:end),2,1), 'b-', 'LineWidth', 2, 'DisplayName', [labelMainOpt ' Mag']);
    
    % Secondary Trace (Red - Magnitude Only)
    if plotSecondaryOptim && hasSec
        aTOptimSecDim_MCMF = (secondaryData.aTOptim' * A_ref)'; 
        plot(tspanOptSec(2:end)*T_ref, vecnorm(aTOptimSecDim_MCMF(:,2:end),2,1), 'r--', 'LineWidth', 2, 'DisplayName', [labelSecOpt ' Mag']);
    end

    plot(tspanOpt(1)*T_ref,norm(nonDimParams.afStarND)*A_ref,'.','MarkerSize',10, 'HandleVisibility', 'off'); 
    plot(tspanOpt(1)*T_ref,A_ref,'x','MarkerSize',10, 'HandleVisibility', 'off'); 
    
    legend('Location', 'best');
    xlabel('Time s'); ylabel('Accel m/s^2'); title(['Planned Accel Profile (MCMF)' titleSuffixOpt]);
    grid on;

    %% 3. Simulation Results Figures (Flown Trajectories)
    if hasSimData
        % Figure: 3D Trajectory (Simulation)
        idx2km = find(alt_m <= 2000, 1, 'first');
        if isempty(idx2km)
            NtotalSim = numel(alt_m);
            idx2km = max(1, NtotalSim - floor(0.20*NtotalSim) + 1);
        end
        idx_sim = idx2km:numel(alt_m);

        ETrajSim = East(idx_sim);
        NTrajSim = North(idx_sim);
        UTrajSim = Up(idx_sim);

        if hasSec
            idx2kmSec = find(alt_m_Sec <= 2000, 1, 'first');
            if isempty(idx2kmSec)
                NtotalSimSec = numel(alt_m_Sec);
                idx2kmSec = max(1, NtotalSimSec - floor(0.20*NtotalSimSec) + 1);
            end
            idx_simSec = idx2kmSec:numel(alt_m_Sec);

            ETrajSimSec = EastSec(idx_simSec);
            NTrajSimSec = NorthSec(idx_simSec);
            UTrajSimSec = UpSec(idx_simSec);
        end

        theta_fun = @(u) min(89.9, 45 + 45*max(0, min(1, (u-250)./250)));
        UminEnv = min(UTrajSim);                       
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

        UPlat = UmaxEnv; RPlat = Rcap;
        azPlat = linspace(0, 2*pi, 360);
        eastPlat  = RPlat * cos(azPlat);
        northPlat = RPlat * sin(azPlat);
        UPlat     = UPlat * ones(size(azPlat));

        figure('Name','Sim 3D'); hold on; grid on; axis equal;
        surf(eastGrid/1000, northGrid/1000, altGrid/1000,altGrid/1000,'EdgeAlpha',0.15,'FaceAlpha',0.4,'MeshStyle','row','LineWidth',0.8);
        plot3(eastPlat/1000, northPlat/1000, UPlat/1000, 'k--', 'LineWidth', 0.8);
        plot3(ETrajSim/1000, NTrajSim/1000, UTrajSim/1000, 'b-', 'LineWidth', 2, 'DisplayName', labelMainSim);
        if hasSec
            plot3(ETrajSimSec/1000, NTrajSimSec/1000, UTrajSimSec/1000, 'r-', 'LineWidth', 2, 'DisplayName', labelSecSim);
        end

        xlabel('East (km)'); ylabel('North (km)'); zlabel('Up (km)');
        title(['3D Trajectory (Sim) Final 2 KM' titleSuffixSim]);
        view(80, 15); camproj orthographic;
        zlim([0, 2]); xlim([-3,3]); ylim([-3,3]); axis square;
        legend('Location','bestoutside');

        % Figure: Range vs Altitude (Simulation)
        rhat      = rDim ./ vecnorm(rDim,2,2);
        central   = acos(rhat*U0);
        arcLength = rMoon * central;
        
        figure('Name','Sim Range vs Altitude'); hold on;
        plot(arcLength/1000, alt_m/1000, 'b-', 'DisplayName', labelMainSim);
        
        if hasSec
            rhatSec      = rDimSec ./ vecnorm(rDimSec,2,2);
            centralSec   = acos(rhatSec*U0);
            arcLengthSec = rMoon * centralSec;
            plot(arcLengthSec/1000, alt_m_Sec/1000, 'r-', 'DisplayName', labelSecSim);
        end

        xlabel('Range (km)'); ylabel('Up (km)');
        title(['Range vs Altitude (Sim)' titleSuffixSim]); grid on;
        legend('Location','best');

        % Figure: TOPO Acceleration (Simulation)
        figure('Name','Sim TOPO Accel'); hold on;
        plot(tTraj*T_ref, aE, 'b-', 'LineWidth', 1.0, 'DisplayName', [labelMainSim ' East']);
        plot(tTraj*T_ref, aN, 'b--', 'LineWidth', 1.0, 'DisplayName', [labelMainSim ' North']);
        plot(tTraj*T_ref, aU, 'b:', 'LineWidth', 1.0, 'DisplayName', [labelMainSim ' Up']);
        plot(tTraj*T_ref, aT_norm_ENU, 'b-', 'LineWidth', 2, 'DisplayName', [labelMainSim ' Mag']);
        
        if hasSec
            plot(tTrajSec*T_ref, aT_norm_ENU_Sec, 'r-', 'LineWidth', 2, 'DisplayName', [labelSecSim ' Mag']);
        end
        
        plot(tTraj(end)*T_ref,norm(nonDimParams.afStarND)*A_ref,'.','MarkerSize',10, 'HandleVisibility','off'); 
        plot(tTraj(end)*T_ref,A_ref,'x','MarkerSize',10, 'HandleVisibility','off'); 
        
        legend('Location', 'best');
        xlabel('Time s'); ylabel('Accel m/s^2'); title(['Commanded Accel Profile (TOPO)' titleSuffixSim]);
        grid on;

        % Figure: MCMF Acceleration (Simulation)
        figure('Name','Sim MCMF Accel'); hold on;
        plot(tTraj*T_ref, aTDim(1,:), 'b-', 'LineWidth', 1.0, 'DisplayName', [labelMainSim ' X']);
        plot(tTraj*T_ref, aTDim(2,:), 'b--', 'LineWidth', 1.0, 'DisplayName', [labelMainSim ' Y']);
        plot(tTraj*T_ref, aTDim(3,:), 'b:', 'LineWidth', 1.0, 'DisplayName', [labelMainSim ' Z']);
        plot(tTraj*T_ref, aTDimNorm, 'b-', 'LineWidth', 2, 'DisplayName', [labelMainSim ' Mag']);
        
        if hasSec
            plot(tTrajSec*T_ref, aTDimNormSec, 'r-', 'LineWidth', 2, 'DisplayName', [labelSecSim ' Mag']);
        end
        
        legend('Location', 'best');
        xlabel('Time s'); ylabel('Accel m/s^2'); title(['Commanded Accel Profile (MCMF)' titleSuffixSim]);
        if flag_thrustGotLimited
            subtitle("Thrust is being Throttled");
        end
        grid on;

        % Figure: Velocity Profile (Simulation)
        figure('Name','Sim Vel'); hold on;
        plot(tTraj*T_ref, vE, 'b-', 'LineWidth', 1.0, 'DisplayName', [labelMainSim ' East']);
        plot(tTraj*T_ref, vN, 'b--', 'LineWidth', 1.0, 'DisplayName', [labelMainSim ' North']);
        plot(tTraj*T_ref, vU, 'b:', 'LineWidth', 1.0, 'DisplayName', [labelMainSim ' Up']);
        plot(tTraj*T_ref, vecnorm(vDim, 2, 2), 'b-', 'LineWidth', 2, 'DisplayName', [labelMainSim ' Mag']);

        if hasSec
            plot(tTrajSec*T_ref, vecnorm(vDimSec, 2, 2), 'r-', 'LineWidth', 2, 'DisplayName', [labelSecSim ' Mag']);
        end

        legend('Location', 'best');
        xlabel('Time s'); ylabel('Velocity m/s'); title(['Velocity Profile (Dim)' titleSuffixSim]);
        grid on;
        
        % Figure: Throttle Profile (Simulation)
        figure('Name','Sim Throttle'); hold on;
        thrustDim = aTDimNorm' .* mDim;
        plot(tTraj*T_ref, thrustDim/maxThrustDim, 'b-', 'DisplayName', labelMainSim);
        if hasSec
            thrustDimSec = aTDimNormSec .* mDimSec;
            plot(tTrajSec*T_ref, thrustDimSec/maxThrustDim, 'r-', 'DisplayName', labelSecSim);
        end
        yline(1.0,'k--','LineWidth',1,'DisplayName','Max Throttle'); 
        yline(minThrustDim/maxThrustDim,'k--','LineWidth',1,'DisplayName','Min Throttle');
        xlabel('Time s'); ylabel('Throttle Fraction'); title(['Time vs Throttle (Sim)' titleSuffixSim]); 
        legend('Location','best'); grid on;

        % Figure: Mass Depletion
        figure('Name','Sim Mass Depletion'); hold on; grid on;
        plot(tTraj*T_ref, mDim, 'b-', 'LineWidth', 2, 'DisplayName', labelMainSim);
        if hasSec
            plot(tTrajSec*T_ref, mDimSec, 'r-', 'LineWidth', 1.5, 'DisplayName', labelSecSim);
        end
        xlabel('Time s'); ylabel('Mass kg'); title('Vehicle Mass vs Time');
        legend('Location','best');
        
        propUsedC  = max(0, mDim(1) - mDim(end));
        if hasSec
            propUsedU = max(0, mDimSec(1) - mDimSec(end));
            subtitle(sprintf(['%s Used: %.1f kg\n%s Used: %.1f kg'], ...
                labelMainSim, propUsedC, labelSecSim, propUsedU));
        else
            subtitle(sprintf('Propellant used = %.1f kg', propUsedC));
        end

        % Figure: Time vs Altitude
        figure('Name','Sim Alt'); hold on;
        plot(tTraj*T_ref, alt_m/1000, 'b-', 'LineWidth', 1.5, 'DisplayName', labelMainSim);
        if hasSec
            plot(tTrajSec*T_ref, alt_m_Sec/1000, 'r-', 'LineWidth', 1.5, 'DisplayName', labelSecSim);
        end
        xlabel('Time s'); ylabel('Altitude m'); title(['Time vs Altitude (Dimensional)' titleSuffixSim]);
        legend('Location','best');
        grid on;

        % Figure: Pointing Angle Analysis
        if optimParams.pointingEnabled
            epsMag = 1e-12;
            aMag   = vecnorm([aE aN aU], 2, 2);
            thrustU_ENU = [aE aN aU] ./ aMag;          
            dotUp       = max(-1, min(1, thrustU_ENU(:,3)));   
            phiSim      = acosd(dotUp);                         
            
            phi0_deg   = optimParams.minPointing;
            phiA_deg   = 0.5 * (optimParams.maxTiltAccel) .* (tgoDim.^2);
            ThetaSim   = min(180, phi0_deg + phiA_deg);
            
            if hasSec
                aMagSec   = max(aT_norm_ENU_Sec, epsMag);
                thrustU_ENU_Sec = [aE_Sec aN_Sec aU_Sec] ./ aMagSec;
                dotUpSec  = max(-1, min(1, thrustU_ENU_Sec(:,3)));
                phiSimSec = acosd(dotUpSec);
            end
            
            figure('Name','Pointing vs Limit (Sim)'); tiledlayout(2,1);
            
            nexttile; hold on; grid on;
            plot(tTraj*T_ref, phiSim, 'b-', 'LineWidth', 1.6, 'DisplayName',[labelMainSim ' \phi']);
            plot(tTraj*T_ref, ThetaSim, 'k--', 'LineWidth', 1.6, 'DisplayName','\Theta limit');
            if hasSec
                plot(tTrajSec*T_ref, phiSimSec, 'r-', 'LineWidth', 1.2, 'DisplayName',[labelSecSim ' \phi']);
            end
            xlabel('Time s'); ylabel('Angle deg');
            title('Pointing angle versus limit');
            legend('Location','best');
            
            nexttile; hold on; grid on;
            plot(tTraj*T_ref, ThetaSim - phiSim, 'b-', 'LineWidth', 1.6, 'DisplayName','Margin (Main)');
            xlabel('Time s'); ylabel('deg');
            title('Pointing margin (positive means within limit)');
            legend('Location','best');
        end
    else
        fprintf('\n=== NOTE: Simulation plots skipped (runSimulation = false) ===\n');
    end

    %% 4. Parameter History (Re-Optimization Only)
    if ~isempty(optHistory) && exist('ICstates', 'var') && ~isempty(ICstates)
        figure('Name', 'ReOpt Parameter History');
        tiledlayout(5, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        
        t_elapsedND = table2array(optHistory(:,1));
        x_opt = t_elapsedND .* refVals.T_ref;
        gamma_hist = table2array(optHistory(:,2));
        gamma2_hist = table2array(optHistory(:,3));
        tgoSolved_ND = table2array(optHistory(:,5));
        tgoSolved = tgoSolved_ND .* refVals.T_ref;
        
        ICstates_array = table2array(ICstates);
        nPoints = size(optHistory, 1);
        c1_norm = zeros(nPoints, 1);
        c2_norm = zeros(nPoints, 1);
        
        rfStar = nonDimParams.rfStarND;
        vfStar = nonDimParams.vfStarND;
        afStar = nonDimParams.afStarND;
        gConst = nonDimParams.gConst;
        
        for i = 1:nPoints
            r_i = ICstates_array(i, 1:3)';
            v_i = ICstates_array(i, 4:6)';
            [c1, c2] = calculateCoeffs(r_i, v_i, tgoSolved_ND(i), gamma_hist(i), gamma2_hist(i), afStar, rfStar, vfStar, gConst);
            c1_norm(i) = norm(c1);
            c2_norm(i) = norm(c2);
        end
        
        nexttile; plot(x_opt, gamma_hist, 'LineWidth', 2); title('\gamma_1'); grid on;
        nexttile; plot(x_opt, gamma2_hist, 'LineWidth', 2); title('\gamma_2'); grid on;
        nexttile; plot(x_opt, tgoSolved, 'LineWidth', 2); title('t_{go} Solved'); grid on;
        nexttile; plot(x_opt, c1_norm, 'LineWidth', 2); title('||c_1||'); grid on;
        nexttile; plot(x_opt, c2_norm, 'LineWidth', 2); title('||c_2||'); grid on;
        
        sgtitle('Optimization Parameters Through Reoptimization');
    end

end