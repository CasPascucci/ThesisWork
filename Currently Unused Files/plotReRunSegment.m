function plotReRunSegment(tSeg, stateSeg, optParams, optCost, aTOptim, mOptim, rdOptim, vdOptim, segAT, refVals, problemParams, nonDimParams, optimParams, flag_thrustGotLimited)
% PLOTRERUNSEGMENT Detailed analysis of a single re-optimization segment.
%   Strictly plots in ENU. 

    %% 1. Data Preparation & Dimensionalization
    L_ref = refVals.L_ref;
    V_ref = refVals.V_ref;
    M_ref = refVals.M_ref;
    A_ref = refVals.A_ref;
    T_ref = refVals.T_ref;
    
    % --- Process The "Flown" Segment (Simulation) ---
    % tSeg comes from ODE45 (N x 1)
    tSegDim = tSeg * T_ref; 
    
    % Simulation State is N x 7. Extract and Transpose to 3 x N for processing
    rSeg = stateSeg(:,1:3)';    % 3 x N
    vSeg = stateSeg(:,4:6)';    % 3 x N
    mSeg = stateSeg(:,7)';      % 1 x N
    
    rSegDim = rSeg * L_ref;     
    vSegDim = vSeg * V_ref;     
    mSegDim = mSeg * M_ref;     
    
    % Guidance Command (segAT) is usually 3 x N from reOptReRun
    if size(segAT, 1) ~= 3
        segAT = segAT'; 
    end
    aSegDim = segAT * A_ref;   % 3 x N
    
    % Convert Simulation to ENU
    % inputs to MCMF2ENU must be (3 x N)
    % isVel = false for Position (enables origin subtract), true for Vectors
    rSegENU = MCMF2ENU(rSegDim, problemParams.landingLatDeg, problemParams.landingLonDeg, true, true);
    vSegENU = MCMF2ENU(vSegDim, problemParams.landingLatDeg, problemParams.landingLonDeg, false, true);
    aSegENU = MCMF2ENU(aSegDim, problemParams.landingLatDeg, problemParams.landingLonDeg, false, true);
    
    EastSeg = rSegENU(1,:); 
    NorthSeg = rSegENU(2,:); 
    UpSeg = rSegENU(3,:);
    
    % --- Process The "Planned" Trajectory (Optimizer) ---
    % The optimizer returns 3 x N usually. We enforce 3 x N.
    tgo0 = optParams(3); % ND
    nodeCount = size(rdOptim, 2); 
    
    % Dimension Check: Ensure 3 x N
    if size(rdOptim, 1) ~= 3; rdOptim = rdOptim'; end
    if size(aTOptim, 1) ~= 3; aTOptim = aTOptim'; end
    if size(vdOptim, 1) ~= 3; vdOptim = vdOptim'; end
    if size(mOptim, 1) > 1; mOptim = mOptim'; end % Ensure 1 x N

    % Dimensionalize
    rdOptimDim = rdOptim * L_ref;
    aTOptimDim = aTOptim * A_ref;
    
    % Convert Plan to ENU
    % isVel = false for Position (X - r0), isVel = true for Accel
    rdOptimENU = MCMF2ENU(rdOptimDim, problemParams.landingLatDeg, problemParams.landingLonDeg, true, true);
    aTOptimENU = MCMF2ENU(aTOptimDim, problemParams.landingLatDeg, problemParams.landingLonDeg, false, true);
    
    EOpt = rdOptimENU(1,:); 
    NOpt = rdOptimENU(2,:); 
    UOpt = rdOptimENU(3,:);
    
    % Time vector for the planned trajectory
    tgospanOpt = linspace(0, tgo0, size(rdOptim,2));
    tPlanRel   = tgo0 - tgospanOpt;           
    tPlanAbs   = tSegDim(1) + tPlanRel*T_ref; 

    %% 2. Figure 1: Spatial Tracking (Zoomed)
    figure('Name', 'Segment: Spatial Tracking');
    tiledlayout(2,2, 'TileSpacing', 'compact');
    
    % -- 3D View --
    nexttile([2,1]); hold on; grid on; axis equal;
    plot3(EOpt/1000, NOpt/1000, UOpt/1000, 'b--', 'LineWidth', 1, 'DisplayName', 'Optimizer Plan');
    plot3(EastSeg/1000, NorthSeg/1000, UpSeg/1000, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Flown Segment');
    plot3(EastSeg(1)/1000, NorthSeg(1)/1000, UpSeg(1)/1000, 'go', 'MarkerFaceColor','g','DisplayName','Seg Start');
    
    xlabel('East (km)'); ylabel('North (km)'); zlabel('Up (km)');
    title(sprintf('Trajectory Segment (t = %.1f to %.1f s)', tSegDim(1), tSegDim(end)));
    legend('Location','best');
    view(3);
    
    % -- Altitude Deviation --
    nexttile; hold on; grid on;
    plot(tPlanAbs, UOpt, 'b.--', 'DisplayName', 'Planned Alt');
    plot(tSegDim, UpSeg, 'r-', 'LineWidth', 2, 'DisplayName', 'Flown Alt');
    ylabel('Altitude (m)'); xticklabels([]);
    title('Vertical Tracking');
    legend('Location','best');
    
    %% 3. Figure 2: Guidance vs Control (Thrust & Throttle)
    figure('Name', 'Segment: Control Dynamics');
    subplot(3,1,1); hold on; grid on;
    
    % Calculate Throttles
    ThrustPlan = vecnorm(aTOptimDim, 2, 1) .* (mOptim * M_ref);
    ThrottlePlan = ThrustPlan / problemParams.maxThrustDim;
    
    ThrustFlown = vecnorm(aSegDim, 2, 1) .* mSegDim; 
    ThrottleFlown = ThrustFlown / problemParams.maxThrustDim;
    
    stairs(tPlanAbs, ThrottlePlan, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Planned Throttle');
    plot(tSegDim, ThrottleFlown, 'r-', 'LineWidth', 2, 'DisplayName', 'Commanded Throttle');
    
    yline(1.0, 'k-');
    yline(problemParams.minThrustDim/problemParams.maxThrustDim, 'k-');
    
    ylabel('Throttle %');
    title('Throttle: Planned vs Commanded');
    legend('Location','best');
    if flag_thrustGotLimited
        subtitle('WARNING: Thrust Saturation Detected in Guidance');
    end
    
    subplot(3,1,2); hold on; grid on;
    % Compare Acceleration Components (ENU)
    plot(tPlanAbs, aTOptimENU(1,:), 'b-', 'LineWidth', 1); 
    plot(tPlanAbs, aTOptimENU(2,:), 'b-', 'LineWidth', 1);
    plot(tPlanAbs, aTOptimENU(3,:), 'b-', 'LineWidth', 1, 'DisplayName', 'Planned Accel (ENU)');
    
    plot(tSegDim, aSegENU(1,:), 'r-', 'LineWidth', 1.5);
    plot(tSegDim, aSegENU(2,:), 'r-', 'LineWidth', 1.5);
    plot(tSegDim, aSegENU(3,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Cmd Accel (ENU)');
    ylabel('Accel (m/s^2)');
    title('Acceleration Vector Components');
    
    subplot(3,1,3); hold on; grid on;
    % Pointing Angle Analysis
    
    % Thrust Vector (Simulation)
    aMagSeg = vecnorm(aSegENU, 2, 1);
    uThrustSeg = aSegENU ./ aMagSeg;
    dotSeg = uThrustSeg(3,:);
    phiSeg = acosd(min(1, max(-1, dotSeg)));
    
    % Thrust Vector (Plan)
    aMagOpt = vecnorm(aTOptimENU, 2, 1);
    uThrustOpt = aTOptimENU ./ aMagOpt;
    dotOpt = uThrustOpt(3,:);
    phiOpt = acosd(min(1, max(-1, dotOpt)));
    
    plot(tPlanAbs, phiOpt, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Planned Tilt');
    plot(tSegDim, phiSeg, 'r-', 'LineWidth', 2, 'DisplayName', 'Flown Tilt');
    
    % Limit Calculation
    tgoSeg = tgo0*T_ref - (tSegDim - tSegDim(1));
    phiLimit = min(180, optimParams.minPointing + 0.5 * (optimParams.maxTiltAccel) * (tgoSeg.^2));
    
    plot(tSegDim, phiLimit, 'k:', 'LineWidth', 2, 'DisplayName', 'Tilt Limit');
    
    ylabel('Tilt Angle (deg)'); xlabel('Time (s)');
    title('Pointing Angle vs Limit');
    legend('Location','best');

end