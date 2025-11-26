function [tTraj, stateTraj, aTList, flag_thrustGotLimited, optHistory, ICstates, exitFlags] = simReOpt(gamma0,gamma20,tgo0, problemParams, nonDimParams, refVals, delta_t, optimizationParams, betaParam, verboseOutput)

kr0 = (gamma20+2)*(gamma0+2);

initialNodeCount = optimizationParams.nodeCount;

% Original IC's
r0 = nonDimParams.r0ND;
v0 = nonDimParams.v0ND;
m0 = nonDimParams.m0ND;
rfStar = nonDimParams.rfStarND;
vfStar = nonDimParams.vfStarND;
afStar = nonDimParams.afStarND;
isp = nonDimParams.ispND;
rMoonND = nonDimParams.rMoonND;
gGuidance = nonDimParams.gConst;

X0 = [r0; v0; m0];

updateFreqND = optimizationParams.updateFreq / refVals.T_ref;
updateStopND = optimizationParams.updateStop / refVals.T_ref;


% Setup Results
t_elapsed = 0;
gamma = gamma0;
gamma2 = gamma20;
kr = kr0;
tgo = tgo0;

tTraj = [];
stateTraj = [];
ICstates = [];
ICstates(:,1) = X0;
ICcounter = 2;
aTList = [];
flag_thrustGotLimited = false;
optHistory = [t_elapsed, gamma, gamma2, kr, tgo];
exitFlags = [];
exitFlags(1) = 99;

odeoptions = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

if verboseOutput
    fprintf('\n=== Starting Re-Optimization Simulation ===\n');
    fprintf('Initial tgo: %.2f s\n', tgo0 * refVals.T_ref);
    fprintf('Update frequency: %.2f s\n', optimizationParams.updateFreq);
    fprintf('Update stop time: %.2f s (%.2f s remaining)\n', ...
        optimizationParams.updateStop, tgo0 * refVals.T_ref);
    fprintf('==========================================\n\n');
end
segmentCount = 0;
minTime = 0.2/refVals.T_ref;

% Main loop
    while tgo > updateStopND
        segmentCount = segmentCount + 1;
        tgo_remaining = tgo - updateStopND;
        t_segment = min(tgo,min(updateFreqND, tgo_remaining));
    
        if t_segment <= 0
            break;
        end
    
        t_start = t_elapsed;
        t_end = t_elapsed + t_segment;
    
        if verboseOutput
            fprintf('--- Segment %d ---\n', segmentCount);
            fprintf('Time: %.2f to %.2f s (elapsed)\n', t_start * refVals.T_ref, t_end * refVals.T_ref);
            fprintf('Tgo at start: %.2f s\n', tgo * refVals.T_ref);
            fprintf('Current params: gamma=%.4f, gamma2=%.4f, kr=%.4f\n', gamma, gamma2, kr);
        end
    
        [tSeg, stateSeg] = ode45(@(t, X) trajectorySegmentODE(t, X, t_elapsed, gamma, kr, tgo, ...
            isp, rMoonND, rfStar, vfStar, afStar, gGuidance, nonDimParams), ...
            [t_start, t_end], X0, odeoptions);
    
        tTraj = [tTraj;tSeg];
        stateTraj = [stateTraj;stateSeg];
    
        numPoints = length(tSeg);
        segAT = zeros(3,numPoints);
        for i = 1:numPoints
            r = stateSeg(i,1:3)';
            v = stateSeg(i,4:6)';
            m = stateSeg(i,7);
            t_since_reopt = tSeg(i) - t_elapsed;
            tgo_now = tgo - t_since_reopt;
            [aTi, limited] = computeAccel(r, v, m, gamma, kr, tgo_now, rfStar, vfStar, afStar, gGuidance, nonDimParams);
            segAT(:, i) = aTi;
            flag_thrustGotLimited = flag_thrustGotLimited || limited;
        end
        aTList = [aTList, segAT];
    
        X0 = stateSeg(end,:)';
        ICstates(:,ICcounter) = X0;
        ICcounter = ICcounter + 1;

        t_elapsed = tSeg(end);
    
        tgo = tgo - t_segment;
    
        if verboseOutput
            fprintf('Segment complete. New tgo: %.2f s\n', tgo * refVals.T_ref);
        end
        if tgo <= updateStopND
            if verboseOutput
                fprintf('Reached update stop time. Exiting re-optimization loop.\n\n');
            end
            break;
        end
    
    
        newNonDimParams = nonDimParams;
        newNonDimParams.r0ND = X0(1:3);
        newNonDimParams.v0ND = X0(4:6);
        newNonDimParams.m0ND = X0(7);
        
        paramsX0 = [gamma, gamma2, tgo];
        
        if verboseOutput
            fprintf('\nRe-optimizing at t=%.2f s (tgo=%.2f s)...\n', t_elapsed * refVals.T_ref, tgo * refVals.T_ref);
        end
        
        % Scale node count proportional to tgo
        newNodeCount = ceil(initialNodeCount * (tgo / tgo0));
        % Ensure odd number
        if mod(newNodeCount, 2) == 0
            newNodeCount = newNodeCount + 1;
        end
        
        optimizationParams.nodeCount = newNodeCount;

        if verboseOutput
            fprintf('Scaled node count to: %d\n', optimizationParams.nodeCount);
        end

        [optParams, ~, ~, ~, ~, ~, exitflag] = optimizationLoop(paramsX0, betaParam, ...
            problemParams, newNonDimParams, optimizationParams, refVals, delta_t, false, false);

        if exitflag > 0
            % Check for unreasonable jumps
            gamma_change = abs(optParams(1) - gamma);
            if gamma_change > 0.5  % More than 30% change
                fprintf('WARNING: Large gamma jump (%.4f -> %.4f). Rejecting.\n', gamma, optParams(1));
                % Don't update parameters
            else
                gamma = optParams(1);
                gamma2 = optParams(2);
                kr = (gamma2+2)*(gamma+2);
                tgo = optParams(3);
            end
        else
            % Keep previous parameters if optimization failed
            if verboseOutput
                fprintf('WARNING: Optimization failed (exitflag=%d). Keeping previous parameters.\n', exitflag);
            end
        end
        
        optHistory = [optHistory; [t_elapsed, gamma, gamma2, kr, tgo]];
        exitFlags = [exitFlags; exitflag];
        
        if verboseOutput
            fprintf('New parameters: gamma=%.4f, gamma2=%.4f, kr=%.4f, tgo=%.2f s (exitflag=%d)\n\n', ...
                gamma, gamma2, kr, tgo * refVals.T_ref, exitflag);
        end
    end
    
    if tgo > 0
        segmentCount = segmentCount + 1;
        
        if verboseOutput
            fprintf('--- Final Segment %d (No Re-opt) ---\n', segmentCount);
            fprintf('Integrating final %.2f s to landing\n', tgo * refVals.T_ref);
        end
        
        t_start = t_elapsed;

        
        
        t_freeze = t_elapsed + (tgo -minTime)
        t_end = t_elapsed + tgo;
        [tSeg1, stateSeg1] = ode45(@(t, X) trajectorySegmentODE(t, X, t_elapsed, gamma, kr, tgo, ...
        isp, rMoonND, rfStar, vfStar, afStar, gGuidance, nonDimParams), ...
        [t_start, t_freeze], X0, odeoptions);
    
        % Compute final thrust acceleration to freeze
        r_freeze = stateSeg1(end, 1:3)';
        v_freeze = stateSeg1(end, 4:6)';
        m_freeze = stateSeg1(end, 7);
        [aT_frozen, ~] = computeAccel(r_freeze, v_freeze, m_freeze, gamma, kr, ...
            minTime, rfStar, vfStar, afStar, gGuidance, nonDimParams);
        
        % Second part: frozen thrust acceleration
        X0_freeze = stateSeg1(end, :)';
        [tSeg2, stateSeg2] = ode45(@(t, X) trajectoryFrozenThrust(t, X, aT_frozen, ...
            isp, rMoonND), [t_freeze, t_end], X0_freeze, odeoptions);
        
        % Combine segments
        tSeg = [tSeg1; tSeg2(2:end)];  % Avoid duplicate point
        stateSeg = [stateSeg1; stateSeg2(2:end, :)];


    
        tTraj = [tTraj;tSeg];
        stateTraj = [stateTraj; stateSeg];
        numPoints = length(tSeg);
        segAT = zeros(3, numPoints);
        for i = 1:numPoints
            r = stateSeg(i, 1:3)';
            v = stateSeg(i, 4:6)';
            m = stateSeg(i, 7);
            t_since_reopt = tSeg(i) - t_elapsed;
            tgo_now = tgo - t_since_reopt;
            
            % Use active guidance until threshold, then freeze
            if tgo_now > minTime
                [aTi, limited] = computeAccel(r, v, m, gamma, kr, tgo_now, rfStar, vfStar, afStar, gGuidance, nonDimParams);
                flag_thrustGotLimited = flag_thrustGotLimited || limited;
            else
                % Use frozen value (already computed above)
                aTi = aT_frozen;
            end
            segAT(:, i) = aTi;
        end
        aTList = [aTList, segAT];
        if verboseOutput
            fprintf('Final segment complete.\n\n');
        end
    end
    if verboseOutput
        fprintf('=== Simulation Complete ===\n');
        fprintf('Total elapsed time: %.2f s\n', tTraj(end) * refVals.T_ref);
        fprintf('Total segments: %d\n', segmentCount);
        fprintf('Re-optimizations: %d\n', size(optHistory, 1) - 1);
        fprintf('Thrust limited: %s\n', string(flag_thrustGotLimited));
        fprintf('===========================\n\n');
    end
end


function [aTi, limited] = computeAccel(r, v, m, gamma, kr, tgo, rfStar, vfStar, afStar, gGuidance, nonDimParams)
    % Compute thrust acceleration with limits
    minAccel = nonDimParams.minThrustND / m;
    maxAccel = nonDimParams.maxThrustND / m;
    
    % Guidance law
    aT1 = gamma * (kr/(2*gamma + 4) - 1) * afStar;
    aT2 = (gamma*kr/(2*gamma+4) - gamma - 1) * gGuidance;
    aT3 = ((gamma+1)/tgo) * (1 - kr/(gamma+2)) * (vfStar - v);
    aT4 = (kr/tgo^2) * (rfStar - r - v*tgo);
    aTi = aT1 + aT2 + aT3 + aT4;
    
    % Apply thrust limits
    normAT = norm(aTi);
    if normAT > maxAccel
        aTi = aTi / normAT * maxAccel;
        limited = true;
    elseif normAT < minAccel
        aTi = aTi / normAT * minAccel;
        limited = true;
    else
        limited = false;
    end
end

function dXdt = trajectorySegmentODE(t, X, t_reopt_start, gamma, kr, tgo_at_reopt, isp, rMoonND, rfStar, vfStar, afStar, gGuidance, nonDimParams)
    % ODE function for trajectory segment
    % t - current time (elapsed from simulation start)
    % t_reopt_start - time when current parameters were set
    % tgo_at_reopt - time-to-go when current parameters were set
    
    r = X(1:3);
    v = X(4:6);
    m = X(7);
    
    % Calculate current time-to-go
    t_since_reopt = t - t_reopt_start;
    tgo = tgo_at_reopt - t_since_reopt;
    
    % Ensure tgo doesn't go negative
    if tgo < 0.001
        tgo = 0.001;
    end
    
    % Gravitational acceleration
    g = -(rMoonND^2) * r / (norm(r)^3);
    
    % Compute thrust acceleration
    [aT, ~] = computeAccel(r, v, m, gamma, kr, tgo, rfStar, vfStar, afStar, gGuidance, nonDimParams);
    
    % Thrust magnitude and mass flow rate
    F_mag = norm(aT) * m;
    dm_dt = -F_mag / isp;
    
    % State derivatives
    dXdt = [v; aT + g; dm_dt];
end

function dXdt = trajectoryFrozenThrust(t, X, aT_frozen, isp, rMoonND)
    % ODE function for trajectory with frozen thrust acceleration
    % Used to avoid singularity when tgo approaches zero
    
    r = X(1:3);
    v = X(4:6);
    m = X(7);
    
    % Gravitational acceleration
    g = -(rMoonND^2) * r / (norm(r)^3);
    
    % Use frozen thrust acceleration
    aT = aT_frozen;
    
    % Thrust magnitude and mass flow rate
    F_mag = norm(aT) * m;
    dm_dt = -F_mag / isp;
    
    % State derivatives
    dXdt = [v; aT + g; dm_dt];
end