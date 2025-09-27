function [tTraj, stateTraj, aTList] = closedLoopSim(gamma,kr,tgo0, problemParams, nonDimParams, refVals)

    r0 = nonDimParams.r0ND;
    v0 = nonDimParams.v0ND;
    m0 = nonDimParams.m0ND;
    
    rMoonND = nonDimParams.rMoonND;
    isp = nonDimParams.ispND;

    rfStar = nonDimParams.rfStarND;
    vfStar = nonDimParams.vfStarND;
    afStar = nonDimParams.afStarND;

    X0 = [r0; v0; m0];
    odeoptions = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
     [tTraj, stateTraj] = ode45(@(t, X) trajectory(t, X, gamma, kr, ...
                          tgo0, isp, rMoonND, rfStar, vfStar, afStar,refVals.T_ref), [0, tgo0], X0, odeoptions);

     rState = stateTraj(:,1:3);
     vState = stateTraj(:,4:6);
     mState = stateTraj(:,7);
     tgoState = tgo0 - tTraj;

     numStates = size(rState,1);
     aTList = zeros(3,numStates);

     for idx = 1:numStates
         Ri = rState(idx, :).';
         Vi = vState(idx, :).';
         Mi = mState(idx);
         tgoi = tgoState(idx);
         gi = -(rMoonND^2) * Ri / (norm(Ri)^3);
         if (tgoi * refVals.T_ref) > 0.5
             aT1 = afStar + (( (gamma*kr) / (2*(gamma + 2)) ) - gamma - 1)*(afStar + gi);
             aT2 = ((gamma + 1) / tgoi)*(1 - kr/(gamma + 2))*(vfStar - Vi);
             aT3 = (rfStar - Ri - Vi*tgoi)* kr / tgoi^2;
             aTi = aT1 + aT2 + aT3;
         else
             [rfVirtual, vfVirtual, afVirtual, tgoVirtual] = computeBeyondTerminationTargeting(Ri, Vi, gamma, kr, rfStar, vfStar, afStar, 1/refVals.T_ref, tgoi, gi);
             aT1 = afVirtual + (( (gamma*kr) / (2*(gamma + 2)) ) - gamma - 1)*(afVirtual + gi);
             aT2 = ((gamma + 1) / tgoVirtual)*(1 - kr/(gamma + 2))*(vfVirtual - Vi);
             aT3 = (rfVirtual - Ri - Vi*tgoVirtual)* kr / tgoVirtual^2;
             aTi = aT1 + aT2 + aT3;
         end
         aTList(:,idx) = aTi;
     end
end


%% Functions
function dXdt = trajectory(t, X, gamma, kr, tgo0, isp, rMoonND, rfStar, vfStar, afStar, T_ref)
    r    = X(1:3);
    v    = X(4:6); 
    mass = X(7);

    g = -(rMoonND^2) * r / (norm(r)^3);
    

    tgo  = tgo0 - t;
    
    if (tgo * T_ref) > 0.5    
        aT1 = afStar + (( (gamma*kr) / (2*(gamma + 2)) ) - gamma - 1)*(afStar + g);
        aT2 = ((gamma + 1) / tgo)*(1 - kr/(gamma + 2))*(vfStar - v);
        aT3 = (rfStar - r - v*tgo)* kr / tgo^2;
        aT = aT1 + aT2 + aT3;
    else
        [rfVirtual, vfVirtual, afVirtual, tgoVirtual] = computeBeyondTerminationTargeting(r, v, gamma, kr, rfStar, vfStar, afStar, 1/T_ref, tgo, g);
        aT1 = afVirtual + (( (gamma*kr) / (2*(gamma + 2)) ) - gamma - 1)*(afVirtual + g);
        aT2 = ((gamma + 1) / tgoVirtual)*(1 - kr/(gamma + 2))*(vfVirtual - v);
        aT3 = (rfVirtual - r - v*tgoVirtual)* kr / tgoVirtual^2;
        aT = aT1 + aT2 + aT3;
    end
    F_mag = norm(aT) * mass;

    dm_dt = -F_mag / (isp);
    dXdt = [v; aT + g; dm_dt];
end

function [rfVirtual, vfVirtual, afVirtual, tgoVirtual] = computeBeyondTerminationTargeting(r, v, gamma, kr, rfStar, vfStar, afStar, delta_t, tgo_true, g)
    r = r(:); v = v(:);
    rfStar = rfStar(:); vfStar = vfStar(:); afStar = afStar(:);

    if gamma < 0, error('gamma must be >= 0'); end
    if kr < 2*(gamma + 2), error('kr must be >= 2*(gamma+2)'); end

    tgoVirtual = tgo_true + delta_t;
    if delta_t < 1e-15
        rfVirtual = rfStar; vfVirtual = vfStar; afVirtual = afStar;
        tgoVirtual = tgo_true;
        return;
    end

    gamma1 = gamma;
    gamma2 = kr/(gamma + 2) - 2;

    if abs(gamma1 - gamma2) < 1e-15
        rfVirtual = rfStar; vfVirtual = vfStar; afVirtual = afStar;
        tgoVirtual = tgo_true;
        return;
    end

    tgo = tgoVirtual;

    phi1_bar = -(1/(gamma1 + 1)) * tgo^(gamma1 + 1);
    phi2_bar = -(1/(gamma2 + 1)) * tgo^(gamma2 + 1);
    phi1_hat = tgo^(gamma1 + 2) / ((gamma1 + 1)*(gamma1 + 2));
    phi2_hat = tgo^(gamma2 + 2) / ((gamma2 + 1)*(gamma2 + 2));
    delta = phi1_hat * phi2_bar - phi2_hat * phi1_bar;
    
    
    %fprintf("Value of tgo^(gamma1+1): %.5f, Value of tgo^(gamma2+1): %.5f\n",tgo^(gamma1+1),tgo^(gamma2+1));

    k1r = -phi2_bar / delta;
    k1v = (phi2_bar * tgo + phi2_hat) / delta;
    k1a = -(0.5 * tgo * phi2_bar + phi2_hat) * tgo / delta;

    k2r =  phi1_bar / delta;
    k2v = -(phi1_bar * tgo + phi1_hat) / delta;
    k2a =  (0.5 * tgo * phi1_bar + phi1_hat) * tgo / delta;

    d1 = ((r - 0.5*g*tgo^2)*phi2_bar - (v + g*tgo)*phi2_hat)/delta;
    d2 = -((r - 0.5*g*tgo^2)*phi1_bar - (v + g*tgo)*phi1_hat)/delta;

    phi1_delta = delta_t^gamma1;
    phi2_delta = delta_t^gamma2;
    phi1_bar_delta = -(1/(gamma1 + 1)) * delta_t^(gamma1 + 1);
    phi2_bar_delta = -(1/(gamma2 + 1)) * delta_t^(gamma2 + 1);
    phi1_hat_delta = (1/((gamma1 + 1)*(gamma1 + 2))) * delta_t^(gamma1 + 2);
    phi2_hat_delta = (1/((gamma2 + 1)*(gamma2 + 2))) * delta_t^(gamma2 + 2);

    m11 = k1r * phi1_hat_delta + k2r * phi2_hat_delta + 1;
    m12 = k1v * phi1_hat_delta + k2v * phi2_hat_delta - delta_t;
    m13 = k1a * phi1_hat_delta + k2a * phi2_hat_delta + 0.5 * delta_t^2;

    m21 = k1r * phi1_bar_delta + k2r * phi2_bar_delta;
    m22 = k1v * phi1_bar_delta + k2v * phi2_bar_delta + 1;
    m23 = k1a * phi1_bar_delta + k2a * phi2_bar_delta - delta_t;

    m31 = k1r * phi1_delta + k2r * phi2_delta;
    m32 = k1v * phi1_delta + k2v * phi2_delta;
    m33 = k1a * phi1_delta + k2a * phi2_delta + 1;

    b1 = phi1_hat_delta * d1 + phi2_hat_delta * d2 + 0.5 * g * delta_t^2;
    b2 = phi1_bar_delta * d1 + phi2_bar_delta * d2 - g * delta_t;
    b3 = phi1_delta * d1 + phi2_delta * d2;

    M = zeros(9, 9);
    M(1:3,1:3) = m11*eye(3); M(1:3,4:6) = m12*eye(3); M(1:3,7:9) = m13*eye(3);
    M(4:6,1:3) = m21*eye(3); M(4:6,4:6) = m22*eye(3); M(4:6,7:9) = m23*eye(3);
    M(7:9,1:3) = m31*eye(3); M(7:9,4:6) = m32*eye(3); M(7:9,7:9) = m33*eye(3);

    b_vec = [rfStar - b1; vfStar - b2; afStar - b3];
    %cond(M)
    solution = M \ b_vec;

    rfVirtual = solution(1:3);
    vfVirtual = solution(4:6);
    afVirtual = solution(7:9);
end