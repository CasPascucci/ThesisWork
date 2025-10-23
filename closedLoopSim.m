function [tTraj, stateTraj, aTList, flag_thrustGotLimited] = closedLoopSim(gamma,kr,tgo0, problemParams, nonDimParams, refVals, delta_t)

    
    r0 = nonDimParams.r0ND;
    v0 = nonDimParams.v0ND;
    m0 = nonDimParams.m0ND;
    landingLat = problemParams.landingLatDeg;
    landingLon = problemParams.landingLonDeg;

    rMoonND = nonDimParams.rMoonND;
    isp = nonDimParams.ispND;

    rfStar = nonDimParams.rfStarND;
    vfStar = nonDimParams.vfStarND;
    afStar = nonDimParams.afStarND;

    gConst = nonDimParams.gConst;
    BTT = false;


    X0 = [r0; v0; m0];
    odeoptions = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
     [tTraj, stateTraj] = ode45(@(t, X) trajectory(t, X, gamma, kr, ...
                          tgo0, isp, rMoonND, rfStar, vfStar, afStar,refVals.T_ref, refVals, gConst, nonDimParams.minThrustND, nonDimParams.maxThrustND, delta_t, BTT), [0, tgo0], X0, odeoptions);

     rState = stateTraj(:,1:3);
     vState = stateTraj(:,4:6);
     mState = stateTraj(:,7);
     tgoState = tgo0 - tTraj;
    
     flag_thrustGotLimited = false;
     

     numStates = size(rState,1);
     aTList = zeros(3,numStates);
     aTNormList = zeros(1,numStates);
    


     for idx = 1:numStates
         r = rState(idx, :).';
         v = vState(idx, :).';
         m = mState(idx);
         tgo = tgoState(idx);
         %g = -(rMoonND^2) * r / (norm(r)^3);
         gGuidance = gConst;

         minAccel = nonDimParams.minThrustND/m;
         maxAccel = nonDimParams.maxThrustND/m;
         % BTT
         
         if BTT
             if (tgo * refVals.T_ref) > (1.5*delta_t)
                 [rfVirtual, vfVirtual, afVirtual, tgoVirt] = computeBeyondTerminationTargeting(r, v, gamma, kr, rfStar, vfStar, afStar, delta_t, tgo, gGuidance);
                 aT1 = gamma*(kr/(2*gamma +4) -1)*afVirtual;
                 aT2 = (gamma*kr/(2*gamma+4)-gamma-1)*gGuidance;
                 aT3 = ((gamma+1)/tgoVirt)*(1-kr/(gamma+2))*(vfVirtual-v);
                 aT4 = (kr/tgoVirt^2)*(rfVirtual-r-v*tgoVirt);
                 aTi = aT1 + aT2 + aT3 + aT4;
                 % if norm(aTi) > maxAccel
                 %     aTi = aTi / norm(aTi) * maxAccel;
                 %     flag_thrustGotLimited = true;
                 % elseif norm(aTi) < minAccel
                 %     aTi = aTi / norm(aTi) * minAccel;
                 %     flag_thrustGotLimited = true;
                 % end
             else
                 aT1 = gamma*(kr/(2*gamma +4) -1)*afVirtual;
                 aT2 = (gamma*kr/(2*gamma+4)-gamma-1)*gGuidance;
                 aT3 = ((gamma+1)/tgoVirt)*(1-kr/(gamma+2))*(vfVirtual-v);
                 aT4 = (kr/tgoVirt^2)*(rfVirtual-r-v*tgoVirt);
                 aTi = aT1 + aT2 + aT3 + aT4;
                 % if norm(aTi) > maxAccel
                 %     aTi = aTi / norm(aTi) * maxAccel;
                 %     flag_thrustGotLimited = true;
                 % elseif norm(aTi) < minAccel
                 %     aTi = aTi / norm(aTi) * minAccel;
                 %     flag_thrustGotLimited = true;
                 % end
             end
         else    % No BTT
             if (tgo * refVals.T_ref) > 0.2
                 aT1 = gamma*(kr/(2*gamma +4) -1)*afStar;
                 aT2 = (gamma*kr/(2*gamma+4)-gamma-1)*gGuidance;
                 aT3 = ((gamma+1)/tgo)*(1-kr/(gamma+2))*(vfStar-v);
                 aT4 = (kr/tgo^2)*(rfStar-r-v*tgo);
                 aTi = aT1 + aT2 + aT3 + aT4;
                 if norm(aTi) > maxAccel
                     aTi = aTi / norm(aTi) * maxAccel;
                     flag_thrustGotLimited = true;
                 elseif norm(aTi) < minAccel
                     aTi = aTi / norm(aTi) * minAccel;
                     flag_thrustGotLimited = true;
                 end
             else
                 aTi = aTList(:,idx-1);
             end

         end
         aTList(:,idx) = aTi;
         aTNormList(idx) = norm(aTi);
     end
end


%% Functions
function dXdt = trajectory(t, X, gamma, kr, tgo0, isp, rMoonND, rfStar, vfStar, afStar, T_ref, refVals, gConst, minThrust, maxThrust, delta_t, BTT)
    r    = X(1:3);
    v    = X(4:6);
    mass = X(7);
    tgo  = tgo0 - t;

    g = -(rMoonND^2) * r / (norm(r)^3);
    gGuidance = gConst;

    minAccel = minThrust/mass;
    maxAccel = maxThrust/mass;

    persistent aT;
    persistent rfVirtual
    persistent vfVirtual
    persistent afVirtual
    persistent tgoVirt
    if BTT
        if (tgo * T_ref) > (1.5*delta_t)
            [rfVirtual, vfVirtual, afVirtual, tgoVirt] = computeBeyondTerminationTargeting(r, v, gamma, kr, rfStar, vfStar, afStar, delta_t, tgo, gGuidance);
            aT1 = gamma*(kr/(2*gamma +4) -1)*afVirtual;
            aT2 = (gamma*kr/(2*gamma+4)-gamma-1)*gGuidance;
            aT3 = ((gamma+1)/tgoVirt)*(1-kr/(gamma+2))*(vfVirtual-v);
            aT4 = (kr/tgoVirt^2)*(rfVirtual-r-v*tgoVirt);
            aT = aT1 + aT2 + aT3 + aT4;
            if norm(aT) > maxAccel
                aT = aT / norm(aT) * maxAccel;
            elseif norm(aT) < minAccel
                aT = aT / norm(aT) * minAccel;
            end
            
        else

            aT1 = gamma*(kr/(2*gamma +4) -1)*afVirtual;
            aT2 = (gamma*kr/(2*gamma+4)-gamma-1)*gGuidance;
            aT3 = ((gamma+1)/tgoVirt)*(1-kr/(gamma+2))*(vfVirtual-v);
            aT4 = (kr/tgoVirt^2)*(rfVirtual-r-v*tgoVirt);  
            aT = aT1 + aT2 + aT3 + aT4;
            if norm(aT) > maxAccel
                aT = aT / norm(aT) * maxAccel;
            elseif norm(aT) < minAccel
                aT = aT / norm(aT) * minAccel;
            end
        end
    else
        if (tgo * T_ref) > 0.2
            aT1 = gamma*(kr/(2*gamma +4) -1)*afStar;
            aT2 = (gamma*kr/(2*gamma+4)-gamma-1)*gGuidance;
            aT3 = ((gamma+1)/tgo)*(1-kr/(gamma+2))*(vfStar-v);
            aT4 = (kr/tgo^2)*(rfStar-r-v*tgo);
            aT = aT1 + aT2 + aT3 + aT4;
            if norm(aT) > maxAccel
                aT = aT / norm(aT) * maxAccel;
            elseif norm(aT) < minAccel
                aT = aT / norm(aT) * minAccel;
            end
        else
            aT=aT;
        end

    end
    F_mag = norm(aT) * mass;
    dm_dt = -F_mag / (isp);
    dXdt = [v; aT + g; dm_dt];
end