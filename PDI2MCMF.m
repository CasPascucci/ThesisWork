function [r_mcmf, v_mcmf] = PDI2MCMF(altitude_km, lonInitDeg, latInitDeg, ...
                                      landingLonDeg, landingLatDeg, ...
                                      inertialVel_mps, flightPathAngleDeg, ...
                                      lunarRadius_km)
%PDIToMCMF Convert PDI inputs to Moon-Centered, Moon-Fixed position and velocity.
% Outputs: r_mcmf [m], v_mcmf [m/s] in MCMF.
%
% Assumptions:
% - Initial position is at (latInit, lonInit) at altitude above a spherical Moon.
% - Initial velocity is directed along the great-circle toward the landing site,
%   with flight-path angle measured positive up from local horizontal.
%
% Inputs:
%   altitude_km           scalar
%   lonInitDeg, latInitDeg, landingLonDeg, landingLatDeg   degrees
%   inertialVel_mps       scalar speed at PDI (m/s)
%   flightPathAngleDeg    degrees (positive up)
%   lunarRadius_km        scalar radius (default 1737.4 km)

    if nargin < 8 || isempty(lunarRadius_km)
        lunarRadius_km = 1737.4;
    end

    % Units
    radMoon = lunarRadius_km * 1000;        % m
    h  = altitude_km  * 1000;          % m
    latLEM = deg2rad(latInitDeg);  lonLEM = deg2rad(lonInitDeg);
    latLand = deg2rad(landingLatDeg); lonLand = deg2rad(landingLonDeg);
    flightPathAngle  = deg2rad(flightPathAngleDeg);
    
    [E1, N1, U1] = enuBasis(latLEM, lonLEM); % Defines what East North and Up are at the initial site

    r_mcmf = (radMoon +h) * U1;
    if nargout > 1
        az = getHeading(latLEM, lonLEM, latLand, lonLand);
        u_h = sin(az)*E1 + cos(az)*N1;
        v_mcmf = inertialVel_mps * (cos(flightPathAngle)*u_h + sin(flightPathAngle)*U1);
    end
end

%% Helpers

function [E, N, U] = enuBasis(lat, lon)
    clat = cos(lat); slat = sin(lat);
    clon = cos(lon); slon = sin(lon);
    U = [clat*clon; clat*slon; slat];
    E = [-slon; clon; 0];
    N = [-slat*clon; -slat*slon; clat];
end

function az = getHeading(lat1, lon1, lat2, lon2)
    deltaLon = lon2 - lon1;
    y = sin(deltaLon) * cos(lat2);
    x = cos(lat1)*sin(lat2) - sin(lat1)*cos(lat2)*cos(deltaLon);
    az = atan2(y, x); % from North, east-positive
end