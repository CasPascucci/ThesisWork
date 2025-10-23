function [r_mcmf, v_mcmf] = PDI2MCMF(altitude_km, lonInitDeg, latInitDeg, ...
                                      landingLonDeg, landingLatDeg, ...
                                      inertialVel_mps, flightPathAngleDeg, azimuth)

    lunarRadius_km = 1737.4; % km
    radMoon = lunarRadius_km * 1000; % m
    h  = altitude_km  * 1000; % m
    latPDI = deg2rad(latInitDeg);  lonPDI = deg2rad(lonInitDeg);
    latLand = deg2rad(landingLatDeg); lonLand = deg2rad(landingLonDeg);
    flightPathAngle  = deg2rad(flightPathAngleDeg);
    
    [E1, N1, U1] = enuBasis(latPDI, lonPDI); % Defines what East North and Up are at the initial site

    r_mcmf = (radMoon +h) * U1;
    if nargout > 1
        %az = getHeading(latPDI, lonPDI, latLand, lonLand);
        az = azimuth;
        u_h = sin(az)*E1 + cos(az)*N1;
        v_mcmf = inertialVel_mps * (cos(flightPathAngle)*u_h + sin(flightPathAngle)*U1);
    end
end

%% Functions

function az = getHeading(lat1, lon1, lat2, lon2)
    deltaLon = lon2 - lon1;
    y = sin(deltaLon) * cos(lat2);
    x = cos(lat1)*sin(lat2) - sin(lat1)*cos(lat2)*cos(deltaLon);
    az = atan2(y, x); % from North, east-positive
end