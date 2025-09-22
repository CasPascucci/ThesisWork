function [localPosition, localVelocity] = PDIToLocalCartesian(altitude, longitudeDeg, latitudeDeg, landingLonDeg, landingLatDeg, inertialVel, flightPathAngleDeg, lunarRadius)



if nargin < 8
    lunarRadius = 1737.4;
end

% Convert deg to rad
longitudeRad = deg2rad(longitudeDeg);
latitudeRad = deg2rad(latitudeDeg);
landingLongRad = deg2rad(landingLonDeg);
landingLatRad = deg2rad(landingLatDeg);
flightPathAngleRad = deg2rad(flightPathAngleDeg);

deltaLon = longitudeRad - landingLongRad;
east = lunarRadius*deltaLon * cos(landingLatRad);

deltaLat = latitudeRad - landingLatRad;
north = lunarRadius*deltaLat;

up = altitude;
localPosition = [east; north; up];
localPosition = localPosition * 1000; % km to m

if nargout > 1
    % unitRadial = [cos(latitudeRad)*cos(longitudeRad), cos(latitudeRad)*sin(longitudeRad), sin(latitudeRad)];
    % unitEast = [-sin(longitudeRad), cos(longitudeRad), 0];
    % unitNorth = [-sin(latitudeRad)*cos(longitudeRad), -sin(latitudeRad)*sin(longitudeRad), cos(latitudeRad)];
    % 
    % velRadial = inertialVel * sin(flightPathAngleRad);
    % velEast = 0;
    % velNorth = -inertialVel * cos(flightPathAngleRad);


    %While flight path angle is zero and the PDI is directly north of the
    %landing site, then velocity can be entire set due South
    
    
    localVelocity = [0; -inertialVel; 0]; %already in m
end
