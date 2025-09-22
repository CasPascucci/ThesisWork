function [enu, alt] = MCMF2ENU(X, landingLatDeg, landingLonDeg, isVector, lunarRad_km)

if nargin < 5 || isempty(lunarRad_km)
    lunarRad_km = 1737.4; % Default value for lunar radius in kilometers
end

Rm = lunarRad_km * 1000;

lat0 = deg2rad(landingLatDeg);
lon0 = deg2rad(landingLonDeg);

[E0, N0, U0] = enuBasis(lat0,lon0);

if ~isVector
    r0 = Rm*U0; % position of landing in MCMF
    X = X - r0;
end

e = E0' * X;
n = N0' * X;
u = U0' * X;
enu = [e; n; u];

if nargout > 1 % Might be unnecessary, but this is only useful for radius calcs
    alt = vecnorm(X, 2, 1) - Rm;
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