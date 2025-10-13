function [enu, alt] = MCMF2ENU(X, landingLatDeg, landingLonDeg, isPosition, isDim)
if nargin < 5 || ~isDim
    Rm = 173.601428; % ND Lunar Radius
else
    Rm = 173.601428 * 10000; % Dim Lunar Radius
end



lat0 = deg2rad(landingLatDeg);
lon0 = deg2rad(landingLonDeg);

[E0, N0, U0] = enuBasis(lat0,lon0);
R = [E0, N0, U0]';

if isPosition
    Xorig = X;
    r0 = Rm * U0; % position of landing in MCMF
    X = X - r0;
    if nargout > 1
        alt = vecnorm(Xorig, 2, 1) - Rm;
    end
end

enu = R * X;
end