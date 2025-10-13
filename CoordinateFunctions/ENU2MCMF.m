function mcmf = ENU2MCMF(enu, anchorLatDeg, anchorLonDeg, isPosition)
    Rm = 173.601428; % ND Lunar Radius

    lat0 = deg2rad(anchorLatDeg);
    lon0 = deg2rad(anchorLonDeg);

    [E0, N0, U0] = enuBasis(lat0, lon0);
    R = [E0, N0, U0]';

    if isPosition
        r0 = Rm * U0;
        mcmf = R.' * enu + r0;
    else
        mcmf = R.' * enu;
    end
end
