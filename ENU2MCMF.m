function r_or_v_mcmf = ENU2MCMF(enu, anchorLatDeg, anchorLonDeg, isVector)
% ENU_to_MCMF_general
% Convert ENU to MCMF for either positions or free vectors.
%
% Usage:
%   % Vector mode (velocity, acceleration, thrust)
%   v_mcmf = ENU_to_MCMF_general(v_enu, latDeg, lonDeg, 'isVector', true);
%
%   % Position mode from surface origin
%   r_mcmf = ENU_to_MCMF_general(enu_pos, latDeg, lonDeg, 'isVector', false, ...
%                                'lunarRadius_km', 1737.4);
%
% Inputs:
%   enu              3xN or Nx3 matrix of ENU components [E; N; U] in meters
%   anchorLatDeg     anchor latitude in degrees
%   anchorLonDeg     anchor longitude in degrees
%
% Output:
%   r_or_v_mcmf      3xN MCMF vectors (positions if isVector=false, vectors if true)



    lunarRadius_km = 1737.4;

    % Ensure ENU is 3xN
    E = enu;
    if size(E,1) ~= 3
        E = E.';  % allow Nx3
        if size(E,1) ~= 3
            error('enu must be 3xN or Nx3');
        end
    end
    

    % Basis at anchor
    lat0 = deg2rad(anchorLatDeg);
    lon0 = deg2rad(anchorLonDeg);
    [E0, N0, U0] = enuBasis(lat0, lon0);

    % Rotate ENU components to MCMF
    v_mcmf = E0*E(1,:) + N0*E(2,:) + U0*E(3,:);

    if isVector
        % Free vector: no origin shift
        r_or_v_mcmf = v_mcmf;
        return;
    end

    % Position mode: add anchor origin on the reference sphere and altitude
    Rm = lunarRadius_km * 1000.0;
    r0 = Rm * U0;  % anchor on the surface

    r_or_v_mcmf = r0 + v_mcmf;
   
end

% ============ Helpers ============
function [E, N, U] = enuBasis(lat, lon)
    clat = cos(lat); slat = sin(lat);
    clon = cos(lon); slon = sin(lon);
    U = [clat*clon; clat*slon; slat];
    E = [-slon; clon; 0];
    N = [-slat*clon; -slat*slon; clat];
end
