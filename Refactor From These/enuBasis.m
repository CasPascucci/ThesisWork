function [E, N, U] = enuBasis(lat, lon)
    clat = cos(lat); slat = sin(lat);
    clon = cos(lon); slon = sin(lon);

    U = [clat*clon;   clat*slon; slat];
    E = [-slon;       clon;      0];
    N = [-slat*clon; -slat*slon; clat];
end