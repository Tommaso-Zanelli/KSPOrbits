function [o] = rawOrbitToStructOrbit(rCells, R)
%
%   function [o] = rawOrbitToStructOrbit(rCells, R)
%
%   Given an orbit reported as:
%   [pe      ap      e       lan     incl	aop
%    pe_s	pe_y	pe_d	pe_hh	pe_mm	pe_ss
%    now_s	now_y	now_d	now_hh	now_mm	now_ss]
%   Derives the corresponding structure
%

    % If no radius was specified
    if nargin < 2
        R = 0;
    end

    % Semi-latus rectum and eccentricity
    [o.p, o.e]	=	apoPeriToPe(rCells(1, 1) + R, rCells(1, 2) + R);

    % Orbital plane angles
    o.lan       =   rCells(1, 4) * 1.7453292519943295076197955e-2;
    o.incl      =   rCells(1, 5) * 1.7453292519943295076197955e-2;
    o.aop       =   rCells(1, 6) * 1.7453292519943295076197955e-2;

    % Time to periapsis
    ttp = rCells(2, 1) * ydhmsTimeToSeconds(rCells(2, 2:6));

    % Now
    now = rCells(3, 1) * ydhmsTimeToSeconds(rCells(3, 2:6));

    % Periapsis time
    o.peT = ttp + now;

end