function [o, ec] = rawOrbitToStruct(ro, R)
%
%   function [o] = rawOrbitToStruct(ro)
%
%   Converts a raw orbit data to an orbit data structure.
%

    % If the radius is not specified
    if nargin < 2
        R = 0;
    end

    % Apoapsis and periapsis to semi-latus rectum and eccentricity
    [o.p, ec] = apoPeriToPe((ro(1, 1) + R), (ro(1, 2) + R));

    % Eccentricity
    if (size(ro, 2) >= 3)
        o.e = ro(1, 3);
    else
        o.e = ec;
    end

    % Longitude of ascending node
    o.lan = 0;
    if (size(ro, 2) >= 4)
        o.lan = ro(1, 4) * 1.745329251994329701382e-2;
        while (o.lan > pi)
            o.lan = o.lan - (2 * pi);
        end
        while (o.lan < -pi)
            o.lan = o.lan + (2 * pi);
        end
    end

    % Inclination
    o.incl = 0;
    if (size(ro, 2) >= 5)
        o.incl = ro(1, 5) * 1.745329251994329701382e-2;
        while (o.incl > pi)
            o.incl = o.incl - (2 * pi);
        end
        while (o.incl < -pi)
            o.incl = o.incl + (2 * pi);
        end
    end

    % Argument of periapsis
    o.aop = 0;
    if (size(ro, 2) >= 6)
        o.aop = ro(1, 6) * 1.745329251994329701382e-2;
        while (o.aop > pi)
            o.aop = o.aop - (2 * pi);
        end
        while (o.aop < -pi)
            o.aop = o.aop + (2 * pi);
        end
    end

    % Periapsis time
    o.peT = 0;

end