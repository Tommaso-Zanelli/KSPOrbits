function [th_1, th_2, p, e, lan, incl, aop] = transferOrbit(x1, x2, th0, incl0)
%
%   function [th_1, th_2, p, e, lan, inc, aop] = transferOrbit(x1, x2, th0, incl0)
%
%   Plots a transfer orbit between two points.
%

    % If the last argument was not specified
    if nargin < 4
        incl0 = 0;
    end

    % sine and cosine of the angle between x1 and x2
    sindth = norm(cross(x1, x2), 2);
    cosdth = dot(x1, x2);

    % Distances
    r1 = norm(x1, 2);
    r2 = norm(x2, 2);

    % Angle between x1 and x2
    dth = atan2(sindth, cosdth);

    % If the two points are aligned
    if sindth < eps(min(abs([x1; x2])))

        % Initializes a frame of reference
        i1 = x1 / r1;
        [~, d3] = min(abs(i1));
        i3 = [0; 0; 0];
        i3(d3(1, 1), 1) = 1;
        i3 = i3 - dot(i1, i3) * i1;
        i3 = i3 / norm(i3, 2);
        i2 = cross(i3, i1);

        % Gets the second direction based on the specified inclination
        v1 = i2 * cos(incl0) + i3 * sin(incl0);

    % If the two points are not aligned
    else

        % Gets the second direction as the difference between the points
        v1 = x2 - x1;

    end

    % Angles identifying the plane of the two points
    [lan, incl, aop_th] = rotAngles(x1, v1);

    % Distance ratio
    xi = r1 / r2;

    % determinant
    d = (cos(dth + th0) - xi * cos(th0));

    % Semi-latus rectum
    p = (r2 * xi * (cos(dth + th0) - cos(th0))) / d;

    % Eccentricity
    e = (xi - 1) / d;

    % Returns error values for degenerate cases
    if (d == 0) || (p < 0) || (e < 0)
        th_1 = NaN;
        th_2 = NaN;
        p    = NaN;
        e    = NaN;
        lan  = NaN;
        incl = NaN;
        aop  = NaN;
        return;
    end

    % Initializes the argument of periapsis
    aop = (aop_th - th0);

    % Checks whether the argument of periapsis is correct
    [x1_ck, ~] = paramsToPosVel(th0, p, e, lan, incl, aop, 0);
    d_aop = acos(max([min([(dot(x1, x1_ck) / (norm(x1, 2) * norm(x1_ck, 2))), 1]), -1]));
    aop = inMPiPiInt(aop - d_aop);

    % Assigns the initial and final angle
    th_1 = inMPiPiInt(th0);
    th_2 = th_1 + dth;

    % Excludes infinite transfers
    if e >= 1

        % Angle limits
        th_l1 = acos(1 / e);
        th_l2 = -th_l1;

        % If any of the points are in the "non-existent" portion of the orbit
        if (th_1 < th_l1) || (th_1 > th_l2) || (th_2 < th_l1) || (th_2 > th_l2)
            th_1 = NaN;
            th_2 = NaN;
            p    = NaN;
            e    = NaN;
            lan  = NaN;
            incl = NaN;
            aop  = NaN;
            return;
        end

        % if
        if (inMPiPiInt(th_2) < th_1) || (th_2 < inMPiPiInt(th_1))
            th_1 = NaN;
            th_2 = NaN;
            p    = NaN;
            e    = NaN;
            lan  = NaN;
            incl = NaN;
            aop  = NaN;
            return;
        end

    end

end