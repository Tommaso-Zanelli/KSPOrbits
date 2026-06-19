function [th, p, e, lan, incl, aop, er] = posVelToParams1p(x, v, mu)
%
% function [th, p, e, lan, incl, aop] = posVelToParams1p(x, v, mu)
%
% Given a point, the velocity at that point, and the gravitational
% constant, computes all orbital parameters, anomaly included

    % If the input does not consist of three element vectors, aborts execution
    if numel(x) > 3 || numel(v) > 3
        fprintf('Runtime error: function posVelToParams1p received wrongly shaped input\n')
    end

    % If the input is not a colum vector, transposes it
    if size(x, 1) == 1 && size(x, 2) == 3
        x = x';
    end
    if size(v, 1) == 1 && size(v, 2) == 3
        v = v';
    end

    % Gets the angles identifying the plane of the orbit
    [lan, incl, aop_th] = rotAngles(x, v);

    % Computes the corresponding rotation matrix
    R = rOrb(lan, incl, aop_th);

    % Rotates the position and velocity
    xp = (R') * x;
    vp = (R') * v;

    % Computes the error
    if nargout > 6
        er = max([(max(abs(xp(2:3, 1))) / abs(xp(1, 1))), ...
                (abs(vp(3, 1)) / sqrt(sum((vp(1:2, 1).^ 2), 1)))]);
    end

    % Computes p, the eccentricity and the anomaly
    ect = ((xp(1, 1) * (vp(2, 1) ^ 2)) / mu) - 1;
    est = (xp(1, 1) * vp(1, 1) * vp(2, 1)) / mu;
    th = atan2(est, ect);
    e = sqrt((est ^ 2) + (ect ^ 2));
    p = xp(1, 1) * (1 + e * cos(th));

    % Gets the argument of the periapsis
    aop = aop_th - th;
    if aop < -pi
        aop = aop + 2 * pi;
    end
    if aop > pi
        aop = aop - 2 * pi;
    end

end