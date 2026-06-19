function [th1, th2, otp, dv1, dv2, er] = getTransferOrbitStatic(prm, o1p, o2p, mu)
%
% function [th1, th2, otp] = getTransferOrbitStatic(prm, o1p, o2p, mu)
%
% Designs a transfer orbit between the orbit of parameters contained in o1p
% and the orbit of parameters contained in o2p. prm contains the values for
% the three possible degrees of freedom
%

    % Gets the orbital parameters in a more readable form for the first orbit
    if size(o1p, 2) > 1
        o1p = o1p';
    end
    [p1, e1, lan1, incl1, aop1] = vectorToParams(o1p);

    % Gets the orbital parameters in a more readable form for the second orbit
    if size(o2p, 2) > 1
        o2p = o2p';
    end
    [p2, e2, lan2, incl2, aop2] = vectorToParams(o2p);

    % Gets the degrees of freedom in a more readable form
    if size(prm, 2) > 1
        prm = prm';
    end
    th1_1 = prm(1, 1);          % The anomaly of the starting point on the first orbit
    th2_2 = prm(2, 1);          % the anomaly of the final point on the second orbit
    th0_3  = prm(3, 1);         % The "anomaly offset" for the transfer orbit

    % Computes the starting point
    if nargout > 3
        [x1, v1] = paramsToPosVel(th1_1, p1, e1, lan1, incl1, aop1, mu);
    else
        [x1, ~] = paramsToPosVel(th1_1, p1, e1, lan1, incl1, aop1, mu);
    end

    % Computes the arrival point
    if nargout > 4
        [x2, v2] = paramsToPosVel(th2_2, p2, e2, lan2, incl2, aop2, mu);
    else
        [x2, ~] = paramsToPosVel(th2_2, p2, e2, lan2, incl2, aop2, mu);
    end

    % Radiuses of the two points
    rho1 = norm(x1, 2);
    rho2 = norm(x2, 2);

    % Cross product between the two, and its norm
    cr  = cross(x1, x2);
    ncr = norm(cr, 2);

    % Vector normal to the transfer orbit plane
    k = cr / ncr;

    % Angle between the two points
    Dth = atan2(ncr, dot(x1, x2));

    % First and last angles of the transfer orbit
    th1 = th0_3;
    th2 = Dth + th0_3;
    if th2 > pi
        th2 = th2 - 2 * pi;
    end
    if th2 < -pi
        th2 = th2 + 2 * pi;
    end

    % Gets p and the eccentricity
    detM = rho2*cos(th2) - rho1*cos(th1);
    p3 = rho1 * rho2 * (cos(th2) - cos(th1)) / detM;
    e3 = (rho1 - rho2) / detM;

    % Aborts if e is less than zero
    if (p3 < 0 || e3 < 0)
        er = Inf;
        dv1 = [Inf; Inf; Inf];
        dv2 = [Inf; Inf; Inf];
        [otp] = paramsToVector(p3, e3, 0, 0, 0);
        return
    end

    % Directional unit vectors
    ir1 = x1 / norm(x1, 2);
    ir2 = x2 / norm(x2, 2);
    ith1 = cross(k, ir1);
    ith1 = ith1 / norm(ith1, 2);
    ith2 = cross(k, ir2);
    ith2 = ith2 / norm(ith2, 2);

    % Gets the velocities at the transfer points in the local references
    if nargout > 5
        [rho1c, v3_1lr] = paramsToPosVel(th1, p3, e3, 0, 0, -th1, mu);
        [rho2c, v3_2lr] = paramsToPosVel(th2, p3, e3, 0, 0, -th2, mu);
    else
        [~, v3_1lr] = paramsToPosVel(th1, p3, e3, 0, 0, -th1, mu);
        [~, v3_2lr] = paramsToPosVel(th2, p3, e3, 0, 0, -th2, mu);
    end

    % Gets these velocities in the absolute frame of teference
    v3_1 = v3_1lr(1, 1) * ir1 + v3_1lr(2, 1) * ith1;
    v3_2 = v3_2lr(1, 1) * ir2 + v3_2lr(2, 1) * ith2;

    % gets the remaining parameters
    if nargout > 5
        [th1c, p3c1, e3c1, lan3, incl3, aop3, er1] = posVelToParams1p(x1, v3_1, mu);
        [th2c, p3c2, e3c2, lan3c, incl3c, aop3c, er2] = posVelToParams1p(x2, v3_2, mu);
    else
        [~, ~, ~, lan3, incl3, aop3, ~] = posVelToParams1p(x1, v3_1, mu);
    end

    % Computes the velocity differences if requested
    if nargout > 3
        dv1 = v3_1 - v1;
    end
    if nargout > 4
        dv2 = v2 - v3_2;
    end

    % Performs error checks if requested
    if nargout > 5
       er = ([ ...
       %er = max([ ...
       er1; ...
       er2; ...
       abs(lan3 - lan3c) / (2 * pi); ...
       abs(incl3 - incl3c) / (2 * pi); ...
       abs(aop3 - aop3c) / (2 * pi); ...
       abs(e3c1 - e3) / max([e3, 1]); ...
       abs(e3c2 - e3) / max([e3, 1]); ...
       abs(p3c1 - p3) / p3; ...
       abs(p3c2 - p3) / p3; ...
       abs(th1c - th1) / (2 * pi); ...
       abs(th2c - th2) / (2 * pi); ...
       abs(v3_1lr(3, 1) / sqrt(sum((v3_1lr(1:2, 1).^ 2), 1))); ...
       abs(v3_2lr(3, 1) / sqrt(sum((v3_2lr(1:2, 1).^ 2), 1))); ...
       abs(sqrt(sum((rho1c(2:3, 1).^ 2), 1)) / rho1c(1, 1)); ...
       abs(sqrt(sum((rho2c(2:3, 1).^ 2), 1)) / rho2c(1, 1)); ...
       abs((rho1c(1, 1) - rho1) / rho1); ...
       abs((rho2c(1, 1) - rho2) / rho2)]);
    end

    % Compresses the output orbital parameters
    [otp] = paramsToVector(p3, e3, lan3, incl3, aop3);

end