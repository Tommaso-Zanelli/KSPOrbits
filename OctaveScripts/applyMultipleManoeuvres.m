function [os, ths, dvs] = applyMultipleManoeuvres(o0, t, dv, phi1, phi2, mu)
%
%   function [os, ths] = applyMultipleManoeuvres(o0, t, dv, phi1, phi2)
%
%   aaa
%

    % Gets the number of manoeuvres
    Nm = size(t, 2);

    % Initializes the delta-vs
    if nargout > 2
        dvs = zeros(3, Nm);
    end

    % Assigns the first element of the output structure
    os(1, 1) = o0;

    % Initializes the container for the arguments
    ths = zeros(2, (Nm + 1));

    % For each manoeuvre
    for idx_1 = 1:Nm

        % Gets argument
        ths(2, idx_1) = timeToArg(t(1, idx_1), os(1, idx_1).p, ...
            os(1, idx_1).e, mu, os(1, idx_1).peT);

        % Gets position and velocity before the manoeuvre
        [xx1, vv1] = paramsToPosVel(ths(2, idx_1), os(1, idx_1).p, ...
            os(1, idx_1).e, os(1, idx_1).lan, os(1, idx_1).incl, ...
            os(1, idx_1).aop, mu);

        % Gets the reference frame at the current manoeuvre point
        [T] = refFrameLocal(xx1, vv1);

        % Builds the velocity variation
        v_vect = T(:, 1) * cos(phi1(1, idx_1)) * cos(phi2(1, idx_1)) + ...
                 T(:, 3) * sin(phi1(1, idx_1)) * cos(phi2(1, idx_1)) + ...
                 T(:, 2) * sin(phi2(1, idx_1));
        c_dv   = dv(1, idx_1) * v_vect;

        % Assigns the delta-vs
        dvs(:, idx_1) = dv(1, idx_1) * [cos(phi1(1, idx_1)) * cos(phi2(1, idx_1)); ...
            sin(phi2(1, idx_1)); sin(phi1(1, idx_1)) * cos(phi2(1, idx_1))];

        % Gets the new velocity
        vv2 = vv1 + c_dv;

        % Gets the new orbit
        [ths(1, (idx_1 + 1)), os(1, (idx_1 + 1)).p, os(1, (idx_1 + 1)).e, ...
         os(1, (idx_1 + 1)).lan, os(1, (idx_1 + 1)).incl, ...
         os(1, (idx_1 + 1)).aop] = posVelToParams(xx1, vv2, mu);

        % Gets the base time for the current position
        Tb = argToTime(ths(1, (idx_1 + 1)), os(1, (idx_1 + 1)).p, ...
                       os(1, (idx_1 + 1)).e, mu, 0);

        % Gets the periapsis time of the new orbit
        os(1, (idx_1 + 1)).peT = t(1, idx_1) - Tb;

    end

end