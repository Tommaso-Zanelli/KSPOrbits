function [c] = preciseOrbitCost(orb0, ref0, times0, refvs, mnvrs, mu)
%
%   function [c] = preciseOrbitCost(orb0, ref0, times0, refvs, mnvrs, mu)
%
%   Given a guess for the actual orbital parameters of an orbit, compares
%   it with the reported reference values for it, then applies different
%   maneouvres and compares the resulting parameters with reported
%   reference values.
%
%   Comparison is done through the cost function cf0p5, which assigns a
%   negligible cost for values under 0.5, and an increasing exponentially
%   one after 0.5
%

    % Number of manoeuvres
    N = size(refvs, 1);

    % Gets cost of main orbit
    c = orbitCompareCost(ref0, orb0);

    % Gets times cost
    c    = c + cf0p5((times0.pe - orb0.peT) / times0.peR);	% Confronts periapsis time
    T    = argToTime(pi, orb0.p, orb0.e, mu, orb0.peT);     % Computes half the period
    c    = c + cf0p5((times0.T - (2 * T)) / times0.TR);     % Confronts period
    Tap1 = orb0.peT + T;                                    % Computes subsequent apoapsis time
    Tap2 = orb0.peT - T;                                    % Computes previous apoapsis time
    if abs(times0.pe - Tap1) < abs(times0.pe - Tap2)        % Selects more appropriate apoapsis time
        Tap = Tap1;
    else
        Tap = Tap2;
    end
    c    = c + cf0p5((times0.ap - Tap) / times0.apR);       % Confronts apoapsis time

    % For each manoeuvre
    for idx_1 = 1:N

        % Applies current manoeuvre
        acm = timeToArg(mnvrs(idx_1, 1).t, orb0.p, orb0.e, mu, orb0.peT);       % Computes the argument of the manoeuvre point
        [xx, vv] = paramsToPosVel(acm, orb0.p, orb0.e, orb0.lan, ...            % Computes position and velocity at the manoeuvre point
                                  orb0.incl, orb0.aop, mu);
        [T] = refFrameLocal(xx, vv);                                            % Computes the local frame of reference at the manoeuvre point
        dv = T(:, 1) * mnvrs(idx_1, 1).v1 + T(:, 2) * mnvrs(idx_1, 1).v2 ...    % Computes the manoeuvre delta-v in the "absolute" frame of reference
           + T(:, 3) * mnvrs(idx_1, 1).v3;
        vvam = vv + dv;                                                         % Computes the velocity after the manoeuvre
        [~, p_am, e_am, lan_am, incl_am, aop_am] = ...                          % Computes the parameters of the orbit after the manoeuvre
            posVelToParams(xx, vvam, mu);
        o_am = paramsToStruct(p_am, e_am, lan_am, incl_am, aop_am, 0);          % Packs the paramenters in a data structure

        % Coomputes and adds cost
        c = c + orbitCompareCost(refvs(idx_1, 1), o_am);

    end

end