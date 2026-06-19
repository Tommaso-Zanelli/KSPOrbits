function [of] = applyManoeuvre(oi, m, mu)
%
%   function [of] = applyManoeuvre(oi, m, mu)
%
%   Applies a manoeuvre to an orbit and returns the subsequent orbit
%

    % Gets position and velocity before the manoeuve
    [xx1, vv1, ~] = structTimeToPosVelArg(m.t, oi, mu);

    % Gets the local frame of reference
    T = refFrameLocal(xx1, vv1);

    % Velocity jump vector in the absolute frame of reference
    dv = T * m.v;

    % Gets the absolute velocity after the manoeuvre
    vv2 = vv1 + dv;

    % Gets the orbital parameters after the manoeuvre
    [~, of] = posVelToStruct(xx1, vv2, mu, m.t, 1);

end