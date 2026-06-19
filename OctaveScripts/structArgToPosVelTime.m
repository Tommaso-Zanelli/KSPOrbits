function [x, v, t] = structArgToPosVelTime(argument, o, mu, att)
%
% function [x, v, t] = structArgToPosVelTime(argument, o, mu)
%
% Computes the position and (if requested) the velocity at a given (or for
% a given set of) argument(s), for the given orbital parameters
%

    % Additional arguments for "argToTime"
    if nargin < 4
    	att = [];
    end
    att = check_att(att);

    % Handles variable input / output
    if nargin < 3
        if nargout > 1
            x = NaN;
            v = NaN;
            return;
        else
            mu = 0;
        end
    end

    % Computes position and velocity
    if nargout > 1
        [x, v] = paramsToPosVel(argument, o.p, o.e, o.lan, o.incl, o.aop, mu);
    else
        [x, ~] = paramsToPosVel(argument, o.p, o.e, o.lan, o.incl, o.aop, mu);
    end

    % Computes the time if requested
    if nargout > 2
        t = argToTime(argument, o.p, o.e, mu, o.peT, att.ep, att.d);
    end

end