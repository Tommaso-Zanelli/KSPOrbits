function [x, v, argument] = structTimeToPosVelArg(t, o, mu, tta)
%
% function [x, v, argument] = structTimeToPosVelArg(t, o, mu, tta)
%
% Computes the position and (if requested) the velocity at a given (or for
% a given set of) argument(s), for the given orbital parameters
%

    % Additional arguments for "timeToArg" if unspecified
	if nargin < 4
        tta = [];
    end
    tta = check_tta(tta, o.e);

    % Gets the argument
    argument = timeToArg(t, o.p, o.e, mu, o.peT, ...
        tta.e_bs, tta.e_nt, tta.n_bs, tta.n_nt, tta.ep, tta.d);

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

end