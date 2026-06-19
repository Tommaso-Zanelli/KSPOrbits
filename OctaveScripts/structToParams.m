function [p, e, lan, incl, aop, peT] = structToParams(o)
%
%   function [p, e, lan, incl, aop, peT] = structToParams(o)
%
%   Splits up an orbit data structure into the individual parameters.
%

    % Assigns semi-latus rectum and eccentricity
    p   = o.p;
    e   = o.e;

    % Assigns each orientation angle requested
    if nargout > 2
        lan  = o.lan;
    end
    if nargout > 3
        incl = o.incl;
    end
    if nargout > 4
        aop  = o.aop;
    end

    % Assigns the periapsis time if requested
    peT = o.peT;

end