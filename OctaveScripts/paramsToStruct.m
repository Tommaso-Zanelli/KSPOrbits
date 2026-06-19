function [o] = paramsToStruct(p, e, lan, incl, aop, peT)
%
%   function [o] = paramsToStruct(p, e, lan, incl, aop, peT)
%
%   Wraps up all orbital parameters into a single structure.
%

    % If the periapsis time is not specified
    if nargin < 6
        peT  = 0;
    end

    % If the orientation angles are not specified
    if nargin < 5
        aop  = 0;
    end
    if nargin < 4
        incl = 0;
    end
    if nargin < 3
        lan  = 0;
    end

    % Fills the structure
    o.p    = p;
    o.e    = e;
    o.lan  = lan;
    o.incl = incl;
    o.aop  = aop;
    o.peT  = peT;

end