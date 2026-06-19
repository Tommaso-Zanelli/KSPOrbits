function [c] = transferOrbitCost(x1, x2, cf, mu, th0, incl0)
%
%   function [c] = transferOrbitCost(x1, x2, cf, mu_K, th0, incl0)
%
%   Gets the cost of a transfer orbit given a cost function in the form
%   cf(v1, v2)
%

    % If the last argument was not specified
    if nargin < 6
        incl0 = 0;
    end

    % Plots the transfer orbit
    [th_1, th_2, p, e, lan, incl, aop] = transferOrbit(x1, x2, th0, incl0);

    % Infinite for degenerate cases
    if isnan(p)
        c = Inf;
        return;
    end

    % Velocities at the manoeuvre nodes
    [~, v1] = paramsToPosVel(th_1, p, e, lan, incl, aop, mu);
    [~, v2] = paramsToPosVel(th_2, p, e, lan, incl, aop, mu);

    % Comptues the cost
    c = cf(v1, v2);

end