function [th_1, th_2, p, e, lan, incl, aop] = optiTransfer2(x1, x2, costf, mu, dx, nPoints, d1Tol, nIterMax)
%
%   function [th_1, th_2, p, e, lan, incl, aop] = optiTransfer2(x1, x2, costf, mu, delta_th, n_s, nm_nt, ep_ntS)
%
%   Given two points around a body, and a cost function, finds the optimal
%   transfer orbit. The cost function must be in the form costf(vt1, vt2)
%

    % Maximum number of iterations, if unspecified
    if nargin < 8
        nIterMax = 1000;
    end

    % First derivative tolerance, if unspecified
    if nargin < 7
        d1Tol = 1e-12;
    end

    % Number of section points, if unspecified
    if nargin < 6
        nPoints = 12;
    end

    % Variable increment, if unspecified
    if nargin < 5
        dx = 4e-6;
    end

    if norm(cross(x1, x2), 2) < eps(min(abs([x1; x2])))

        % Re-defines the cost function for a given inclination
        cfi = @(incl0) costIncl(x1, x2, costf, mu, incl0, dx, nPoints, d1Tol, nIterMax);

        % Defines the search interval
        incl_inc = 2 * pi / (nPoints - 3);
        incli = [-pi - incl_inc, pi + incl_inc];

        % Finds the optimum inclination
        [incl0, ~, ~] = optiFun(incli, cf, dx, nPoints, d1Tol, nIterMax);

        % Finds the corresponding optimum th0
        th0 = findOptIncl(x1, x2, costf, mu, incl0, dx, nPoints, d1Tol, nIterMax);

    else

        % Performs the search for the optimum
        incl0 = 0;
        th0 = findOptIncl(x1, x2, costf, mu, incl0, dx, nPoints, d1Tol, nIterMax);

    end

    % Finds he orbital parameters
    [th_1, th_2, p, e, lan, incl, aop] = transferOrbit(x1, x2, th0, incl0);

end

% Performs the search for the optimum
function th0 = findOptIncl(x1, x2, costf, mu, incl0, dx, nPoints, d1Tol, nIterMax)

        % Re-defines the const function as a function of the sole th0
        cf = @(th0) transferOrbitCost(x1, x2, costf, mu, th0, incl0);

        % Defines the search interval
        th_inc = 2 * pi / (nPoints - 3);
        thi = [-pi - th_inc, pi + th_inc];

        % Finds the optimum
        [th0, ~, ~] = optiFun(thi, cf, dx, nPoints, d1Tol, nIterMax);
end

% Cost given an inclination
function c = costIncl(x1, x2, costf, mu, incl0, dx, nPoints, d1Tol, nIterMax)

    % Finds the optimum for the given inclination
    th0 = findOptIncl(x1, x2, costf, mu, incl0, dx, nPoints, d1Tol, nIterMax);

    % Computes the corresponding cost
    c =  transferOrbitCost(x1, x2, costf, mu, th0, incl0);

end

