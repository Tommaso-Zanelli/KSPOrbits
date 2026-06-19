function [x0, xs, xn] = optiFun(xInterval, costFun, dx, nPoints, d1Tol, nIterMax)
%
%   function [x0, xs, xn] = optiFun(xInterval, costFun, dx, nPoints, d1Tol, nIterMax)
%
%

    % Maximum number of iterations, if unspecified
    if nargin < 6
        nIterMax = 1000;
    end

    % First derivative tolerance, if unspecified
    if nargin < 5
        d1Tol = 1e-12;
    end

    % Number of section points, if unspecified
    if nargin < 4
        nPoints = 12;
    end

    % Variable increment, if unspecified
    if nargin < 3
        dx = 4e-6;
    end

    % Initializes first and second order derivatives
    d1 = Inf;
    d2 = -Inf;

    % Initializes the number of iterations
    nIter = 0;

    % Initializes the container for the cost
    cc = zeros(1, nPoints);

    % If a second output was required
    if nargout > 1
        xs = [];
    end

    % Initializes the final guess
    x0 = 1e308;

    % Divides the interval until a minimum in a convex region is found
    while (abs(d1) > d1Tol) && (nIter < nIterMax) && (d2 <= 0)

        % Search points
        xx = linspace(xInterval(1, 1), xInterval(1, 2), nPoints);

        % Values of the cost function
        for idx_1 = 1:nPoints
            cc(1, idx_1) = costFun(xx(1, idx_1));
        end

        %cc'

        % Finds the minimum
        [~, ii] = min(cc);

        % Finds the difference
        Dx = xx(1, ii) - x0;

        % Updates the current minimum
        x0 = xx(1, ii);

        % Updates the dx used for the finite differences
        dx = min([dx, max([0.125 * abs(Dx), (eps(x0) * (2 ^ 26))])]);

        % If a second output was required
        if nargout > 1
            xs = [xs; x0];
        end

        % Finds the first and second order derivatives at the minimum
        [d1, d2] = d12(x0, dx, costFun);

        % Redefines the search interval for the subsequent iteration
        xInterval = [xx(1, max([1, (ii - 1)])), xx(1, min([nPoints, (ii + 1)]))];

        % Increases the iteration counter
        nIter = nIter + 1;

    end

    % If a third output was required
    if nargout > 2
        xn = [];
    end

    % Performs a certain number of Newton iterations to refine the minimum
    while (abs(d1) > d1Tol) && (nIter < nIterMax)

        % Finds the difference
        Dx = - 0.5 * dx * d1 / d2;

        % Newton update
        x0 = x0 + Dx;

        % Updates the dx used for the finite differences
        dx = min([dx, max([0.125 * abs(Dx), (eps(x0) * (2 ^ 26))])]);

        % If a third output was required
        if nargout > 2
            xn = [xn; x0];
        end

        % Finds the first and second order derivatives at the minimum
        [d1, d2] = d12(x0, dx, costFun);

        % Increases the iteration counter
        nIter = nIter + 1;

    end

end

% Finite differences for the first and second order derivatives
function [d1, d2] = d12(x, dx, cf)

    c = cf(x);
    cp1 = cf(x + dx);
    cm1 = cf(x - dx);
    d1 = cp1 - cm1;
    d2 = cp1 + cm1 - 2 * c;

end

