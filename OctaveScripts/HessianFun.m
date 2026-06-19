function [H] = HessianFun(fun, x, dx)
%
%   function [H] = HessianFun(fun, x, dx)
%
%   Computes the Hessian of a multi-variable function
%

    % If x is not a column vector
    transposed = false;
    if size(x, 2) > 1
        x = x';             % Transposes it
        transposed = true;
    end

    % Gets the number of elements in x
    N = size(x, 1);

    % If no dx was provided
    if nargin < 3
        dx = 2.^ (ceil(log2(abs(x))) - 18);     % Assigns a default value
        dx(dx == 0) = (eps(x(dx == 0))).^ 0.4;
    end

    % If less than the required number of dx-es was provided
    if numel(dx) < numel(x)
        dx = ones(N, 1) * dx(1, 1);      % Assigns the first (or only) value to all
    end

    % Initializes the Hessian
    H = zeros(N);

    % Cycles through all possible combinations of values
    for idx_1 = 1:N

        % Computes the diagonal value
        cxp1 = x;
        cxp1(idx_1, 1) = cxp1(idx_1, 1) + dx(idx_1, 1);     % Value increased along the current direction
        cxm1 = x;
        cxm1(idx_1, 1) = cxm1(idx_1, 1) - dx(idx_1, 1);     % Value decreased along the current direction
        % Computes the second order derivative
        H(idx_1, idx_1) = (fun(cxp1) + fun(cxm1) - 2 * fun(x)) / ...
                          (dx(idx_1, 1) ^ 2);

        % Computes all extra-diagonal values
        for idx_2 = (idx_1 + 1):N

            cxp1p2 = x;
            cxp1p2(idx_1, 1) = cxp1p2(idx_1, 1) + dx(idx_1, 1);
            cxp1p2(idx_2, 1) = cxp1p2(idx_2, 1) + dx(idx_2, 1);

            cxm1p2 = x;
            cxm1p2(idx_1, 1) = cxm1p2(idx_1, 1) - dx(idx_1, 1);
            cxm1p2(idx_2, 1) = cxm1p2(idx_2, 1) + dx(idx_2, 1);

            cxp1m2 = x;
            cxp1m2(idx_1, 1) = cxp1m2(idx_1, 1) + dx(idx_1, 1);
            cxp1m2(idx_2, 1) = cxp1m2(idx_2, 1) - dx(idx_2, 1);

            cxm1m2 = x;
            cxm1m2(idx_1, 1) = cxm1m2(idx_1, 1) - dx(idx_1, 1);
            cxm1m2(idx_2, 1) = cxm1m2(idx_2, 1) - dx(idx_2, 1);

            % Computes the second order derivative
            H(idx_1, idx_2) = (fun(cxp1p2) - fun(cxp1m2) ...
                             - fun(cxm1p2) + fun(cxm1m2)) / ...
                             (4 * dx(idx_1, 1) * dx(idx_2, 1));
            H(idx_2, idx_1) = H(idx_1, idx_2);

        end

    end

    % Transposes the results
    if transposed
        H = H';
    end

end