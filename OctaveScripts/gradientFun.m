function [g] = gradientFun(fun, x, dx)
%
%   function [g] = gradientFun(fun, x, dx)
%
%   Computes the gradient of a multi-variable function
%

    % If no dx was provided
    if nargin < 3
        dx = 2.^ (ceil(log2(abs(x))) - 18);     % Assigns a default value
        dx(dx == 0) = eps(x(dx == 0)).^ 0.4;
    end

    % If less than the required number of dx-es was provided
    if numel(dx) < numel(x)
        dx = ones(size(x)) * dx(1, 1);      % Assigns the first (or only) value to all
    end

    % Initializes the gradient
    g = zeros(size(x));

    % For each value
    for idx_1 = 1:size(x, 1)
        for idx_2 = 1:size(x, 2)

            cxp1 = x;
            cxp1(idx_1, idx_2) = cxp1(idx_1, idx_2) + dx(idx_1, idx_2);     % Value increased along the current direction

            cxm1 = x;
            cxm1(idx_1, idx_2) = cxm1(idx_1, idx_2) - dx(idx_1, idx_2);     % Value decreased along the current direction

            % Computes the current component of the gradient
            g(idx_1, idx_2) = (fun(cxp1) - fun(cxm1)) / (2 * dx(idx_1, idx_2));

        end
    end

end 