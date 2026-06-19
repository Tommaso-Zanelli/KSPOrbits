function [x] = optiGrad(fun, x, delta, dx, nN)
%
%   function [x] = optiGrad(fun, x, dx, nN)
%
%   Performs optimization on function fun using the gradient method,
%   starting from x.
%   "dx" contains initial values for the "delta" values of
%   nN contains the number of iterations after which a test is performed on
%   the Hessian of the function. If the function is found to be convex at
%   that point, all further iterations perform Newton's method. if nN is
%   null, no check of the Hessian is performed and the ordinary gradient
%   method is performed until convergence.
%

    % If nN was not specified
    if nargin < 5
        nN = 0;
    end

    % If the initial vector is not a column vector
    transposed = false;
    if size(x, 2) > 1
        transposed = true;
        x = x';
    end

    % Number of variables
    N = size(x, 1);

    % If no initial dx was specified
    if nargin < 4 || (numel(dx) == 0)
        dx = 2.^ (ceil(log2(abs(x))) - 18);          % Assigns a default value
        dx(dx == 0) = (eps(x(dx == 0))).^ 0.4;
    end

    % If less than the required number of dx-es was provided
    if numel(dx) < numel(x)
        dx = ones(N, 1) * dx(1, 1);      % Assigns the first (or only) value to all
    end

    % Initializes the iteration switches
    iterate = true;
    newton = false;

    % Computes the cost of the reference values
    cr = fun(x);

    % Initializes the iteration counter
    c_int = 0;

    % Main cycle
    while iterate

        % Computes the gradient
        %tic
        g = gradientFun(fun, x, dx);
        %toc
        % Stores the current value
        x_old = x;

        % for the classic gradient method
        if ~newton

            % Normalizes the gradient
            g = g / norm(g, 2);

            %  --- Looks for the minimum along the gradient ---------------

            % Whether the
            doubler = true;

            % Initializes the current minimum and the scalar increment
            cc = cr;
            ddx = 1;

            % Computes the increment
            q = -g * ddx;

            % Initializes the new cost
            cl = Inf;

            % Tests the increment
            while (cl > cc) %max(abs(q)) > min(abs(dx))

                % Computes the increment
                q = -g * ddx;

                % Computes the cost at the incremented position
                cl = fun(x + q);

                % If the cost is greater than the current minimum
                if cl > cc

                    % Halves the scalar increment and disables doubling
                    ddx = 0.5 * ddx;
                    doubler = false;

                % If a new minimum was found
                else

                    % Increases the current position
                    x = x + q;

                    % Re-assigns the current minimum
                    cc = cl;

                    % doubles the increment if requested
                    if doubler
                        ddx = 2 * ddx;
                    end

                end

            end

            % -------------------------------------------------------------

            % Increases the iteration counter
            c_int = c_int + 1;

        % for the Newton method
        else

            % Performs Newton step
            inz = find(g ~= 0);
            x = x - H \ g;
            c_int = c_int + 1;
            cc = fun(x);

        end

        % Checks the current decrement
        %fprintf('%g\n%g\n\n', abs(cr - cc), max([delta * abs(cr), eps(cr)]))
        fprintf('%g\n', cc)
        if (abs(cr - cc) < max([delta * abs(cr), eps(cr)])) || (cr < cc)
            iterate = false;
        else
            if cr < cc
                fprintf('Error!\n')
                return;
            end
            cr = cc;
        end

        % --- Re-evaluates dx ---------------------------------------------

        % Increment for the current step
        x_inc = x - x_old;

        % minimum acceptable values
        dx_min = 2.^ (ceil(log2(abs(x))) - 24);
        dx_min(dx_min == 0) = (eps(x(dx_min == 0))).^ 0.4;

        % Re-assigns each value
        for idx_1 = 1:N
            dx(idx_1, 1) = max([min([x_inc(idx_1, 1), dx(idx_1, 1)]), dx_min(idx_1, 1)]);
        end

        % -----------------------------------------------------------------

        % Checks the Hessian if requested
        if (nN && ~mod(c_int, nN)) || newton
            H = HessianFun(fun, x, dx);
            h = eig(H);
            newton = all(h > eps(max(max(abs(H)))));
        end

    end

    % Transposes the result if a transposition had taken place
    if transposed
        x = x';
    end

end