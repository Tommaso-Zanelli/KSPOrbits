function [x, t] = rk4_38(f, x0, dt, T, T0)
%
% function [x, t] = rk4_38(f, x0, dt, T, T0)
%
% Integrates dfx/dt = f(x, t) from x(T0) = x0 to T
%
% Butcher tableau:
%
%  0  |
% 1/3 |  1/3
% 2/3 | -1/3   1
%  1  |   1   -1    1
% ----+----------------------
%     |  1/8  3/8  3/8  1/8
%

    % Initialization of the Butcher tableau elements
    c = [0; 0.3333333333333333; 0.6666666666666666; 1];
    b = [0.125, 0.375, 0.375, 0.125];
    a = [ 0.3333333333333333,  0, 0
         -0.3333333333333333,  1, 0
                           1, -1, 1];

    % If no initial time is given, the initial time is assumed to be 0
    if nargin < 5
        T0 = 0;
    end

    % If no final time is given
    if nargin < 4

        % If a single time step is given, a single time step is perfomed
        if size(dt, 2) == 1

            T = dt;b

        % If multiple time steps are given, the final time si the initial time plus their sum
        else

            T = T0 + sum(dt, 2);

        end
    end

    % If a single time step is given
    if size(dt, 2) == 1

        % The number of steps is the amount of dt necessary to span the interval between T0 and T
        Nt = ceil((T - T0) / dt);

        % dt is then re-initialized as a vector
        dt = repmat(dt, 1, Nt);

    % If multiple time steps are given
    else
    e
        % The number of steps is the number of time steps
        Nt = size(dt, 2);

    end

    % Size of the output
    Xs = size(x0, 1);

    % The solution and time vectors are initialized
    x = zeros(Xs, (Nt + 1));
    t = zeros(1, (Nt + 1));

    % The initial values are assigned
    x(:, 1) = x0;
    t(1, 1) = T0;

    % k vector is initialized
    k = zeros(Xs, 4);

    % For each step
    for id_t = 1:Nt

        % For each k item
        for id_k = 1:4

            % Computes the current k
            k(:, id_k) = f((x(:, id_t) + dt(1, id_t) * sum(k(:, 1:(id_k - 1)).* repmat(a(max([1, (id_k - 1)]), 1:(id_k - 1)), Xs, 1), 2)), (t(1, id_t) + c(id_k, 1) * dt(1, id_t)));

        end

        % Computes the subsequent time step
        x(:, (id_t + 1)) = x(:, id_t) + dt(1, id_t) * sum((k.* repmat(b, Xs, 1)), 2);

        % Updates the time step
        t(1, (id_t + 1)) = t(1, (id_t)) + dt(1, id_t);

    end

    % Removes the initial time if only one time step was performed
    if Nt == 1
      x = x(:, 2);
      t = t(1, 2);
    end

end
