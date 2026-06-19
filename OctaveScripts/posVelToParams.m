function [th, p, e, lan, incl, aop, er] = posVelToParams(x, v, mu, mn)
%
% function [th, p, e, lan, incl, aop] = posVelToParams(x, v, mu, mn)
%
% The

    if nargin < 4
        mn = 1;             % 1: Gets all parameters from the first point, uses the others to compute the anomaly
                            % 2: Computes all parameters from each point, averages all except the anomaly, which is re-computed
                            % 3: Computes all parameters from each point, averages all except the anomaly, which is passed as is
                            % 4: Computes all parameters and anomalies, returning them as vectors
    end

    % Checks wether the input is oriented in the right way
    if size(x, 2) == 3 && size(x, 1) ~= 3
        x = x';
    end
    if size(v, 2) == 3 && size(v, 1) ~= 3
        v = v';
    end

    % Checks wether the input has the correct size, returns errors otherwise
    if size(x, 2) ~= size(v, 2) || size(x, 1) ~= 3
        fprintf('Runtime error: function posVelToParams received input with inconsistent sizes.\n')
        th = NaN;
        p = NaN;
        e = NaN;
        lan    = NaN;
        incl   = NaN;
        aop = NaN;
        return;
    end

    % Sets the number of points
    np = size(x, 2);

    % Gets all parameters from the first point, uses the others to compute the anomaly
    if mn == 1

        % Gets all parameters except for the anomaly
        [~, p, e, lan, incl, aop, ~] = posVelToParams1p(x(:, 1), v(:, 1), mu);

        % Initializes the output container
        th     = zeros(1, np);
        if nargout > 6
            er   = zeros(1, np);
        end

        % Computes the corresponding rotation matrix
        R = rOrb(lan, incl, aop);

        % Rotates the position
        xp = (R') * x;

        % Computes the angle
        th = atan2(xp(2, :), xp(1, :));

        % Computes the error
        if nargout > 6

            rho = (p./ (1 + e * cos(th)));

            er = max([abs(xp(3, :))./ sqrt(sum((xp(1:2, :).^ 2), 1)); ...
                      abs(vp(3, :))./ sqrt(sum((vp(1:2, :).^ 2), 1)); ...
                      abs((rho - (sqrt(sum((xp(1:2, :).^ 2), 1))))./ rho); ...
                      abs((p_v - p) / p); abs(e_v - e); ...
                      abs((lan_v - lan) / (2 * pi)); abs((incl_v - incl) / (2 * pi)); ...
                      abs((aop_v - aop) / (2 * pi))]);

        end

    % Computes all parameters from each point, averages all except the anomaly, which is re-computed
    elseif mn == 2

        % Initializes the output containers
        th     = zeros(1, np);
        p_v    = zeros(1, np);
        e_v    = zeros(1, np);
        lan_v  = zeros(1, np);
        incl_v = zeros(1, np);
        aop_v  = zeros(1, np);
        if nargout > 6
            er   = zeros(1, np);
        end

        % For each point
        for idx_1 = 1:np

            % Computes all parameters
            [~, p_v(1, idx_1), e_v(1, idx_1), lan_v(1, idx_1), incl_v(1, idx_1), aop_v(1, idx_1), ~] = posVelToParams1p(x(:, idx_1), v(:, idx_1), mu);

        end

        % "Normalizes" the rotational angles
        [lan_v, incl_v, aop_v] = normalizeRotAngles(lan_v, incl_v, aop_v);

        % Averages all parameters except for the anomaly
        p    = mean(p_v);
        e    = mean(e_v);
        lan  = mean(lan_v);
        incl = mean(incl_v);
        aop  = mean(aop_v);

        % Computes the corresponding rotation matrix
        R = rOrb(lan, incl, aop);

        % Rotates the position
        xp = (R') * x;

        % Computes the angle
        th = atan2(xp(2, :), xp(1, :));

        % Computes the error
        if nargout > 6

            rho = (p./ (1 + e * cos(th)));

            er = max([abs(xp(3, :))./ sqrt(sum((xp(1:2, :).^ 2), 1)); ...
                      abs(vp(3, :))./ sqrt(sum((vp(1:2, :).^ 2), 1)); ...
                      abs((rho - (sqrt(sum((xp(1:2, :).^ 2), 1))))./ rho); ...
                      abs((p_v - p) / p); abs(e_v - e); ...
                      abs((lan_v - lan) / (2 * pi)); abs((incl_v - incl) / (2 * pi)); ...
                      abs((aop_v - aop) / (2 * pi))]);

        end

    % Computes all parameters from each point, averages all except the anomaly, which is passed as is
    elseif mn == 3

        % Initializes the output containers
        th     = zeros(1, np);
        p_v    = zeros(1, np);
        e_v    = zeros(1, np);
        lan_v  = zeros(1, np);
        incl_v = zeros(1, np);
        aop_v  = zeros(1, np);
        if nargout > 6
            er   = zeros(1, np);
        end

        % For each point
        for idx_1 = 1:np

            % Computes all parameters
            if nargout > 6
                [th(1, idx_1), p_v(1, idx_1), e_v(1, idx_1), lan_v(1, idx_1), ...
                 incl_v(1, idx_1), aop_v(1, idx_1), er(1, idx_1)] = ...
                    posVelToParams1p(x(:, idx_1), v(:, idx_1), mu);
            else
                [th(1, idx_1), p_v(1, idx_1), e_v(1, idx_1), lan_v(1, idx_1), ...
                 incl_v(1, idx_1), aop_v(1, idx_1), ~] = ...
                    posVelToParams1p(x(:, idx_1), v(:, idx_1), mu);
            end

        end

        % "Normalizes" the rotational angles
        [lan_v, incl_v, aop_v] = normalizeRotAngles(lan_v, incl_v, aop_v);

        % Averages all parameters except for the anomaly
        p    = mean(p_v);
        e    = mean(e_v);
        lan  = mean(lan_v);
        incl = mean(incl_v);
        aop  = mean(aop_v);

        % Error estimates
        if nargout > 6

            er = max([er; abs((p_v - p) / p); abs(e_v - e); ...
                      abs((lan_v - lan) / (2 * pi)); abs((incl_v - incl) / (2 * pi)); ...
                      abs((aop_v - aop) / (2 * pi))]);

        end

    % Computes all parameters and anomalies, returning them as vectors
    elseif mn == 4

        % Initializes the output containers
        th   = zeros(1, np);
        p    = zeros(1, np);
        e    = zeros(1, np);
        lan  = zeros(1, np);
        incl = zeros(1, np);
        aop  = zeros(1, np);
        if nargout > 6
            er   = zeros(1, np);
        end

        % For each point
        for idx_1 = 1:np

            % Computes all parameters
            if nargout > 6
                [th(1, idx_1), p(1, idx_1), e(1, idx_1), lan(1, idx_1), incl(1, idx_1), aop(1, idx_1), er(1, idx_1)] = posVelToParams1p(x(:, idx_1), v(:, idx_1), mu);
            else
                [th(1, idx_1), p(1, idx_1), e(1, idx_1), lan(1, idx_1), incl(1, idx_1), aop(1, idx_1), ~] = posVelToParams1p(x(:, idx_1), v(:, idx_1), mu);
            end

        end

        % "Normalizes" the rotational angles
        [lan, incl, aop] = normalizeRotAngles(lan, incl, aop);

    % For any other number, returns an error message and aborts
    else
        fprintf('Runtime error: invalid option.\n')
        th = NaN;
        p = NaN;
        e = NaN;
        lan    = NaN;
        incl   = NaN;
        aop = NaN;
        return;

    end

end
