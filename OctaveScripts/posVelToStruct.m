function [th, o, t2] = posVelToStruct(x, v, mu, t, mn, att)
%
% [th, o] = posVelToStruct(x, v, mu, t, mn)
%
% Computes the orbital parameters given a position, velocity, time and
% gravitational constant of the attractor. Can take multiple positions,
% velocities and times. If no time is specified, periapsis time is assigned
% as null. Parameter mn determines behaviour when multiple points are
% given:
%       - 1:    Gets all parameters from the first point, uses the others to
%               compute the anomaly
%       - 2:    Computes all parameters from each point, averages all except
%               the anomaly, which is re-computed
%       - 3:    Computes all parameters from each point, averages all except
%               the anomaly, which is passed as is
%

    % if time was not provided
    if nargin < 4
        t = [];     % Default is no time
    end

    % Default mn parameters. Also option 4 from "posVelToParams" is not available
    if nargin < 5
        mn = 1;             % 1: Gets all parameters from the first point, uses the others to compute the anomaly
                            % 2: Computes all parameters from each point, averages all except the anomaly, which is re-computed
                            % 3: Computes all parameters from each point, averages all except the anomaly, which is passed as is
    else
        if mn > 3
            mn = 1;
        end
    end

    % Additional arguments for "argToTime"
    if nargin < 6
        att = [];
    end
    att = check_att(att);

    % Gets all the parameters
    [th, o.p, o.e, o.lan, o.incl, o.aop] = posVelToParams(x, v, mu, mn);

    % Computes the output time if requested
    if numel(t) > 0 || nargout > 2

         % Computes the time for (t(pe) = 0)
        t_pe0 = argToTime(th, o.p, o.e, mu, 0, att.ep, att.d);

    end

    % If time is to be computed
    if numel(t) > 0

        % Computes the periapsis time
        if mn == 1
            o.peT = t(1, 1) - t_pe0(1, 1);
        else
            o.peT = mean(t - t_pe0);
        end

    else

        % Default value for the periapsis time is zero
        o.peT = 0;

    end

    % Time recomputed if requested
    if nargout > 2

        if size(t, 2) == size(x, 2)
            t2 = t;
        else
            t2 = t_pe0 + o.peT;
        end

    end

end