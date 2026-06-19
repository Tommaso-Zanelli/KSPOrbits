function [o] = join4DPoints(x1, t1, x2, t2, mu, d1, ep, d, N)
%
%   function [o] = join4DPoints(x1, t1, x2, t2, mu, d1, ep, d)
%
%   Given two points at two instants, around an attractor of gravitational
%   parameter mu, finds the transfer orbit.
%

	if nargin < 6
        d1 = x2;
	end
	if nargin < 7
        ep = 1e-4;                      % The size of the interval arount pi
                                        % which is to be linearized
    end
	if nargin < 8
        d = 16;                         % The negative esponent of the power of
                                        % two representing the interval around
                                        % e = 1 which is to be linearized
    end
	if nargin < 9
        N = 101;                        % The number of points initially used
                                        % to find the initial guess
    end

    % First point radius
	rho1 = norm(x1, 2);

    % Second point radius
    rho2 = norm(x2, 2);

    % Angular distance between the two points
    dth = atan2(norm(cross(x1, x2), 2), dot(x1, x2));

    % Corrects for preferred direction
    if dot(cross(x1, x2), cross(x1, d1)) < 0
        dth = 2 * pi - dth;
    end

    % Defines function to be "zeroed"
    fz = @(th0) th02PlanarTransferTime(rho1, rho2, th0, dth, mu, ep, d) ...
                - t2 + t1;

    % Finds the best "initial guess" for zeroing the function



end