function [dt] = th02PlanarTransferTime(rho1, rho2, th0, dth, mu, ep, d)
%
%   function [dt] = th02PlanarTransferTime(rho1, rho2, th0, dth, mu, ep, d)
%
%   Given two radii and their angular distance (identifying two points),
%   an initial argument and a gravitational parameter, computes the time
%   between them on a transfer orbit.
%

	if nargin < 6
        ep = 1e-4;                   	% The size of the interval arount pi
                                        % which is to be linearized
    end
	if nargin < 7
        d = 16;                         % The negative esponent of the power of
                                        % two representing the interval around
                                        % e = 1 which is to be linearized
    end

	% Determinant
	dd = rho2 * cos(th0 + dth) - rho1 * cos(th0);

	% Inverse of the determinant
	dm1 = 1 / dd;

	% Semi-latus rectum
	p = rho1 * rho2 * (cos(th0 + dth) - cos(th0)) * dm1;

	% Eccentricity
	e = (rho1 - rho2) * dm1;

	% Times of the two points
	[t1] = argToTime(th0, p, e, mu, 0, ep, d);
	[t2] = argToTime((th0 + dth), p, e, mu, 0, ep, d);

	% Time interval
	dt = t2 - t1;

end