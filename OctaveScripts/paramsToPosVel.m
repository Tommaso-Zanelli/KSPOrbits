function [x, v] = paramsToPosVel(argument, p, e, lan, incl, aop, mu)
%
% function [x, v] = paramsToPosVel(argument, p, e, lan, incl, aop, mu)
%
% Computes the position and (if requested) the velocity at a given (or for
% a given set of) argument(s), for the given orbital parameters

	% Gets the sizes of the input array in order to reshape it
	nrows = size(argument, 1);
	ncolumns = size(argument, 2);
	n_el = nrows * ncolumns;

    argument = reshape(argument, 1, n_el);

    % Computes the distance
    rho = p./ (1 + e * cos(argument));

    % Computes the directional unit vectors in the local reference frame
    uvect_r =  [ cos(argument); sin(argument); zeros(1, n_el)];
    uvect_th = [-sin(argument); cos(argument); zeros(1, n_el)];

    % Computes the rotational matrix for the local reference frame
    R = rOrb(lan, incl, aop);

    % Computes the position in the local frame of reference
    x = repmat(rho, 3, 1).* uvect_r;

    % Rotates the position
    x = R * x;

    % If the velocity was requeste
    if nargout > 1

        % Computes the velocity in the local frame of reference
        v = sqrt(mu / p) * (repmat(e * sin(argument), 3, 1).* uvect_r + repmat((1 + e * cos(argument)), 3, 1).* uvect_th);

        % Rotates the velocity
        v = R * v;

    end

end