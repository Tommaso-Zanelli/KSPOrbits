function [t] = argToTimeEllipse(th, e, ep)
%
% function [t] = argToTimeEllipse(th, e, ep)
%
% computes the time given the argument and the eccentricityfor an elliptic
% orbit, linearizing the interval around pi. It also corrects for multiple
% full rounds.
%

	if nargin < 3
		ep = 1e-4;                      % The size of the interval arount pi
                                        % which is to be linearized
    end

    % Initializes the return vector
	t   =  zeros(size(th));

    % Computes the orbital period (halved)
	T   =  pi / sqrt((1 - (e ^ 2)) ^ 3);

    % Reduces the argument to the interval [-pi, pi]
	n   =  floor((th + pi) / (2 * pi));
	th2 =  th - 2 * pi * n;

	m1 = th2 < (ep - pi);                           % indexes of arguments too close to -pi
	m2 = (th2 >= (ep - pi)) & (th2 <= (pi - ep));   % indexes of arguments surely between -pi and pi
	m3 = th2 > (pi - ep);                           % indexes of arguments too close to pi

	t(m1) = -T + (T + argToTimeeMin1((ep - pi), e)) * ((th2(m1) + pi) / ep);    % Linearized expression
	t(m2) = argToTimeeMin1(th2(m2), e);                                         % Analytical expression
	t(m3) = T + (T - argToTimeeMin1((pi - ep), e)) * ((th2(m3) - pi) / ep);     % Linearized expression

    % Corrects the period
	T = 2 * T;

    % Adds the period to the computed time
	t = t + T * n;

end
