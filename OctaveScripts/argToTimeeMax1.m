function [t] = argToTimeeMax1(th, e)
%
% function [t] = argToTimeeMax1(th, e)
%
% Contains the  analytical expression to compute the time given the
% argument for an hyperbolic orbit

    % Analytical expression
	t = -2 * atanh(sqrt((e - 1) / (1 + e)) * tan(0.5 * th)) / sqrt(((e ^ 2) - 1) ^ 3);
	t = t + e * sin(th)./ (((e ^ 2) - 1) * (1 + e * cos(th)));

end
