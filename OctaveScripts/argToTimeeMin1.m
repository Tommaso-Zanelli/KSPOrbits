function [t] = argToTimeeMin1(th, e)
%
% function [t] = argToTimeeMin1(th, e)
%
% Contains the analytical expression to compute the time given the argument
% for an elliptical orbit

    % Analytical expression
	t = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(0.5 * th)) / sqrt((1 - (e ^ 2)) ^ 3);
	t = t - e * sin(th)./ ((1 - (e ^ 2)) * (1 + e * cos(th)));

end
