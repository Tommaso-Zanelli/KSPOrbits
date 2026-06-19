function [t] = argToTimeParabola(th)
%
% function [t] = argToTimeParabola(th)
%
% computes the time given the argument for a parabolic orbit

    % Analytic expression
    t = (sin(th).* (2 + cos(th)))./ (3 * ((1 + cos(th)).^ 2));

end