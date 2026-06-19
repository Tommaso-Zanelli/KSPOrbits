function [t] = argToTimeHyperbola(th, e)
%
% function [t] = argToTimeHyperbola(th, e)
%
% Computes the time given the argumet and the eccentricity for an
% hyperbolic orbit, assigning Inf or NaN to invalid values
%

    % Initializes the return vector
    t = zeros(size(th));

    % Computes the limit argument(s)
    th_l = acos(((-1) / e));

    m1 = (th < -(th_l + eps(th_l)));                                % Indexes of arguments below the lower limit
    m2 = (th >= -(th_l + eps(th_l)) & th <= -(th_l - eps(th_l)));   % Indexes of arguments at the lower limit
    m3 = (th > -(th_l - eps(th_l)) & th < (th_l - eps(th_l)));      % Indexes of arguments in the valid interval
    m4 = (th >= (th_l - eps(th_l)) & th <= (th_l + eps(th_l)));     % Indexes of arguments at the upper limit
    m5 = (th > (th_l + eps(th_l)));                                 % Indexes of arguments above the lower limit

    t(m1) = NaN;                        % Assigns NaN to values below the lower limit
    t(m2) = -Inf;                       % Assigns -Inf to values at the lower limit
    t(m3) = argToTimeeMax1(th(m3), e);  % Computes the times for arguments in the valid interval
    t(m4) = Inf;                        % Assigns Inf to values at the upper limit
    t(m5) = NaN;                        % Assigns NaN to values above the upper limit

end
