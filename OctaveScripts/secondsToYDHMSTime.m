function [ydhms, seconds] = secondsToYDHMSTime(seconds, toLower)
%
%   function [ydhms, seconds] = secondsToYDHMSTime(seconds, toLower)
%
%   Converts a time given in seconds to (Kerbin) years, days, hours,
%   minutes and seconds, each rounded to the smaller integer, except for
%   the seconds. These are rounded to the smaller integer if toLower is
%   set to true or unspecified, rounded to the closest integer otherwise
%   The remaining fractions of seconds are returned separately
%

    % If toLower was not specified, it is assumed true
    if nargin < 2
        toLower = true;
    end

    % years, days, hours, minutes and seconds... in seconds
    q     = [9203400, 21600, 3600, 60, 1];

    % Initializes the output
    ydhms = zeros(1, 5);

    % For each element in q
    for idx_1 = 1:5
        ydhms(1, idx_1) = floor(seconds / q(1, idx_1));
        seconds = seconds - q(1, idx_1) * ydhms(1, idx_1);
    end

    % Corrects the seconds in case the closest integer is requested
    if (~toLower && (seconds > 0.5))
        q(1, 5) = q(1, 5) + 1;
        if nargout > 1
            seconds = seconds - 1;
        end
    end

end