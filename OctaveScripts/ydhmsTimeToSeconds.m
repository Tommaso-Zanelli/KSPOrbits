function [seconds, seconds_precision] = ydhmsTimeToSeconds(ydhms)
%
%   function [seconds, seconds_precision] = ydhmsTimeToSeconds(ydhms)
%
%   Converts a time given in (Kerbin) years, days, hours, minutes and
%   seconds to seconds, interpreting NaN values as beyiond precision,
%   and computing a confidence value.
%

    % Adds missing fields
    while size(ydhms, 2) < 5

        ydhms = [0, ydhms];

    end

    % years, days, hours, minutes and seconds... in seconds
    q = [9203400, 21600, 3600, 60, 1];

    % Time in secods
    seconds = dot(ydhms(1, ~isnan(ydhms(1, :))), q(1, ~isnan(ydhms(1, :))));

    if nargout > 1

        % Finds the precision
        seconds_precision = q(1, (find(isnan([ydhms(1, :), NaN]), 1, 'first') - 1));

    end

end