function [T] = argStructToTime(argument, o, mu, att)
%
% function [T] = argStructToTime(argument, o, mu, att)
%
% Computes the time given the argument and the orbital parameters
%

    % Additional arguments for "argToTime"
    if nargin < 4
    	att = [];
    end
    att = check_att(att);

    % Computes the time
    [T] = argToTime(argument, o.p, o.e, mu, o.peT, att.ep, att.d);

end