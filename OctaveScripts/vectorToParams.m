function [p, e, lan, incl, aop] = vectorToParams(v)
%
% function [p, e, lan, incl, aop] = vectorToParams(v)
%
% Extracts all orbital parameters from a column vector

    if size(v, 1) < 5
        fprintf('Runtime Error: wrong vector size!\n')
        p    = NaN;
        e    = NaN;
        lan  = NaN;
        incl = NaN;
        aop  = NaN;
        return
    end

    p    = v(1, 1);
    e    = v(2, 1);
    lan  = v(3, 1);
    incl = v(4, 1);
    aop  = v(5, 1);

end