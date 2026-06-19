function [v] = paramsToVector(p, e, lan, incl, aop)
%
% function [v] = paramsToVector(p, e, lan, incl, aop)
%
% Compresses all orbital parameters in a single column vector

    v = [p; e; lan; incl; aop];

end
