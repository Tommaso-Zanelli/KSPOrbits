function [p, e] = apoPeriToPe(apo, peri)
%
%   function [p, e] = apoPeriToPe(apo, peri)
%
%   Computes the semi-latus rectum and the eccentricity from the apoapsis
%   and the periapsis
%

    detpe = peri + apo;
    p     = 2 * peri.* apo./ detpe;
    e     = (apo - peri)./ detpe;

end 