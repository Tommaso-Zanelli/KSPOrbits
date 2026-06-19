function [c] = orbitCompareCost(refv, orb)
%
%   function [c] = orbitCompareCost(refv, orb)
%
%   Compares a precisely defined orbit's parameters with reference values.
%
%   Comparison is done through the cost function cf0p5, which assigns a
%   negligible cost for values under 0.5, and an increasing exponentially
%   one after 0.5
%

    % Gets apoapsis and periapsis of the orbit
    [ap_o, pe_o] = peToApoPeri(orb.p, orb.e);

    % Computes the cost
    c = cf0p5((refv.ap - ap_o) * (10 ^ (5 - floor(log10(refv.ap))))) + ...  % Apoapsis comparison. The significant decimal places are 6.
        cf0p5((refv.pe - pe_o) * (10 ^ (5 - floor(log10(refv.pe))))) + ...  % Periapsis comparison. The significant decimal places are 6.
        cf0p5((refv.e - orb.e) * 1e4) + ...                                 % Eccentricity comparison. The significant decimal places are 4, assuming elliptic orbit
        cf0p5((refv.lan - orb.lan * 57.295779513082323481365) * 10) + ...   % Longitude of ascending node comparison. The last significant decimal place is 0.1 degrees.
        cf0p5((refv.incl - orb.incl * 57.295779513082323481365) * 10) + ... % Inclination comparison. The last significant decimal place is 0.1 degrees.
        cf0p5((refv.aop - orb.aop * 57.295779513082323481365) * 10);        % Argument of periapsis comparison. The last significant decimal place is 0.1 degrees.

end