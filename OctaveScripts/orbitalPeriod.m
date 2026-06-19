function [T] = orbitalPeriod(p, e, mu)

    T = 2 * pi * sqrt((p ^ 3) / (mu * ((1 - (e ^ 2)) ^ 3)));

end    