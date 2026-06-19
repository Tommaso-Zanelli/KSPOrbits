function [] = showOrbit(o, r)
    [ap, pe] = peToApoPeri(o.p, o.e);
    fprintf('Periapsis: %7.5e\n', ap - r)
    fprintf('Apoapsis : %7.5e\n', pe - r)
    fprintf('e        : %9.4f\n', o.e)
    fprintf('lan      : %5.1f\n', o.lan * 180 / pi)
    fprintf('incl     : %5.1f\n', o.incl * 180 / pi)
    fprintf('aop      : %5.1f\n', o.aop * 180 / pi)


end