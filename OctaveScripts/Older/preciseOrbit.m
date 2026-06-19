function orb0 = preciseOrbit(ref0, times0, refvs, mnvrs, mu)

    [a(1, 1), a(1, 2)] = apoPeriToPe(refv.ap, refv.pe);
    a(1, 3) = refv.lan * 0.01745329251994329701382;
    a(1, 4) = refv.incl * 0.01745329251994329701382;
    a(1, 5) = refv.aop * 0.01745329251994329701382;
    a(1, 6) = times0.pe;

    cf = @(a) preciseOrbitCost(paramsToStruct(a(1, 1), a(1, 2), ...
                               a(1, 3), a(1, 4), a(1, 5), a(1, 6)), ...
                               ref0, times0, refvs, mnvrs, mu);

    % Place optimization cycle here

	orb0 = paramsToStruct(a(1, 1), a(1, 2), a(1, 3), a(1, 4), a(1, 5), a(1, 6));

end