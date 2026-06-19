function cc = MinmusOrbitCost01(b)

load MinmusOrbitData02.mat

    thEM = -acos(((oOM.p + aOMp * a(1, 1)) / SOIr_M - 1) / (oOM.e + 0.0001 * a(2, 1)));
    eT = argToTime(thEM, (oOM.p + aOMp * a(1, 1)), (oOM.e + 0.0001 * a(2, 1)), ...
    mu_M, (oOM.peT + 61 * a(6, 1)));
    [xxEM, vvEM] = paramsToPosVel(thEM, (oOM.p + aOMp * a(1, 1)), ...
            (oOM.e + 0.0001 * a(2, 1)), (oOM.lan + aq * a(3, 1)), ...
            (oOM.incl + aq * a(4, 1)), (oOM.aop + aq * a(5, 1)), mu_M);


    [thEK] = timeToArg(eT, (oOK.p + aOKp * b(1, 1)), (oOK.e + 0.0001 * b(2, 1)), ...
             mu_K, (oOK.peT + 61 * b(6, 1)));
    [xxEK, vvEK] = paramsToPosVel(thEK, (oOK.p + aOKp * b(1, 1)), ...
            (oOK.e + 0.0001 * b(2, 1)), (oOK.lan + aq * b(3, 1)), ...
            (oOK.incl + aq * b(4, 1)), (oOK.aop + aq * b(5, 1)), mu_K);

    [thMK] = timeToArg(eT, oMK.p, oMK.e, mu_K, oMK.peT);
    [xxMK, vvMK] = paramsToPosVel(thMK, oMK.p, oMK.e, oMK.lan, ...
            oMK.incl, oMK.aop, mu_K);

    cc = ((norm(xxEK - xxMK, 2) - SOIr_M)) ^ 2 + ...
            (norm(((xxEK - xxMK) - xxEM), 2)) ^ 2  + ...
            (SOIr_M * (norm(((vvEK - vvMK) - vvEM), 2) / norm(vvMK, 2))) ^ 2;

    %for idx_1 = 1:size(a, 1)
    %    cc = cc + ((a(idx_1, 1)) ^ 4);
    %end

end

