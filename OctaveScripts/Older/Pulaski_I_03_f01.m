function [r] = Pulaski_I_03_f01(o0, ttm_offs, dvm, Nm1, Nm2, mu_K, A, T, da1, da2, h)

    [os, ths] = applyMultipleManoeuvres(o0, T + ttm_offs(1, 1:(Nm1 + Nm2)), ...
                dvm(1, 1:(Nm1 + Nm2)), (A * ones(1, (Nm1 + Nm2)) + da1), ...
                zeros(1, (Nm1 + Nm2)) + da2, mu_K);
    ths(2, (Nm1 + Nm2 + 1)) = timeToArg(ttm_offs(1, (Nm1 + Nm2 + 1)), ...
            os(1, (Nm1 + Nm2 + 1)).p, os(1, (Nm1 + Nm2 + 1)).e, mu_K, ...
            os(1, (Nm1 + Nm2 + 1)).peT);
    [x, v] = paramsToPosVel(ths(2, (Nm1 + Nm2 + 1)), ...
               os(1, (Nm1 + Nm2 + 1)).p, os(1, (Nm1 + Nm2 + 1)).e, ...
               os(1, (Nm1 + Nm2 + 1)).lan, os(1, (Nm1 + Nm2 + 1)).incl, ...
               os(1, (Nm1 + Nm2 + 1)).aop, mu_K);
    [lan, incl, aop_th] = rotAngles(x, [-x(2, 1); x(1, 1); 0]);
    [xc, vc] = paramsToPosVel(0, norm(x, 2), 0, lan, incl, aop_th, mu_K);

    % Tests for
    c0 = 0;
    for idx_1 = 1:(Nm1 + Nm2 + 1)
        [xp, ~] = paramsToPosVel(ths(2, idx_1), os(1, idx_1).p, ...
                  os(1, idx_1).e, os(1, idx_1).lan, os(1, idx_1).incl, ...
                  os(1, idx_1).aop, mu_K);
        c0 = max([c0, (0.5e16 + 0.5e16 *tanh(0.15*(670000 - norm(xp, 2))))]);
    end

    if norm(x - xc, 2) > 1e-5
       fprintf('AIUTO!!!!\n')
    end
    c1 = ((norm(x(1:2, 1), 2) - h)) ^ 4;
    c2 = (x(3, 1)) ^ 4;
    c3 = norm((v - vc), 2);

    r = c0 + c1 + c2 + c3;

end