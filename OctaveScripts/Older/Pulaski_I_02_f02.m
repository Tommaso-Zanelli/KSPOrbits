function A = Pulaski_I_02_f02(cfb1, T, da1, da2, h)

    cfbL = @(A) cfb1(A, T, da1, da2, h);
    np = 16;
    cc = 1e308;
    while ((cc > 1e-11) && np < (2 ^ 16))

        A = optiFun([-pi, pi], cfbL, 4e-6, np, 1e-12, 1000);
        cc = cfbL(A);
        np = np * 2;

    end

end