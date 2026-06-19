function [g] = Pulaski_I_03_grad(cf, a, T, da)
    N = size(a, 2);
    g = zeros(1, N);
    for idx_1 = 1:N
        q = zeros(1, N);
        q(1, idx_1) = da;
        g(1, idx_1) = (cf(a + q, T) - cf(a - q, T)) / (2 * da);
    end
    g = g / norm(g, 2);
end

