function [g] = Pulaski_I_03_grad2(cf, a, da, nlize)
    if nargin < 4
        nlize = true;
    end
    N = size(a, 2);
    g = zeros(1, N);
    for idx_1 = 1:N
        q = zeros(1, N);
        q(1, idx_1) = da;
        g(1, idx_1) = (cf(a(1, 2:end) + q(1, 2:end), 382 * (a(1, 1) + ...
                      q(1, 1))) - cf(a(1, 2:end) - q(1, 2:end), 382 * ...
                      (a(1, 1) - q(1, 1)))) / (2 * da);
    end
    if nlize
        g = g / norm(g, 2);
    end
end

