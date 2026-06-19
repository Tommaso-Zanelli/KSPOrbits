%%

fclose all;
close all
format long g
more off
clear
clc

%%

N = 2;
Nt = 6;

M = rand(N);
M = 0.5 * (M + M');
[R, lmb] = eig(M);
lmb = -abs(lmb);
M = R * lmb * R';

x0 = rand(N, 1);

xi0 = (R') * x0;

T0 = 0;
T1 = 1;

er = zeros(1 + N, Nt);

ee = zeros(1 + N, Nt);

erk2 = zeros(1 + N, Nt);

for idx_1 = 1:Nt

    fprintf('Step %g of %g\n', idx_1, Nt)

    Ntt = 2 ^ (idx_1 + 4);

    dt = 1 / (Ntt - 1);

    f = @(xx, t) M * xx;

    [xr, tr] = rk4_opt(f, x0, dt, T1, T0);

    xe = zeros(N, Ntt);
    xrk2 = zeros(N, Ntt);
    xe(:, 1) = x0;
    xrk2(:, 1) = x0;
    for idx_2 = 2:Ntt
      xe(:, idx_2) = xe(:, (idx_2 - 1)) + dt * f(xe(:, (idx_2 - 1)), tr(1, (idx_2 - 1)));
      xrk2(:, idx_2) = xrk2(:, (idx_2 - 1)) + dt * f(xrk2(:, (idx_2 - 1)), tr(1, (idx_2 - 1)));
      xrk2(:, idx_2) = xrk2(:, (idx_2 - 1)) + dt * f(0.5 * (xrk2(:, (idx_2 - 1)) + xrk2(:, idx_2)), 0.5 * (tr(1, (idx_2 - 1)) + tr(1, (idx_2 - 1))));
    end

    xi = exp(repmat(diag(lmb), 1, Ntt).* repmat(tr, N, 1)).* repmat(xi0, 1, Ntt);

    x = R * xi;

    er(1, idx_1) = dt;

    er(2:(N + 1), idx_1) = max((abs(x - xr)./ repmat(max(abs(x)')', 1, Ntt))')';

    ee(1, idx_1) = dt;

    ee(2:(N + 1), idx_1) = max((abs(x - xe)./ repmat(max(abs(x)')', 1, Ntt))')';

    erk2(1, idx_1) = dt;

    erk2(2:(N + 1), idx_1) = max((abs(x - xrk2)./ repmat(max(abs(x)')', 1, Ntt))')';

end

pp = zeros(N, 2);

for idx_1 = 1:N

    pq = polyfit(log(er(1, :)), log(er((idx_1 + 1), :)), 1);
    pp(idx_1, :) = pq;

end

pe = zeros(N, 2);

for idx_1 = 1:N

    pq = polyfit(log(ee(1, :)), log(ee((idx_1 + 1), :)), 1);
    pe(idx_1, :) = pq;

end

prk2 = zeros(N, 2);

for idx_1 = 1:N

    pq = polyfit(log(erk2(1, :)), log(erk2((idx_1 + 1), :)), 1);
    prk2(idx_1, :) = pq;

end