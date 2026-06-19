fclose all;
close all
format long g
more off
clear
clc

%% ------------------------------------------------------------------------
%  Minmus orbit

R_M    = 60000;
mu_M   = 1.7658e9;
secs   = [21600, 3600, 60, 1, 0];

roM0   = [-9.74286e4,-3.33373e4, 5.9532, 88.2, 65.1, 279.4];
csroM0 = [(10.^ (floor(log10(abs(roM0(1, 1:2)))) - 5)), 0.0001, 0.1 * ones(1, 3)];
ro0peT = [  2,   4,  40, NaN,   +1];
ro0cmT = [139,   0,  31,   6,   +1];

cs0peT = secs(1, min([(find(isnan(ro0peT)) - 1), 4])) + secs(1, min([(find(isnan(ro0cmT)) - 1), 4]));

moM0 = [-565.83,   -41.98,    13.784, 3062368.521,  2.13469e6,  9.72169e5,   0.3603,  90.1, 136.5, 228.3
        -607.058,   -8.063,  -22.188, 3062773.779,  1.22825e6,  9.25680e4,   0.7882, 267.4, 134,   206
        -637.327,   45.298,  150.805, 3063139.719, -8.31200e5,  5.59070e5,   9.1387,  87.3,  48.1, 157.1
        -478.238,  -21.771,  -37.485, 3061678.447, -9.60052e5,  4.42155e5,   3.524,  267.3, 146.5,  80.4
        -503.556,  -68.774, -105.424, 3060958.418, -1.90737e6,  1.60177e6,  18.9078, 267.5, 148.3,  30.4
        -579.447,  -35.277,  -16.871, 3062211.968,  1.38614e6,  1.08948e6,   0.1143, 212.3, 178.6, 241.9
        -576.017,   21.576,   24.529, 3061372.629,  1.79892e6,  1.57264e6,   0.0648,  86.9,  23.7, 113.9
        -590.395,   35.425,  -68.538, 3062285.351, -2.32597e6,  1.09067e6,   3.0634, 268.8,  87.7, 346
         388.108,  481.912,  -93.243, 3062285.351, -5.89891e5,  4.66863e5, 348.9618, 270.6,  32.9,  61.2
         -51.82,  -247.087,  205.262, 3063401.242, -3.60774e5,  2.31330e5,  62.6954,  90.7, 111.3, 240.2];
N      = size(moM0, 1);
csmoM0 = [0.001 * ones(N, 4), (10.^ (floor(log10(abs(moM0(:, 5:6)))) - 5)), 0.0001 * ones(N, 1), 0.1 * ones(N, 3)];

[oM0, ecM0] = rawOrbitToStruct(roM0, R_M);
oM0.peT     = ro0peT(1, 5) * dot(ro0peT(1, ~isnan(ro0peT)), secs(1, ~isnan(ro0peT))) + ...
              ro0cmT(1, 5) * dot(ro0cmT(1, ~isnan(ro0cmT)), secs(1, ~isnan(ro0cmT)));
for idx_1 = 1:N
    [oMm(idx_1, 1), ecMm(idx_1, 1)] = rawOrbitToStruct(moM0(idx_1, 5:size(moM0, 2)), R_M);
end


roM1   = [-9.74286e4,-3.33373e4, 5.9532, 88.2, 65.1, 279.4, oM0.peT];

%{
[oM1, ecM1] = rawOrbitToStruct(roM1(1, 1:6), R_M);
oM1.peT = roM1(1, 7);
%emx = 0;
for idx_1 = 1:N
    thm1(idx_1, 1) = timeToArg(moM0(idx_1, 4), oM1.p, oM1.e, mu_M, oM1.peT);
    [xm(:, idx_1), vm1(:, idx_1)] = paramsToPosVel(thm1(idx_1, 1), oM1.p, oM1.e, ...
                                                   oM1.lan, oM1.incl, oM1.aop, mu_M);
    [T] = refFrameLocal(xm(:, idx_1), vm1(:, idx_1));
    dvm(:, idx_1) = [0; 0; 0];
    for idx_2 = 1:3
        dvm(:, idx_1) = dvm(:, idx_1) + T(:, idx_2) * moM0(idx_1, idx_2);
    end
    vm2(:, idx_1) = vm1(:, idx_1) + dvm(:, idx_1);
    [thm2(idx_1, 1), oMm2(idx_1, 1).p, oMm2(idx_1, 1).e, oMm2(idx_1, 1).lan, ...
     oMm2(idx_1, 1).incl, oMm2(idx_1, 1).aop] = posVelToParams(xm(:, idx_1), ...
                                                    vm2(:, idx_1), mu_M);
    oMm2(idx_1, 1).peT = 0;
    [cap, cpe] = peToApoPeri(oMm2(idx_1, 1).p, oMm2(idx_1, 1).e);
    moMc(idx_1, :) = [moM0(idx_1, 1:4), (cpe - R_M), (cap - R_M), ecMm(idx_1, 1), ...
        (57.29577951308232017687) * [oMm2(idx_1, 1).lan, oMm2(idx_1, 1).incl, ...
        oMm2(idx_1, 1).aop]];
    %{
    fprintf('\nManoeuvre %g:\n', idx_1)
    fprintf('p: | %g - %g | / | %g | = %g\n', oMm(idx_1, 1).p, ...
        oMm2(idx_1, 1).p, oMm(idx_1, 1).p, abs(oMm(idx_1, 1).p - ...
        oMm2(idx_1, 1).p) / abs(oMm(idx_1, 1).p))
    emx = max([emx, abs(oMm(idx_1, 1).p - ...
        oMm2(idx_1, 1).p) / abs(oMm(idx_1, 1).p)]);
    fprintf('e: | %g - %g | = %g\n', oMm(idx_1, 1).e, oMm2(idx_1, 1).e, ...
        abs(oMm(idx_1, 1).e - oMm2(idx_1, 1).e) / max([1, abs(oMm(idx_1, 1).e), ...
        abs(oMm2(idx_1, 1).e)]))
    emx = max([emx, abs(oMm(idx_1, 1).e - oMm2(idx_1, 1).e) / max([1, abs(oMm(idx_1, 1).e), ...
        abs(oMm2(idx_1, 1).e)])]);
    fprintf('lan : | %g - %g | / (2 * pi) = %g\n', oMm(idx_1, 1).lan, ...
        oMm2(idx_1, 1).lan, abs(oMm(idx_1, 1).lan - oMm2(idx_1, 1).lan) / (2 * pi))
    emx = max([emx, abs(oMm(idx_1, 1).lan - oMm2(idx_1, 1).lan) / (2 * pi)]);
    fprintf('incl: | %g - %g | / (2 * pi) = %g\n', oMm(idx_1, 1).incl, ...
        oMm2(idx_1, 1).incl, abs(oMm(idx_1, 1).incl - oMm2(idx_1, 1).incl) / (2 * pi))
    emx = max([emx, abs(oMm(idx_1, 1).incl - oMm2(idx_1, 1).incl) / (2 * pi)]);
    fprintf('aop : | %g - %g | / (2 * pi) = %g\n', oMm(idx_1, 1).aop, ...
        oMm2(idx_1, 1).aop, abs(oMm(idx_1, 1).aop - oMm2(idx_1, 1).aop) / (2 * pi))
        emx = max([emx, abs(oMm(idx_1, 1).aop - oMm2(idx_1, 1).aop) / (2 * pi)]);
    %}
end

e_raw1 = moM0 - moMc;
for idx_1 = 1:N
   for idx_2 = 8:10
       while e_raw1(idx_1, idx_2) > 180
           e_raw1(idx_1, idx_2) = e_raw1(idx_1, idx_2) - 360;
       end
       while e_raw1(idx_1, idx_2) < -180
           e_raw1(idx_1, idx_2) = e_raw1(idx_1, idx_2) + 360;
       end
   end
end
e_raw1 = (e_raw1(:, 5:end)) ./ csmoM0(:, 5:end);
e_raw2 = (roM1 - [roM0, oM0.peT])./ [csroM0, cs0peT];

ee = sum(sum(c05(e_raw1), 1), 2) + sum(c05(e_raw2), 2);
%}

%{
roM1_2   = [-97428.9916412538, -33338.0394041494, 5.95281065790811, 90.2572119828415, ...
                65.0164413072883, 279.405444594618, 3064300.96602154];

cf = @(x) Bashir_2_costf1(x, true);
dx = [0.00390625, 0.001953125, 2.38418579101562e-07, 3.814697265625e-06, 3.814697265625e-06, 1.52587890625e-05, 0.125]' * 128;
cc = optiGrad(cf, roM1_2', 1e-12, dx, 15);
%}



x0 = [roM1'; reshape(moM0(:, 1:4), 40, 1)];
c0 = [csroM0'; cs0peT; reshape(csmoM0(:, 1:4), 40, 1)];
%u0 = [zeros(size(roM1, 2), 1); 1; zeros(40, 1)];
u0 = zeros(size(x0));

%cc(:, 1) = x0;

%load Bashir_2_IntRes.mat
%cc(:, 1) = x;

%for idx_1 = 1:25
%    idx_2 = 2 * idx_1;
%    idx_3 = 2 * idx_1 + 1;
%    n = 10 ^ (-3 - floor((idx_1 - 1) / 3));
%    cf = @(x) Bashir_2_costf1(x, 0);
%    cc(:, idx_2) = optiGrad(cf, cc(:, (idx_2 - 1)), n, [], 30);
%    cf = @(x) Bashir_2_costf1(x, 1);
%    cc(:, idx_3) = optiGrad(cf, cc(:, (idx_3 - 1)), n, [], 30);
%end

%% Now for some Monte Carlo fun

N = size(x0, 1);

x_max = x0 + (0.5 * (u0 + 1) - eps(0.5 * (u0 + 1))).* c0;
x_min = x0 + 0.5 * (u0 - 1).* c0;


if exist('Bashir_2_cx.mat', 'file')
    load Bashir_2_cx.mat
else
    cx = x0;
end
ci = Bashir_2_costf1(cx, 0);
fprintf('%g\n', ci)

n1 = 0;
n2 = 0;
rr = 1;

N1 = 1000;
N2 = 10;

while true

    n1 = n1 + 1;
    x_rand = cx + (rand(N, 1) - 0.5).* c0 * rr;
    x_rand(x_rand > x_max) = x_max(x_rand > x_max);
    x_rand(x_rand < x_min) = x_max(x_rand < x_min);
    cc = Bashir_2_costf1(x_rand, 0);
    if cc < ci
        fprintf('%g\n', cc)
        n2 = n2 + 1;
        n1 = 0;
        ci = cc;
        cx = x_rand;
    end

    if n1 > N1
        fprintf('Current rr: %g\n', rr)
        rr = 0.5 * rr;
        n1 = 0;
    end

    if rr < (2 ^ -08)
        rr = 1;
    end

end

Bashir_2_costf1(x_min, 1)
Bashir_2_costf1(x0, 1)
Bashir_2_costf1(x_max, 1)
