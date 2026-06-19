%% ------------------------------------------------------------------------

fclose all;
close all
format long g
more off
clear
clc

cyc = true;

%{
so  =   [-0.00960582483485628
        -0.213036492396545
         0.342829239016859
        -0.211581662397275
        -0.397413742283452
        0.0557158645286978
         0.249107766765123
         0.469446517656909
         0.499064728652236
         0.290928099907113
         0.222008384988038
        -0.140032273607815
        -0.293199998484641
        -0.456524796195334
         0.497435702314548
        -0.112465746536519];
%}

cn = 0;
load optDunaShot001.mat
so = cs;
%cs = so;
%cc = 0.000266984986423081;
cp = cc;
gain = 0;
minD = 4;
last5gains = zeros(5, 1);
gainActive = false;
opt = true;
cqo = ones(7, 1);

%so = zeros(15, 1);
d = 2;
fprintf('%gm, %gs (d = %g); %g\n', cc * 6000, cc * 2000, d, gain)

while opt

if cyc
    s = so + (2 ^ (-d)) * (rand(16, 1) - 0.5);
    s(s >= 0.5) = (0.5 - eps(0.5));
    s(s <= -0.5) = (eps(0.5) - 0.5);
else
    s = so;
end

% -------------------------------------------------------------------------

degToRad = (pi / 180);
radToDeg = 1 / degToRad;

R_K  = 600000;
mu_K = 3.5316e12;

apo_1 = 153326 + s(1, 1);
per_1 = 152275 + s(2, 1);
inc_1 = (0     + 0.1 * s(3, 1)) * degToRad;
lan_1 = (195.3 + 0.1 * s(4, 1)) * degToRad;
aop_1 = (332.8 + 0.1 * s(5, 1)) * degToRad;

apT_1  = [10, 46 + s(6, 1)];
peT_1  = [28, 58 + s(7, 1)];
peT_2  = [55, 36 + s(8, 1)];
T = [0, 122, 0, 46, 50 + s(9, 1)];
T2 = [0, 122, 0, 46, 51 + s(10, 1)];
Tp_1 = [36, 23 + s(11, 1)];

apo_2 = -804139 + s(12, 1);
per_2 = -467932 + s(13, 1);
inc_2 = (0.4   + 0.1 * s(14, 1)) * degToRad;
lan_2 = (236.1 + 0.1 * s(15, 1)) * degToRad;
aop_2 = (274.0 + 0.1 * s(16, 1)) * degToRad;

T_M = 2641278.000;
v_M = [0; -12.4; 10133.8];

% -------------------------------------------------------------------------

apT_1s = dot(apT_1, [60, 1]);
peT_1s = dot(peT_1, [60, 1]);
peT_2s = dot(peT_2, [60, 1]);
Tp_1s  = dot(Tp_1, [60, 1]);
T_s    = dot(T, [0, 21600, 3600, 60, 1]);
T2_s    = dot(T2, [0, 21600, 3600, 60, 1]);

Tpe_1 = T_s + peT_1s;
Tpe_2 = T2_s + peT_2s;

[p_1, e_1] = apoPeriToPe(apo_1 + R_K, per_1 + R_K);
[p_2, e_2] = apoPeriToPe(apo_2 + R_K, per_2 + R_K);

Tp_1sc1 = orbitalPeriod(p_1, e_1, mu_K);
Tp_1sc2 = 2 * abs(apT_1s - peT_1s);

% -------------------------------------------------------------------------

M_in_1  = T_M - T_s;
M_in_2  = T_M - T2_s;
T_M_1 = T_M - Tpe_1;
T_M_2 = T_M - Tpe_2;

th_1M = timeToArg(T_M_1, p_1, e_1, mu_K);
th_2M = timeToArg(T_M_2, p_2, e_2, mu_K);
[x_1M, v_1M] = paramsToPosVel(th_1M, p_1, e_1, lan_1, inc_1, aop_1, mu_K);
[x_2M, v_2M] = paramsToPosVel(th_2M, p_2, e_2, lan_2, inc_2, aop_2, mu_K);


i_green = v_1M / norm(v_1M, 2);
i_r = x_1M / norm(x_1M, 2);
i_cyan  = i_r - dot(i_r, i_green) * i_green;
i_cyan  = i_cyan / norm(i_cyan, 2);
i_purple = cross(i_cyan, i_green);

x_2Mc = x_1M;
v_2Mc = v_1M + v_M(1, 1) * i_green - v_M(3, 1) * i_cyan + v_M(2, 1) * i_purple;

[th_2c, p_2c, e_2c, lan_2c, inc_2c, aop_2c] = posVelToParams(x_2Mc, v_2Mc, mu_K);

% -------------------------------------------------------------------------

if ~cyc

    tth1 = linspace(-pi, pi, 5001);
    l2 = 0.55;
    tth2 = linspace(-l2 * pi, l2 * pi, 5001);
    [xx1, ~] = paramsToPosVel(tth1, p_1, e_1, lan_1, inc_1, aop_1, mu_K);
    [xx2, ~] = paramsToPosVel(tth2, p_2, e_2, lan_2, inc_2, aop_2, mu_K);

    [xx2o, ~] = paramsToPosVel(tth2, p_2c, e_2c, lan_2c, inc_2c, aop_2c, mu_K);

    figure(1)
    plot3(xx1(1, :), xx1(2, :), xx1(3, :), 'b', 'LineWidth', 2)
    hold on
    plot3(x_2M(1, :), x_2M(2, :), x_2M(3, :), 'ro', 'LineWidth', 2)
    plot3(xx2(1, :), xx2(2, :), xx2(3, :), 'r', 'LineWidth', 2)
    plot3(x_1M(1, :), x_1M(2, :), x_1M(3, :), 'b*', 'LineWidth', 2)
    plot3(xx2o(1, :), xx2o(2, :), xx2o(3, :), 'm', 'LineWidth', 1)

    [Msph_X, Msph_Y, Msph_Z] = genSphere(R_K, 51);
    surf(Msph_X, Msph_Y, Msph_Z, zeros(size(Msph_X)));

    axis equal
    hold off
    opt = false;

end

cq = ([(1e0) * abs((Tp_1s - Tp_1sc1) / Tp_1sc1), ...
       (1e0) * abs((Tp_1sc2 - Tp_1sc1) / Tp_1sc1), ...
     (1e2) * norm((x_2M - x_1M), 2) / R_K, ...
     (1e4) * abs((p_2c - p_2) / p_2), (1e4) * abs((e_2c - e_2) / e_2), ...
     abs(((inc_2 - 2 * pi * (inc_2 > pi)) - (inc_2c - 2 * pi * (inc_2c > pi))) / (inc_2 - 2 * pi * (inc_2 > pi))), ...
     abs(((lan_2 - 2 * pi * (lan_2 > pi)) - (lan_2c - 2 * pi * (lan_2c > pi))) / (lan_2 - 2 * pi * (lan_2 > pi))), ...
     abs(((aop_2 - 2 * pi * (aop_2 > pi)) - (aop_2c - 2 * pi * (aop_2c > pi))) / (aop_2 - 2 * pi * (aop_2 > pi)))]');

c = max(cq);% + mean(cq);

if cyc

    if c < cc
        cp = cc;
        cc = c;
        cqo = cq;
        gain = (abs(cp - cc) / cc);
        if ~gainActive && gain
            last5gains = repmat(gain, 5, 1);
            minD = d;
            gainActive = true;
        else
            last5gains(1:4, 1) = last5gains(2:5, 1);
            last5gains(5, 1) = gain;
            if max(last5gains) < 1e-4
               d = max([2, (minD - 2)]);
               gainActive = false;
            end
        end
        fprintf('%gm, %gs (d = %g); %g\n', cc * 6000, cc * Tp_1sc1, d, gain)
        cs = s;
        cn = 0;
    else
        cn = cn + 1;
    end

    if cn > 2500
       cn = 0;
       so = cs;
       d = d + 1;
       fprintf('%gm, %gs (d = %g); %g\n', cc * 6000, cc * Tp_1sc1, d, gain)
    end

end

end

if cyc
    save optDunaShot001.mat cs cc
    cs
    cc
end