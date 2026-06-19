%% - Header ---------------------------------------------------------------

fclose all;
close all
format long g
more off
clear
clc

%% - Veichle capacity -----------------------------------------------------

dvt1 = [0 10 20 30 40 50 60 70 76.5900000000000;0 25 51 77 103 129 155 180 197]';
dvt2 = [0 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200 210 220 230 240 250 260 270 280 290 300 310 320 330 340 350 360 370 380 390 400 410 420 430 440 450 460 470 480 490 500 510 520 530 540 550 560 570 580 590 600 610 620 630 630.540000000000;0 22 44 66 88 110 132 154 176 198 220 241 263 285 307 328 350 372 394 415 437 458 480 501 523 544 566 587 609 630 651 672 694 715 736 757 779 800 821 842 863 884 905 926 947 968 989 1009 1030 1051 1072 1093 1113 1134 1155 1175 1196 1217 1237 1258 1278 1299 1319 1340 1341]';
dvt3 = [50 100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300;4 9 14 18 23 27 31 36 40 44 48 52 56 60 64 67 71 75 78 82 85 88 92 95 98 101]';

g0    = 9.81;
Isp   = 800;
Isp3  = 330;
dmdt  = 1.53 * 5;
dmdt3 = 18.54 * (100 / 9);

mr1 = exp((dvt1(2:end, 1) - dvt1(1:(end - 1), 1)) / (Isp * g0));
mr2 = exp((dvt2(2:end, 1) - dvt2(1:(end - 1), 1)) / (Isp * g0));
mr3 = exp((dvt3(2:end, 1) - dvt3(1:(end - 1), 1)) / (Isp3 * g0));


M1  = 154788;
dM2 =  20442;
dM3 = 56166;
c2m = 1e308;

mm1 = zeros(size(dvt1, 1), 3);
mm2 = zeros(size(dvt2, 1), 3);
mm3 = zeros(size(dvt3, 1), 3);

mm1(1, 1) = M1;

for idx_1 = 2:size(dvt1, 1)
    mm1(idx_1, 1) = mm1((idx_1 - 1), 1) / mr1((idx_1 - 1), 1);
    mm1(idx_1, 2) = (mm1((idx_1 - 1), 1) - mm1(idx_1, 1)) / dmdt;
    mm1(idx_1, 3) = mm1((idx_1 - 1), 3) + mm1(idx_1, 2);
end

mm2(1, 1) = mm1(end, 1) - dM2;

for idx_1 = 2:size(dvt2, 1)
    mm2(idx_1, 1) = mm2((idx_1 - 1), 1) / mr2((idx_1 - 1), 1);
    mm2(idx_1, 2) = (mm2((idx_1 - 1), 1) - mm2(idx_1, 1)) / dmdt;
    mm2(idx_1, 3) = mm2((idx_1 - 1), 3) + mm2(idx_1, 2);
end

c1 = abs(floor(mm1(:, 3) + 0.5) - dvt1(:, 2));
c2 = abs(floor(mm2(:, 3) + 0.5) - dvt2(:, 2));

c2t = sum(c2);
c1t = sum(c1);

%dM3s =56165:0.0001:56175;
%ctt = zeros(size(dM3s));
%idx_2 = 1;

%for dM3 = dM3s

mm3(1, 1) = mm2(end, 1) - dM3;

for idx_1 = 2:size(dvt3, 1)
    mm3(idx_1, 1) = mm3((idx_1 - 1), 1) / mr3((idx_1 - 1), 1);
    mm3(idx_1, 2) = (mm3((idx_1 - 1), 1) - mm3(idx_1, 1)) / dmdt3;
    mm3(idx_1, 3) = mm3((idx_1 - 1), 3) + mm3(idx_1, 2);
end

c3 = abs(floor(mm3(:, 3) + 0.5) - dvt3(:, 2));
%c3 = abs((mm3(:, 3)) - dvt3(:, 2));


c3t = sum(c3);
%ctt(1, idx_2) = c3t;
%idx_2 = idx_2 + 1;

%end

%% - Timings --------------------------------------------------------------

Dt = 5;
Dt2 = 55;

tt1 = [49 49 49 50];
tt2 = [53 53 53 53 53 53 53 53 54 54 54 54 54 54 54 54 54 54 54 54 54 54 54 54 54];
tt3 = [46, 46];

Nm = (size(tt1, 2) + size(tt2, 2) + size(tt3, 2));

dtss = 0.5 * [(tt1(1, 2:end) + tt1(1, 1:(end - 1))), ...
              (tt1(1, end) + tt2(1, 1) + Dt2 - Dt), ...
              (tt2(1, 2:end) + tt2(1, 1:(end - 1))), ...
              (tt2(1, end) + tt3(1, 1) + Dt2 - Dt), ...
              (tt3(1, 2:end) + tt3(1, 1:(end - 1)))] + Dt;

tts = zeros(1, (size(dtss, 2) + 1));
dts = [tt1, tt2, tt3];
mmm = zeros(1, Nm);
MM = zeros(1, (Nm + 1));
Dvs = zeros(1, Nm);

for idx_1 = 1:size(dtss, 2)
    tts(1, (idx_1 + 1)) = tts(1, idx_1) + dtss(1, idx_1);
end

stage = [zeros(1, size(tt1, 2)), ones(1, size(tt2, 2)), 2 * ones(1, size(tt3, 2))];
mmm(1, (stage < 2)) = dts(1, (stage < 2)) * dmdt;
mmm(1, (stage >= 2)) = dts(1, (stage >= 2)) * dmdt3;

MM(1, 1) = M1;
Isc = Isp;
for idx_1 = 2:(Nm + 1)
    MM(1, idx_1) = MM(1, (idx_1 - 1)) - mmm(1, (idx_1 - 1));
    Dvs(1, (idx_1 - 1)) = Isc * g0 * log(MM(1, (idx_1 - 1)) / MM(1, idx_1));
    if (idx_1 == 5)
        MM(1, idx_1:end) = MM(1, idx_1:end) - dM2;
    end
    if (idx_1 == 30)
        MM(1, idx_1:end) = MM(1, idx_1:end) - dM3;
        Isc = Isp3;
    end
end

%% - Orbit details --------------------------------------------------------

mu_K  = 3.5316e12;
R_K = 600000;

ap0   = -300.710e6 + R_K;
pe0   = 74.9998e3  + R_K;
incl0 = 7.0   * (pi / 180);
lan0  = 147.0 * (pi / 180);
aop0  = 34.4  * (pi / 180);
pe0T  = 2856546.578;

apf   = 75.0e3 + R_K;
pef   = 75.0e3 + R_K;
inclf = 0.0 * (pi / 180);
lanf  = 0.0 * (pi / 180);
aopf  = 0.0 * (pi / 180);

[p0, e0] = apoPeriToPe(ap0, pe0);
[pf, ef] = apoPeriToPe(apf, pef);

[thl1] = timeToArg((pe0T - 4320), p0, e0, mu_K, pe0T);
[thl2] = timeToArg((pe0T + 4320), p0, e0, mu_K, pe0T);
th0s = linspace(thl1, thl2, 5001);

thfs = linspace(-pi, pi, 5001);

[x0, v0] = paramsToPosVel(th0s, p0, e0, lan0, incl0, aop0, mu_K);


%[xf, vf] = paramsToPosVel(thfs, pf, ef, lanf, inclf, aopf, mu_K);

th0fM = -0.600393262686049;
[x0fM, v0fM] = paramsToPosVel(th0fM, p0, e0, lan0, incl0, aop0, mu_K);
rf = norm(x0fM);
ir = x0fM / rf;
i1 = v0fM / norm(v0fM);
i3 = ir - i1 * dot(ir, i1);
i3 = i3 / norm(i3, 2);
i2 = cross(i3, i1);
[xffM, vffM] = paramsToPosVel(0, rf, 0, 0, 0, 0, mu_K);
ir(3, 1) = 0;
ith = cross([0;0;1], ir);
vffM = vffM(2, 1) * ith;
%dv = vffM - v0fM;

Tm = 2856413.395 - max(tts);
thm = timeToArg(Tm, p0, e0, mu_K, pe0T);
[x0m, v0m] = paramsToPosVel(thm, p0, e0, lan0, incl0, aop0, mu_K);

[Ksph_X, Ksph_Y, Ksph_Z] = genSphere(R_K, 51);

%% - Manoeuvre design -----------------------------------------------------

%{
%[tm] = argToTime(th0fM, p0, e0, mu_K, pe0T);
%[pf, ef, lanf, inclf, aopf, T_apf] = manoeuvre(tm, dv, p0, e0, lan0, incl0, aop0, mu_K, pe0T);

%Dvs(1, (end - 1):end) = 0.5 * Dvs(1, (end - 1):end);
Tm = 2856413.395;
if exist('intResPI.mat', 'file')
    fprintf('Loading file...\n');
    load intResPI.mat
else
    av_q = [ones(1, 31), 0.5 * ones(1, 31), 0.5, 0.6]';
end
tms = Tm - max(tts) + tts;
%vms = Dvs;
tt_q = 0;

n1 = 0;
n2 = 1250;
cnt = 0;
cnt2 = 0;
cm = cost(av_q, tms, Dvs, p0, e0, lan0, incl0, aop0, pe0T, mu_K, x0fM, vffM);
av_qm = av_q;
%av_qb = av_q;

while (cnt2 < 25)

    av_q = av_qm + (2 ^ -n1) * (rand(64, 1) - 0.5);
    av_q(av_q > 1) = 1;
    av_q(av_q < 0) = 0;

    [c] = cost(av_q, tms, Dvs, p0, e0, lan0, incl0, aop0, pe0T, mu_K, x0fM, vffM);

    if c < cm
       cm = c;
       av_qm = av_q;
       fprintf('n1 : %g; c = %23.16e\n', n1, cm);
       cnt = 0;
    end

    cnt = cnt + 1;

    if cnt >= n2
       n1 = n1 + 1;
       %av_qb = av_qm;
       cnt = 0;
       fprintf('n1 : %g; c = %23.16e\n', n1, cm);
    end

    if n1 > 20
        n1 = 1;
        cnt2 = cnt2 + 1;
        fprintf('n1 : %g; c = %23.16e\n', n1, cm);
    end

end

save intResPI.mat av_q

%tic
vms = Dvs;
    vms(1, (end - 1):end) = av_q(64, 1) * vms(1, (end - 1):end);
av = 2 * [av_q(1:31, 1)'; av_q(32:62, 1)'] * pi - pi;
tt = av_q(63, 1) * 7200 - 3600;
[th1s, th2s, ps, es, lans, incls, aops, peTs, dvs] = multiMan(tms + tt, vms, av(1, :), av(2, :), p0, e0, lan0, incl0, aop0, pe0T, mu_K);
[xf, vf] = paramsToPosVel(th1s(1, end), ps(1, end), es(1, end), lans(1, end), incls(1, end), aops(1, end), mu_K);
err = ((norm((xf - xffM), 2) / norm(xffM, 2)) + (norm((vf - vffM), 2) / norm(vffM, 2)));
%toc
%}

%% ------------------------------------------------------------------------

T = refFrameLocal(x0m, v0m);
t = Tm;
Dval = @(al, T, i) T' * (Dvs(1, i) * T * [cos(al); 0; sin(al)]);

al = 0;
viM = v0m + Dval(al, T, 1);
[thi, pi, ei, lani, incli, aopi] = posVelToParams(x0m, viM, mu_K);

thf0 = thi - thm + th0fM;


%% - Graphic output -------------------------------------------------------

figure(1)
plot3(x0(1, :), x0(2, :), x0(3, :), 'b', 'LineWidth', 2)
hold on
plot3(x0m(1, :), x0m(2, :), x0m(3, :), 'bo')
plot3(x0fM(1, :), x0fM(2, :), x0fM(3, :), 'bo')

%plot3(xf(1, :), xf(2, :), xf(3, :), 'r', 'LineWidth', 2)
surf(Ksph_X, Ksph_Y, Ksph_Z, zeros(size(Ksph_X)));
%{
for idx_1 = 1:size(th1s, 2)
    [xx, ~] = paramsToPosVel([th1s(1, idx_1), th2s(1, idx_1)], ps(1, idx_1), es(1, idx_1), lans(1, idx_1), incls(1, idx_1), aops(1, idx_1), mu_K);
    plot3(xx(1, :), xx(2, :), xx(3, :), 'ko')
    %if idx_1 < size(th1s, 2)
    th_c = linspace(th1s(1, idx_1), th2s(1, idx_1), 501);
    %else

    %end
    [xx, ~] = paramsToPosVel(th_c, ps(1, idx_1), es(1, idx_1), lans(1, idx_1), incls(1, idx_1), aops(1, idx_1), mu_K);
    plot3(xx(1, :), xx(2, :), xx(3, :), 'k')
end
%}
axis equal
hold off

%% - "contains": ----------------------------------------------------------

% Direction vector given two angles
function [v] = vdir(th, phi)

    v = zeros(3, size(th, 2));
    v(1, :) = cos(th).* cos(phi);
    v(2, :) = sin(th).* cos(phi);
    v(3, :) = sin(phi);

end

% Applies multiple manoeuvres and returns all intermediate states
function [th1s, th2s, ps, es, lans, incls, aops, peTs, dvs] = multiMan(tms, vms, a1s, a2s, p0, e0, lan0, incl0, aop0, peT0, mu)

    % Number of manoeuvres considered
    Nm    = size(tms, 2);

    % Output initialization
    th1s  = zeros(1, Nm);
    th2s  = zeros(1, Nm);
    ps    = zeros(1, Nm);
    es    = zeros(1, Nm);
    lans  = zeros(1, Nm);
    incls = zeros(1, Nm);
    aops  = zeros(1, Nm);
    peTs  = zeros(1, Nm);
    dvs   = zeros(3, Nm);

    % Initial values
    pi    = p0;
    ei    = e0;
    lani  = lan0;
    incli = incl0;
    aopi  = aop0;
    peTi  = peT0;

    for idx_1 = 1:Nm

        % Gets the argument of the manoeuvre
        [th_m] = timeToArg(tms(1, idx_1), pi, ei, mu, peTi);

        % (this is the last argument for the current orbit)
        if idx_1 > 1
            th2s(1, (idx_1 - 1)) = th_m;
        end

        % Gets position and velocity at the manoeuvre
        [xm, vm] = paramsToPosVel(th_m, pi, ei, lani, incli, aopi, mu);

        % Gets the frame of reference for the current point
        T = refFrameLocal(xm, vm);

        % Builds the velocity jump requested
        v_uv = vdir(a1s(1, idx_1), a2s(1, idx_1));
        v_uv = T * v_uv;

        % Builds the velocity
        dv = vms(1, idx_1) * v_uv;

        % Applies the manoeuvre
        [pf, ef, lanf, inclf, aopf, peTf] = manoeuvre(tms(1, idx_1), dv, pi, ei, lani, incli, aopi, mu, peTi);

        % Velocity jump
        dvs(:, idx_1) = T' * dv;

        % Gets the argument after the manoeuvre
        [th_f] = timeToArg(tms(1, idx_1), pf, ef, mu, peTf);
        th1s(1, idx_1) = th_f;

        % (Re-)assigns data
        pi = pf;
        ei = ef;
        lani = lanf;
        incli = inclf;
        aopi = aopf;
        peTi = peTf;

        ps(1, idx_1)    = pi;
        es(1, idx_1)    = ei;
        lans(1, idx_1)  = lani;
        incls(1, idx_1) = incli;
        aops(1, idx_1)  = aopi;
        peTs(1, idx_1)  = peTi;

    end

end

function [c] = cost(av_q, tms, Dvs, p0, e0, lan0, incl0, aop0, pe0T, mu_K, xffM, vffM)

    vms = Dvs;
    vms(1, (end - 1):end) = av_q(64, 1) * vms(1, (end - 1):end);
    av = 2 * [av_q(1:31, 1)'; av_q(32:62, 1)'] * pi - pi;
    tt = av_q(63, 1) * 7200 - 3600;
    [th1s, th2s, ps, es, lans, incls, aops, peTs, dvs] = multiMan(tms + tt, vms, av(1, :), av(2, :), p0, e0, lan0, incl0, aop0, pe0T, mu_K);
    mh = zeros(size(ps, 2), 1);
    for idx_1 = 1:(size(ps, 2) - 1)
        [xf, ~] = paramsToPosVel([th1s(1, idx_1), th2s(1, idx_1)], ps(1, idx_1), es(1, idx_1), lans(1, idx_1), incls(1, idx_1), aops(1, idx_1), mu_K);
        mh(idx_1, 1) = min(sqrt(sum((xf.* xf), 1)));
    end
    mh(end, 1) = ps(1, end) / (es(1, end) + 1);
    mhm = min(mh) - 600000;
    aemh = -0.000222318560702873;
    [xf, vf] = paramsToPosVel(th1s(1, end), ps(1, end), es(1, end), lans(1, end), incls(1, end), aops(1, end), mu_K);
    c = 100000 * (exp((norm((xf - xffM), 2) / norm(xffM, 2)) + (norm((vf - vffM), 2) / norm(vffM, 2))) - 1) / (exp(1) - 1) + av_q(64, 1) + 100000 * exp(aemh * mhm);

end

