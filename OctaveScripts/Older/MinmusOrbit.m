fclose all;
close all
format long g
more off
clear
clc

%% ------------------------------------------------------------------------

N      = 5001;

mu_K   = 3.5316e12;
mu_M   = 1.7658e9;
R_K    = 600000;
R_M    = 60000;
SOIr_M = 2247428.4;
hours  = [9203400, 21600, 3600, 60, 1];

degToRad = 1.745329251994329618642e-2;
aq = degToRad / 10;

%% ------------------------------------------------------------------------

% Object - Kerbin
timesOK  = [  0,   0,   2,   1,  46
            0, 139,   0,  28,  29];
paramsOK = [-5.61761e7, 8.30646e4, 1.0249, 309, 10.7, 0.8];

[oOK.p, oOK.e] = apoPeriToPe((paramsOK(1, 1) + R_K), (paramsOK(1, 2) + R_K));
oOK.lan  = paramsOK(1, 4) * degToRad;
oOK.incl = paramsOK(1, 5) * degToRad;
oOK.aop  = paramsOK(1, 6) * degToRad;
oOK.peT  = dot(timesOK(2, :), hours) - dot(timesOK(1, :), hours);

% Object - Minmus
timesOM  = [  0,   2,   4,  40,   0
            0, 139,   0,  31,   6];
paramsOM = [-9.74286e4, -3.33373e4, 5.9532, 88.2, 65.1, 279.4];

[oOM.p, oOM.e] = apoPeriToPe((paramsOM(1, 1) + R_M), (paramsOM(1, 2) + R_M));
oOM.lan  = paramsOM(1, 4) * degToRad;
oOM.incl = paramsOM(1, 5) * degToRad;
oOM.aop  = paramsOM(1, 6) * degToRad;
oOM.peT  = dot(timesOM(2, :), hours) + dot(timesOM(1, :), hours);

% Minmus encounter
timeseT  = [  0, 139,   0,  30,  36
            0,   2,   3,  36,   0];
eTc    = dot(timeseT(2, :), hours) + dot(timeseT(1, :), hours);

% Minmus - Kerbin
paramsMK = [47000000, 47000000, 0, 78, 6, 38];

[oMK.p, oMK.e] = apoPeriToPe(paramsMK(1, 1), paramsMK(1, 2));
oMK.lan  = paramsMK(1, 4) * degToRad;
oMK.incl = paramsMK(1, 5) * degToRad;
oMK.aop  = paramsMK(1, 6) * degToRad;
oMK.peT  = -argToTime(pi, oMK.p, oMK.e, mu_K, 0) * (0.9 / pi);

%% ------------------------------------------------------------------------

aOMp = 10 ^ (floor(log10(oOM.p)) - 5);
aSOI = 10 ^ (floor(log10(SOIr_M)) - 5);
aOKp = 10 ^ (floor(log10(oOK.p)) - 5);
aMKp = 10 ^ (floor(log10(oMK.p)) - 5);

%tt = linspace((eTc - 8 * 7200), (eTc + 8 * 7200), 50001);

%[th_MK] = timeToArg(tt, oMK.p, oMK.e, mu_K, oMK.peT);
%[th_OK] = timeToArg(tt, oOK.p, oOK.e, mu_K, oOK.peT);

%[xxOK, vvOK] = paramsToPosVel(th_OK, oOK.p, oOK.e, oOK.lan, oOK.incl, oOK.aop, mu_K);
%[xxMK, vvMK] = paramsToPosVel(th_MK, oMK.p, oMK.e, oMK.lan, oMK.incl, oMK.aop, mu_K);

%{
NT = 51;
tt1 = linspace(-28.071195342766081637,    -28.071195322777601433, NT);
tt2 = linspace(-47.896026888601603844,    -47.896026868613120087, NT);
[TT1, TT2] = meshgrid(tt1, tt2);
dd = zeros(size(TT1));

% Encounter
%[thEM] = timeToArg(eT, oOM.p, oOM.e, mu_M, oOM.peT);
for idx_1 = 1:NT
    for idx_2 = 1:NT
%}

save MinmusOrbitData.mat

%a = zeros(12, 1);

        %a(6, 1) = -28.071195327375 / 61; %TT1(idx_1, idx_2);
        %a(12, 1) = -47.8960268728107 / 61; %TT2(idx_1, idx_2);

cf = @(xx)  MinmusOrbitCost01(xx);

x0 = [     -0.000832728822884588
       -0.0211873382153458
         0.269165299197424
      0.000179376563343825
         0.386899953248897
         0.152897223620252
      -0.00112928021251684
      0.000284520665238659
         0.193205626723198
        0.0736324471247531
         0.166010669510419
        -0.413775597164129];

%a = optiGrad(cf, x0, 1e-15, 1e-9, 15);
a = x0;

        thEM = -acos(((oOM.p + aOMp * a(1, 1)) / SOIr_M - 1) / (oOM.e + 0.0001 * a(2, 1)));
        eT   =  argToTime(thEM, (oOM.p + aOMp * a(1, 1)), (oOM.e + 0.0001 * a(2, 1)), ...
            mu_M, (oOM.peT + 61 * a(6, 1)));
        [xxEM, vvEM] = paramsToPosVel(thEM, (oOM.p + aOMp * a(1, 1)), ...
            (oOM.e + 0.0001 * a(2, 1)), (oOM.lan + aq * a(3, 1)), ...
            (oOM.incl + aq * a(4, 1)), (oOM.aop + aq * a(5, 1)), mu_M);


        [thEK] = timeToArg(eT, (oOK.p + aOKp * a(7, 1)), (oOK.e + 0.0001 * a(8, 1)), ...
            mu_K, (oOK.peT + 61 * a(12, 1)));
        [xxEK, vvEK] = paramsToPosVel(thEK, (oOK.p + aOKp * a(7, 1)), ...
            (oOK.e + 0.0001 * a(8, 1)), (oOK.lan + aq * a(9, 1)), ...
            (oOK.incl + aq * a(10, 1)), (oOK.aop + aq * a(11, 1)), mu_K);

        [thMK] = timeToArg(eT, oMK.p, oMK.e, mu_K, oMK.peT);
        [xxMK, vvMK] = paramsToPosVel(thMK, oMK.p, oMK.e, oMK.lan, ...
                                        oMK.incl, oMK.aop, mu_K);

    cc = ((norm(xxEK - xxMK, 2) - SOIr_M)) ^ 2 + ...
            (norm(((xxEK - xxMK) - xxEM), 2)) ^ 2  + ...
            (SOIr_M * (norm(((vvEK - vvMK) - vvEM), 2) / norm(vvMK, 2))) ^ 2;


%save MinmusOrbitData02.mat

%x0 = x0(7:12, 1);

%cf = @(xx)  MinmusOrbitCost02(xx);

%b = optiGrad(cf, x0, 1e-18, 1e-9, 15);


        %dd(idx_1, idx_2) = norm(xxEK - xxMK, 2) - SOIr_M;

%{
    end
end

[i1, i2] = find(abs(dd)==min(min(abs(dd))));
mv = dd(i1, i2);
t1i = [TT1(i1, i2-1), TT1(i1, i2+1)];
t2i = [TT2(i1-1, i2), TT2(i1+1, i2)];
fprintf('%g\n', mv);
fprintf('%25.18f, %25.18f\n', min(min(t1i)), max(max(t1i)));
fprintf('%25.18f, %25.18f\n', min(min(t2i)), max(max(t2i)));
%}

%MT = argToTime(pi, oMK.p, oMK.e, mu_K, 0);
%ff = @(tt) (norm((xxEK - paramsToPosVel(timeToArg(eT, oMK.p, oMK.e, mu_K, tt), oMK.p, oMK.e, oMK.lan, oMK.incl, oMK.aop, mu_K)), 2) - SOIr_M) ^ 2;

%tv = linspace(-6 * MT, 6 * MT, N);
%cv = zeros(1, N);
%for idx_1 = 1:N
%    cv(1, idx_1) = ff(tv(1, idx_1));
%end

%[~, pm] = min(cv);
%tv = linspace(tv(1, max([(pm(1, 1) - 1), 1])), tv(1, min([(pm(1, 1) + 1), N])), N);
%[x0, xs, xn] = optiFun([-MT, MT], ff, 1e-9, 2501, 1e-16)

figure(1)
plot3(0,0,0,'bo')
hold on
axis equal
plot3(xxMK(1, 1), xxMK(2, 1),xxMK(3, 1),'go')
m=max(abs(xxMK));
axis([-m, m, -m, m, -m, m])
%plot3(xxEM(1, 1)+xxMK(1, 1), xxEM(2, 1)+xxMK(2, 1),xxEM(3, 1)+xxMK(3, 1),'m*')
plot3(xxEK(1, 1), xxEK(2, 1),xxEK(3, 1),'ro')
hold off


%figure(2)
%hsrf2 = surf(TT1, TT2, dd);
%hsrf2.EdgeColor = 'none';

%{
plot3(0,0,0,'bo')
hold on
axis equal
m = 1.25 * oMK.p;
axis([-m, m, -m, m, -m, m])
plot3(xxOK(1, :), xxOK(2, :), xxOK(3, :), 'r')
plot3(xxOK(1, 1), xxOK(2, 1), xxOK(3, 1), 'ro')
plot3(xxOK(1, end), xxOK(2, end), xxOK(3, end), 'ro')
plot3(xxMK(1, :), xxMK(2, :), xxMK(3, :), 'g')
plot3(xxMK(1, 1), xxMK(2, 1), xxMK(3, 1), 'go')
plot3(xxMK(1, end), xxMK(2, end), xxMK(3, end), 'go')
hold off
%}