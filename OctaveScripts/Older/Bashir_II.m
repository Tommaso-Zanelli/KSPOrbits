fclose all;
close all
format long g
more off
clear
clc

NSphere = 256;
vstr = version;
vnum = str2double(vstr(1, (end - 5):(end - 2)));

R_M     = 60000;
mu_M    = 1.7658e9;
R_M_soi = 2247428.4;

ap = -97428.6;
pe = -33337.3;
[o0.p, o0.e] = apoPeriToPe((ap + R_M), (pe + R_M));
o0.peT = dot([21600, 3600, 60, 1], ([3, 0, 21, 0] + [138, 4, 49, 47]));
o0.lan  = pi * 88.2  / 180;
o0.incl = pi * 65.1  / 180;
o0.aop  = pi * 279.4 / 180;

% --
o1.p    = 2097152; %900000;
o1.e    = 0.0;
o1.lan  = 0;
o1.incl = pi / 6;
o1.aop  = 0;
o1.peT  = o0.peT;

o2.p    = 2097152; %900000;
o2.e    = 0.0;
o2.lan  = 0.5 * pi;
o2.incl = pi / 6;
o2.aop  = pi;
o2.peT  = o0.peT;

o3.p    = 2097152; %900000;
o3.e    = 0.0;
o3.lan  = pi;
o3.incl = pi / 6;
o3.aop  = 0;
o3.peT  = o0.peT;

o4.p    = 2097152; %900000;
o4.e    = 0.0;
o4.lan  = -0.5 * pi;
o4.incl = pi / 6;
o4.aop  = pi;
o4.peT  = o0.peT;

ttm0 = 120;
dtm  = 300;

th_l     = acos((o0.p / (R_M_soi * 1) - 1) / o0.e);
th_v     = linspace(-th_l, th_l, 501);
th_v2    = linspace(-pi, pi, 501);
[x_v, ~] = paramsToPosVel(th_v, o0.p, o0.e, o0.lan, o0.incl, o0.aop, mu_M);
tse = argToTime(-th_l, o0.p, o0.e, mu_M, o0.peT);
t1 = tse + ttm0;
t2 = tse + ttm0 + 1 * dtm;
t3 = tse + ttm0 + 2 * dtm;
t4 = tse + ttm0 + 3 * dtm;
theta1 = timeToArg(t1, o0.p, o0.e, mu_M, o0.peT);
theta2 = timeToArg(t2, o0.p, o0.e, mu_M, o0.peT);
theta3 = timeToArg(t3, o0.p, o0.e, mu_M, o0.peT);
theta4 = timeToArg(t4, o0.p, o0.e, mu_M, o0.peT);
[x1m, ~] = paramsToPosVel(theta1, o0.p, o0.e, o0.lan, o0.incl, o0.aop, mu_M);
[x2m, ~] = paramsToPosVel(theta2, o0.p, o0.e, o0.lan, o0.incl, o0.aop, mu_M);
[x3m, ~] = paramsToPosVel(theta3, o0.p, o0.e, o0.lan, o0.incl, o0.aop, mu_M);
[x4m, ~] = paramsToPosVel(theta4, o0.p, o0.e, o0.lan, o0.incl, o0.aop, mu_M);

[x1, ~] = paramsToPosVel(th_v2, o1.p, o1.e, o1.lan, o1.incl, o1.aop, mu_M);
[x2, ~] = paramsToPosVel(th_v2, o2.p, o2.e, o2.lan, o2.incl, o2.aop, mu_M);
[x3, ~] = paramsToPosVel(th_v2, o3.p, o3.e, o3.lan, o3.incl, o3.aop, mu_M);
[x4, ~] = paramsToPosVel(th_v2, o4.p, o4.e, o4.lan, o4.incl, o4.aop, mu_M);

[MX, MY, MZ] = genSphere(R_M, NSphere);
 MC          = MinmusCMap(NSphere);

figure(1)
hold on
axis equal
if vnum >= 2019
	addToolbarExplorationButtons(gcf)
end
plot3(x_v(1, :), x_v(2, :), x_v(3, :), '-', 'color', [1, 0.5, 0], 'LineWidth', 2)
plot3(x1m(1, :), x1m(2, :), x1m(3, :), 'o', 'color', [0.5, 0.25, 0], 'LineWidth', 1)
plot3(x2m(1, :), x2m(2, :), x2m(3, :), 'o', 'color', [0.5, 0.25, 0], 'LineWidth', 1)
plot3(x3m(1, :), x3m(2, :), x3m(3, :), 'o', 'color', [0.5, 0.25, 0], 'LineWidth', 1)
plot3(x4m(1, :), x4m(2, :), x4m(3, :), 'o', 'color', [0.5, 0.25, 0], 'LineWidth', 1)

plot3(x1(1, :), x1(2, :), x1(3, :), 'r')
plot3(x2(1, :), x2(2, :), x2(3, :), 'g')
plot3(x3(1, :), x3(2, :), x3(3, :), 'b')
plot3(x4(1, :), x4(2, :), x4(3, :), 'y')



hsf = surf(MX, MY, MZ, MC);
hsf.EdgeColor = 'none';
hold off

