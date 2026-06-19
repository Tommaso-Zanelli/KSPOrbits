fclose all;
close all
format long g
more off
clear
clc

load('Bashir_2_OrbitData.mat');
load('Bashir_2_TargetsData.mat');

SOIr_K = 84159286;
N = 251;

% Initial time
t0 = 3000000;

% Initial argument
thiK = timeToArg(t0, oOK.p, oOK.e, mu_K, oOK.peT);

% Exit argument
thOM = -thEM;

% Exit time
oT   = argToTime(thOM, oOM.p, oOM.e, mu_M, oOM.peT);

% Minmus position and velocity at exit
[xxOMK, vvOMK] = paramsToPosVel(timeToArg(oT, oMK.p, oMK.e, mu_K, oMK.peT), ...
    oMK.p, oMK.e, oMK.lan, oMK.incl, oMK.aop, mu_K);

% Object position and velocity relative to Minmus at exit
[xxOOM, vvOOM] = paramsToPosVel(thOM, oOM.p, oOM.e, oOM.lan, oOM.incl, oOM.aop, mu_M);

% Object position and velocity relative to Kerbin at exit
xxOOK = xxOOM + xxOMK;
vvOOK = vvOOM + vvOMK;

% New orbit around Kerbin after exit
[thOK, oOK2.p, oOK2.e, oOK2.lan, oOK2.incl, oOK2.aop] = posVelToParams(xxOOK, vvOOK, mu_K);

% Periapsis time for the new orbit
oOK2.peT = oT - argToTime(thOK, oOK2.p, oOK2.e, mu_K, 0);

% Kerbin exit argument for new orbit
thFK2 = acos((oOK2.p / SOIr_K - 1) / (oOK2.e));

% Kerbin exit time for new orbit
fT = argToTime(thFK2, oOK2.p, oOK2.e, mu_K, oOK2.peT);

%
eTi = eT - 10 * 3600;
oTf = oT + 10 * 3600;

% Builds Bashir 2 default trajectory (through Minmus...)
thB1 = linspace(thiK, thEK,  N);
thB2 = linspace(thEM, thOM,  N);
thM2B = linspace(timeToArg(eT, oMK.p, oMK.e, mu_K, oMK.peT), timeToArg(oT, oMK.p, oMK.e, mu_K, oMK.peT),  N);
thM2M = linspace(timeToArg(eTi, oMK.p, oMK.e, mu_K, oMK.peT), timeToArg(oTf, oMK.p, oMK.e, mu_K, oMK.peT),  N);
thB3 = linspace(thOK, thFK2, N);
[XXB1, ~] = paramsToPosVel(thB1, oOK.p, oOK.e, oOK.lan, oOK.incl, oOK.aop, mu_K);
[XXB2M, ~] = paramsToPosVel(thB2, oOM.p, oOM.e, oOM.lan, oOM.incl, oOM.aop, mu_M);
[XXM2K, ~] = paramsToPosVel(thM2B, oMK.p, oMK.e, oMK.lan, oMK.incl, oMK.aop, mu_K);
[XXM3K, ~] = paramsToPosVel(thM2M, oMK.p, oMK.e, oMK.lan, oMK.incl, oMK.aop, mu_K);
XXB2 = XXB2M + XXM2K;
[XXB3, ~] = paramsToPosVel(thB3, oOK2.p, oOK2.e, oOK2.lan, oOK2.incl, oOK2.aop, mu_K);

% Target orbits around Minmus
thT1 = linspace(timeToArg(eTi, oTM(1).p, oTM(1).e, mu_K, oTM(1).peT), timeToArg(oTf, oTM(1).p, oTM(1).e, mu_K, oTM(1).peT),  N);
[XXT1M, ~] = paramsToPosVel(thT1, oTM(1).p, oTM(1).e, oTM(1).lan, oTM(1).incl, oTM(1).aop, mu_M);
XXT1K = XXT1M + XXM3K;
thT2 = linspace(timeToArg(eTi, oTM(2).p, oTM(2).e, mu_K, oTM(2).peT), timeToArg(oTf, oTM(2).p, oTM(2).e, mu_K, oTM(2).peT),  N);
[XXT2M, ~] = paramsToPosVel(thT2, oTM(2).p, oTM(2).e, oTM(2).lan, oTM(2).incl, oTM(2).aop, mu_M);
XXT2K = XXT2M + XXM3K;
thT3 = linspace(timeToArg(eTi, oTM(3).p, oTM(3).e, mu_K, oTM(3).peT), timeToArg(oTf, oTM(3).p, oTM(3).e, mu_K, oTM(3).peT),  N);
[XXT3M, ~] = paramsToPosVel(thT3, oTM(3).p, oTM(3).e, oTM(3).lan, oTM(3).incl, oTM(3).aop, mu_M);
XXT3K = XXT3M + XXM3K;
thT4 = linspace(timeToArg(eTi, oTM(4).p, oTM(4).e, mu_K, oTM(4).peT), timeToArg(oTf, oTM(4).p, oTM(4).e, mu_K, oTM(4).peT),  N);
[XXT4M, ~] = paramsToPosVel(thT4, oTM(4).p, oTM(4).e, oTM(4).lan, oTM(4).incl, oTM(4).aop, mu_M);
XXT4K = XXT4M + XXM3K;




%% ------------------------------------------------------------------------
%  Graphic output

figure(1)
hold on
axis equal
plot3(0, 0, 0, 'bo', 'LineWidth', 2)
plot3(XXM3K(1, :), XXM3K(2, :), XXM3K(3, :), 'color', [0, 1, 0], 'LineWidth', 2)
plot3(XXB1(1, :), XXB1(2, :), XXB1(3, :), 'color', [0.75, 0, 0])
plot3(XXB1(1, 1), XXB1(2, 1), XXB1(3, 1), '*', 'color', [0.75, 0, 0])
plot3(XXB1(1, end), XXB1(2, end), XXB1(3, end), '*', 'color', [0.75, 0, 0])
plot3(XXB2(1, :), XXB2(2, :), XXB2(3, :), 'color', [1, 0, 0])
plot3(XXB2(1, 1), XXB2(2, 1), XXB2(3, 1), 'o', 'color', [1, 0, 0])
plot3(XXB2(1, end), XXB2(2, end), XXB2(3, end), 'o', 'color', [1, 0, 0])
plot3(XXB3(1, :), XXB3(2, :), XXB3(3, :), 'color', [1, 0.25, 0.25])
plot3(XXB3(1, 1), XXB3(2, 1), XXB3(3, 1), '*', 'color', [1, 0.25, 0.25])
plot3(XXB3(1, end), XXB3(2, end), XXB3(3, end), '*', 'color', [1, 0.25, 0.25])

plot3(XXT1K(1, :), XXT1K(2, :), XXT1K(3, :), 'color', [1, 1, 0])
plot3(XXT2K(1, :), XXT2K(2, :), XXT2K(3, :), 'color', [1, 0, 1])
plot3(XXT3K(1, :), XXT3K(2, :), XXT3K(3, :), 'color', [0, 1, 1])
plot3(XXT4K(1, :), XXT4K(2, :), XXT4K(3, :), 'color', [1, 0.5, 0])
hold off
