%% ------------------------------------------------------------------------

fclose all;
close all
format long g
more off
clear
clc

%% ------------------------------------------------------------------------

mu_1 = 1e12;
mu_2 = 1e9;

oA.p    = 1e8;
oA.e    = 0;
oA.lan  = 0;
oA.incl = 0;
oA.aop  = 0;
oA.peT  = 0;

%% ------------------------------------------------------------------------

rSOI = (oA.p / (1 - (oA.e ^ 2))) * ((mu_2 / mu_1) ^ 0.4);

o2.p    = 2570279.468258769251406192779541015625;
o2.e    = 2.8784772791787087697912284056656062602996826171875;
o2.lan  = 0.97854767267054254542557600871077738702297210693359375;
o2.incl = -2.91720955947772520033822729601524770259857177734375;
o2.aop  = 2.193644124940508088883461823570542037487030029296875;
o2.peT  = 238655.02265783547773025929927825927734375;

th_2e   = -acos((o2.p / rSOI - 1) / o2.e);
tt_e    = argStructToTime(th_2e, o2, mu_2);

[xxAe, vvAe, th_Ae] = structTimeToPosVelArg(tt_e, oA, mu_1);
[xx2e, vv2e, th_2e] = structTimeToPosVelArg(tt_e, o2, mu_2);
xx2 = xx2e + xxAe;
vv2 = vv2e + vvAe;
[th_1e, o1] = posVelToStruct(xx2, vv2, mu_1, tt_e);

%% ------------------------------------------------------------------------

rTd = 5.729577951308232286464772187173366546630859375e1;

[ap1, pe1] = peToApoPeri(o1.p, o1.e);
[ap2, pe2] = peToApoPeri(o2.p, o2.e);
lan1  = (o1.lan  * rTd);
incl1 = (o1.incl * rTd);
aop1  = (o1.aop  * rTd);
lan2  = (o2.lan  * rTd);
incl2 = (o2.incl * rTd);
aop2  = (o2.aop  * rTd);
[peT1, r_peT1] = secondsToYDHMSTime(o1.peT);
[peT2, r_peT2] = secondsToYDHMSTime(o2.peT);

f = @(x, n) floor(x / (10 ^ floor(log10(abs(x)) - (n - 1))) + 0.5) * (10 ^ floor(log10(abs(x)) - (n - 1)));

ap1 = f(ap1, 6);
pe1 = f(pe1, 6);
ap2 = f(ap2, 6);
pe2 = f(pe2, 6);

lan1  = floor(lan1  * 10 + 0.5) * 0.1;
incl1 = floor(incl1 * 10 + 0.5) * 0.1;
aop1  = floor(aop1  * 10 + 0.5) * 0.1;
lan2  = floor(lan2  * 10 + 0.5) * 0.1;
incl2 = floor(incl2 * 10 + 0.5) * 0.1;
aop2  = floor(aop2  * 10 + 0.5) * 0.1;

r_peT1 = 0;
r_peT2 = 0;

ref = [pe1, ap1, o1.e, lan1, incl1, aop1, 0; ...
       1, peT1, mu_1; ...
       1, 0, 0, 0, 0, r_peT1, 0; ...
       1, 0, 0, 0, 0, tt_e, mu_2; ...
       1, 0, 0, 0, 0, 0, rSOI; ...
       pe2, ap2, o2.e, lan2, incl2, aop2, 0; ...
       1, peT2, 0; ...
       1, 0, 0, 0, 0, r_peT2, 0];

x0 = zeros(12, 1);
fc = @(x) refineEncounter_cost(x, ref, oA, [], true);
x = optiGrad(fc, x0, 1e-6, 1e-8, 15);

[cc, o1c, o2c] = refineEncounter_cost(x, ref, oA, [], true);

fprintf('%g\n%g\n%g\n%g\n%g\n%g\n%g\n%g\n%g\n%g\n%g\n%g\n', ...
    abs((o1c.p    - o1.p)    / o1.p), ...
    abs((o1c.e    - o1.e)    / o1.e), ...
    abs((o1c.lan  - o1.lan)  / o1.lan), ...
    abs((o1c.incl - o1.incl) / o1.incl), ...
    abs((o1c.aop  - o1.aop)  / o1.aop), ...
    abs((o1c.peT  - o1.peT)  / o1.peT), ...
    abs((o2c.p    - o2.p)    / o2.p), ...
    abs((o2c.e    - o2.e)    / o2.e), ...
    abs((o2c.lan  - o2.lan)  / o2.lan), ...
    abs((o2c.incl - o2.incl) / o2.incl), ...
    abs((o2c.aop  - o2.aop)  / o2.aop), ...
    abs((o2c.peT  - o2.peT)  / o2.peT))
