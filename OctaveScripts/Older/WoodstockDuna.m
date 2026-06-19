fclose all;
close all
format long g
more off
clear
clc

R_K     = 600000;
mu_K    = 3.5316e12;
soiR_Kc = 84159286;
R_S     = 261600000;
mu_S    = 1.1723328e18;

now = 2934540;
apt = 853;
pet = 1785;
T_r = 1864;

aph = 79.9332e3;
peh = 74.9380e3;
[o.p, o.e] = apoPeriToPe((aph + R_K), (peh + R_K));

e_r = 0.0037;

o.lan = 206.3 * 1.745329251994329625305e-2;
o.incl = 0    * 1.745329251994329625305e-2;
o.aop = 140.4 * 1.745329251994329625305e-2;

o.peT = pet + now;

oK.p = 13599840256;
oK.e = 0;
oK.lan = 0;
oK.incl = 0;
oK.aop = 0;

soiR_K = oK.p * (mu_K / mu_S) ^ 0.4;
fprintf('%23.16e\n', abs(soiR_K - soiR_Kc))

m.v1 = 1122.270;
m.v2 = 0;
m.v3 = 0;
m.t  = 2936063.29;

[xm, vm] = paramsToPosVel(timeToArg(m.t, o.p, o.e, mu_K, o.peT),  o.p,  ...
                          o.e,  o.lan,  o.incl,  o.aop, mu_K);
[rfl] = refFrameLocal(xm, vm);
dv    = rfl(:, 1) * m.v1 + rfl(:, 2) * m.v2 + rfl(:, 3) * m.v3;
vm2 = vm + dv;
[thm2, o2.p, o2.e, o2.lan, o2.incl, o2.aop] = posVelToParams(xm, vm2, mu_K);
tm2 = argToTime(thm2, o2.p, o2.e, mu_K, 0);
o2.peT = m.t - tm2;
T2 =  2 * argToTime(pi, o2.p, o2.e, mu_K, 0);
pet2 = o2.peT;
apt2 = o2.peT - 0.5 * T2;
while (pet2 < m.t)
    pet2 = pet2 + T2;
end
while (apt2 < m.t)
    apt2 = apt2 + T2;
end

[ap2, pe2] = peToApoPeri(o2.p, o2.e);
ap2 = ap2 - R_K;
pe2 = pe2 - R_K;
fprintf('Apoapsis    : %g [m], at: %g [s]\n', ap2, apt2)
fprintf('Periapsis   : %g [m], at: %g [s]\n', pe2, pet2)
fprintf('Eccentricity: %g\n', o2.e)
fprintf('LAN         : %g [deg]\n', o2.lan * 57.2957795130823210042763)
fprintf('INCL        : %g [deg]\n', o2.incl * 57.2957795130823210042763)
fprintf('AOP         : %g [deg]\n', o2.aop * 57.2957795130823210042763)
