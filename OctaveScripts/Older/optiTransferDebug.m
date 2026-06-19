%% -----------------------------------------------------------------------------

fclose all;
close all
format long g
more off
clear
clc

%% -----------------------------------------------------------------------------

mu_K = 3.5316e12;

% Orbit 1:
p1 = 1400000;
e1 = 0.39;
lan1 = 0.5;
incl1 = 3;
aop1 = 1.9;
th_p1 = 1.0;

% Orbit 2:
p2 = 700000;
e2 = 0.02;
lan2 = 0;
incl2 = -1;
aop2 = 0.3;
th_p2 = 0.9;

% First point
[x1, v1] = paramsToPosVel(th_p1, p1, e1, lan1, incl1, aop1, mu_K);
[xx1, ~] = paramsToPosVel(linspace(-pi, pi, 501), p1, e1, lan1, incl1, aop1, mu_K);

% Second point
[x2, v2] = paramsToPosVel(th_p2, p2, e2, lan2, incl2, aop2, mu_K);
[xx2, ~] = paramsToPosVel(linspace(-pi, pi, 501), p2, e2, lan2, incl2, aop2, mu_K);

k = cross(x1, x2);
k = k / norm(k, 2);
[lant, inclt, aop_tht] = rotAngles(x1, x2);
r1    = norm(x1, 2);
r2    = norm(x2, 2);
xi    = r1 / r2;
dth   = atan2(norm(cross(x2, x1)), dot(x1, x2));
phi_d = atan2(sin(dth), (cos(dth) - xi));
phi_p = atan2(sin(dth), (cos(dth) - 1));
[pp1, pn1, pp2, pn2] = posNegAngles(phi_d, phi_p);
if xi > 1
    thi = [pp1, pp2];
else
    thi = [pn1, pn2];
end
p = @(th0) (r2 * xi * (cos(dth + th0) - cos(th0)))./ (cos(dth + th0) - xi*cos(th0));
e = @(th0) (xi - 1)./ (cos(dth + th0) - xi*cos(th0));

th0v = linspace(thi(1, 1) + 0.04, thi(1, 2) - 0.04, 501);
figure(2)
hold on
plot(th0v, p(th0v) / r2, 'b')
plot(th0v, e(th0v), 'r')
hold off

th00 = mean(thi);

%{
l = 0;
if (abs(inclt) > 0.5 * pi)
    l = pi;%inclt;
end
l
%}

pt = p(th00);
et = e(th00);
[x1t, ~] = paramsToPosVel(th00, pt, et, lant, inclt, (aop_tht - th00), mu_K);
aopt = aop_tht - th00 + acos((dot(x1, x1t) / (norm(x1, 2) * norm(x1t, 2))));

[x1t, ~] = paramsToPosVel(th00, pt, et, lant, inclt, aopt, mu_K);
[x2t, ~] = paramsToPosVel((th00 + dth), pt, et, lant, inclt, aopt, mu_K);
[xxt, ~] = paramsToPosVel(linspace(th00, th00 + dth, 501), pt, et, lant, inclt, aopt, mu_K);


%{
% Defines the cost function
cf = @(v1t, v2t) norm((v1 - v1t), 2) + norm((v2 - v2t), 2);

% Looks for the optimal transfer orbit
[th_t1, th_t2, p_t, e_t, lan_t, incl_t, aop_t] = optiTransfer(x1, x2, cf, mu_K);


% Builds the transfer orbit
th_tt = linspace(th_t1, th_t2, 501);
[xxt, ~] = paramsToPosVel(th_tt, p_t, e_t, lan_t, incl_t, aop_t, mu_K);
[~, vtt] = paramsToPosVel([th_t1, th_t2], p_t, e_t, lan_t, incl_t, aop_t, mu_K);
fprintf('%g\n', cf(vtt(:, 1), vtt(:, 2)))
%}

% Graphic output
figure(1)
%addToolbarExplorationButtons(gcf)
hold on
plot3(0,0,0, 'sk')
plot3(xx1(1, :), xx1(2, :), xx1(3, :), 'b')
plot3(x1(1, 1), x1(2, 1), x1(3, 1), 'bo')
plot3(xx2(1, :), xx2(2, :), xx2(3, :), 'r')
plot3(x2(1, 1), x2(2, 1), x2(3, 1), 'ro')
%plot3(x1t(1, 1), x1t(2, 1), x1t(3, 1), 'k*')
%plot3(x2t(1, 1), x2t(2, 1), x2t(3, 1), 'k*')
%plot3(xxt(1, :), xxt(2, :), xxt(3, :), 'k')
%plot3(xxt(1, 1), xxt(2, 1), xxt(3, 1), 'k*')
%plot3(xxt(1, end), xxt(2, end), xxt(3, end), 'k*')

for th00 = linspace(-pi, pi, 144)
    %{
        pt = p(th00);
        et = e(th00);
        [x1t, ~] = paramsToPosVel(th00, pt, et, lant, inclt, (aop_tht - th00), mu_K);
        aopt = aop_tht - th00 + acos(min([max([-1, (dot(x1, x1t) / (norm(x1, 2) * norm(x1t, 2)))]), 1]));
acos(min([max([-1, (dot(x1, x1t) / (norm(x1, 2) * norm(x1t, 2)))]), 1]));

        [x1t, ~] = paramsToPosVel(th00, pt, et, lant, inclt, aopt, mu_K);
        [x2t, ~] = paramsToPosVel((th00 + dth), pt, et, lant, inclt, aopt, mu_K);
        [xxt, ~] = paramsToPosVel(linspace(th00, th00 + dth, 501), pt, et, lant, inclt, aopt, mu_K);
%}
[th_1, th_2, p, e, lan, incl, aop] = transferOrbit(x1, x2, th00);

if ~isnan(p)

    [x1t, ~] = paramsToPosVel(th_1, p, e, lan, incl, aop, mu_K);
        [x2t, ~] = paramsToPosVel(th_2, p, e, lan, incl, aop, mu_K);
        [xxt, ~] = paramsToPosVel(linspace(th_1, th_2, 501), p, e, lan, incl, aop, mu_K);


plot3(x1t(1, 1), x1t(2, 1), x1t(3, 1), 'k*')
plot3(x2t(1, 1), x2t(2, 1), x2t(3, 1), 'k*')
isclose = sqrt(sum((xxt.^ 2), 1)) < 2400000000;
plot3(xxt(1, isclose), xxt(2, isclose), xxt(3, isclose), 'k')


end
end

cf = @(v1t, v2t) dot((v1 - v1t), (v1 - v1t)) + dot((v2 - v2t), (v2 - v2t));
[th_1, th_2, p, e, lan, incl, aop] = optiTransfer2(x1, x2, cf, mu_K);
[xxto, ~] = paramsToPosVel(linspace(th_1, th_2, 251), p, e, lan, incl, aop, mu_K);
plot3(xxto(1, :), xxto(2, :), xxto(3, :), 'm', 'LineWidth', 2)

axis([-2350000, 2350000, -2350000, 2350000, -2350000, 2350000])

axis equal
hold off

