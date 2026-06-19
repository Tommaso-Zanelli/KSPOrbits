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
p1 = 840000;
e1 = 0;
lan1 = 0;
incl1 = 0;
aop1 = 0;

% Orbit 2:
p2 = 672500;
e2 = 0;
lan2 = 0;
incl2 = 0;
aop2 = 0;

% First point
[x1, v1] = paramsToPosVel((pi / 6), p1, e1, lan1, incl1, aop1, mu_K);
[xx1, ~] = paramsToPosVel(linspace(-pi, pi, 501), p1, e1, lan1, incl1, aop1, mu_K);

% Second point
[x2, v2] = paramsToPosVel((2 * pi / 3), p2, e2, lan2, incl2, aop2, mu_K);
[xx2, ~] = paramsToPosVel(linspace(-pi, pi, 501), p2, e2, lan2, incl2, aop2, mu_K);








% Graphic output
figure(1)
hold on
plot(xx1(1, :), xx1(2, :), 'b')
plot(x1(1, 1), x1(2, 1), 'bo')
plot(xx2(1, :), xx2(2, :), 'r')
plot(x2(1, 1), x2(2, 1), 'ro')
axis equal
hold off

