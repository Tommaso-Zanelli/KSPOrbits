fclose all;
close all
format long g
more off
clear
clc

vstr = version;
vnum = str2double(vstr(1, (end - 5):(end - 2)));

g12 = false;

R   = 60000;
d   = 2097152; %(40000 + R) * 9;

th1 = -pi / 6;
ph1 = pi * 0.25;

i1 = -[cos(th1) * cos(ph1); sin(th1) * cos(ph1); sin(ph1)];
x = -d * i1;

q = sqrt(d ^ 2 - R ^ 2);
i2 = [1; 0; 0];
if abs(abs(dot(i2, i1)) - 1) < 1e-12
    i2 = [0; 1; 0];
end
while abs(dot(i2, i1)) >= 1e-12
    i2 = i2 - i1 * dot(i2, i1);
    i2 = i2 / norm(i2, 2);
end
i3 = cross(i1, i2);
i3 = i3 / norm(i3, 2);
alpha = atan2(R, q);

c = @(beta) x + q * (i1 * cos(alpha) + (i2 * sin(beta) + i3 * cos(beta)) * sin(alpha));
ii = @(alp, beta) (i1 * cos(alp) + (i2 * sin(beta) + i3 * cos(beta)) * sin(alp));
qq = @(alp, beta) -dot(x, ii(alp, beta)) - real(sqrt((dot(x, ii(alp, beta)) ^ 2) - sum((x.^ 2), 1) + (R ^ 2)));
c2 = @(alp, beta) x + qq(alp, beta) * ii(alp, beta);


th = linspace(-pi, pi, 5001);
th3 = linspace(-pi, pi, 101);
al3 = alpha * (sqrt(linspace(0, 1, 101)));
cc = zeros(3, size(th, 2));
cc2 = zeros(3, size(th3, 2) * size(al3, 2));
ths = zeros(size(th));
phs = zeros(size(th));
th2 = zeros(1, size(th3, 2) * size(al3, 2));
ph2 = zeros(1, size(th3, 2) * size(al3, 2));

for idx_1 = 1:size(cc, 2)
    cc(:, idx_1) = c(th(1, idx_1));
    [ths(1, idx_1), phs(1, idx_1)] = polarCoords(cc(:, idx_1));
end
for idx_1 = 1:size(th3, 2)
    for idx_2 = 1:size(al3, 2)
        idx_3 = size(al3, 2) * (idx_1 - 1) + idx_2;
        cc2(:, idx_3) = c2(al3(1, idx_2), th3(1, idx_1));
        [th2(1, idx_3), ph2(1, idx_3)] = polarCoords(cc2(:, idx_3));
    end
end

[MX, MY, MZ] = genSphere(R, 51);

if g12

    figure(1)
    hold on
    axis equal
    if vnum >= 2019
        addToolbarExplorationButtons(gcf)
    end
    plot3([0, x(1, 1)], [0, x(2, 1)], [0, x(3, 1)], 'r', 'LineWidth', 2)
    plot3(x(1, 1), x(2, 1), x(3, 1), 'ro', 'LineWidth', 2)
    plot3(x(1, 1), x(2, 1), x(3, 1), 'k*', 'LineWidth', 2)
    plot3(cc(1, :), cc(2, :), cc(3, :), 'b', 'LineWidth', 3)
    plot3(cc2(1, :), cc2(2, :), cc2(3, :), 'c.', 'LineWidth', 1)
    surf(MX, MY, MZ, zeros(size(MX)));
    hold off

    figure(2)
    hold on
    if vnum >= 2019
        addToolbarExplorationButtons(gcf)
    end
    axis([-180, 180, -100, 100])
    plot([-180, 180], [-90, -90], 'k--')
    plot([-180, 180], [90, 90], 'k--')
    plot(180 * th2 / pi, 180 * ph2 / pi, 'c.')
    plot(180 * th1 / pi, 180 * ph1 / pi, 'ro', 'LineWidth', 3)
    plot(180 * th1 / pi, 180 * ph1 / pi, 'k*')
    plot(180 * ths / pi, 180 * phs / pi, 'b.')
    hold off

end

%% --

lan = 0;
incl = pi / 6;
aop = 0;

l = zeros(1, 4);
a = 0:0.5:1.5;
in = zeros(1, 4);

l = [0. 0.5 1.0 1.5];
a = [0. 1. 0. 1.];
in = [0, 0, 0, 0] / 3;

[x1, ~] = paramsToPosVel(th, d, 0, lan + l(1, 1) * pi, incl + in(1, 1) * pi, aop + a(1, 1) * pi, 0);
[x2, ~] = paramsToPosVel(th, d, 0, lan + l(1, 2) * pi, incl + in(1, 1) * pi, aop + a(1, 2) * pi, 0);
[x3, ~] = paramsToPosVel(th, d, 0, lan + l(1, 3) * pi, incl + in(1, 1) * pi, aop + a(1, 3) * pi, 0);
[x4, ~] = paramsToPosVel(th, d, 0, lan + l(1, 4) * pi, incl + in(1, 1) * pi, aop + a(1, 4) * pi, 0);

for idx_1 = 1:4
    oTM(idx_1).p    = d;
    oTM(idx_1).e    = 0;
    oTM(idx_1).peT  = 0;
    oTM(idx_1).lan  = lan  +  l(1, idx_1) * pi;
    oTM(idx_1).incl = incl + in(1, idx_1) * pi;
    oTM(idx_1).aop  = aop  +  a(1, idx_1) * pi;
end

TTM = 2 * argToTime(pi, d, 0, 1.7658e9, 0);

save Bashir_2_TargetsData.mat oTM TTM

[theta1, phi1, ~] = polarCoords(x1);
[theta2, phi2, ~] = polarCoords(x2);
[theta3, phi3, ~] = polarCoords(x3);
[theta4, phi4, ~] = polarCoords(x4);

figure(3)
hold on
if vnum >= 2019
    addToolbarExplorationButtons(gcf)
end
axis([-180, 180, -100, 100])
plot([-180, 180], [-90, -90], 'k--')
plot([-180, 180], [90, 90], 'k--')

[~, ~, th11, ph11, th21, ph21] = buildCone(x1(:, 1), R, 4096, 64, 64);
h31 = plot(180 * th11 / pi, 180 * ph11 / pi, '.', 'color', [1, 0.25, 0.25]);
h41 = plot(180 * th21 / pi, 180 * ph21 / pi, '.', 'color', [1, 0.75, 0.75]);

[~, ~, th11, ph11, th21, ph21] = buildCone(x2(:, 1), R, 4096, 64, 64);
h32 = plot(180 * th11 / pi, 180 * ph11 / pi, '.', 'color', [0.25, 1, 0.25]);
h42 = plot(180 * th21 / pi, 180 * ph21 / pi, '.', 'color', [0.75, 1, 0.75]);

[~, ~, th11, ph11, th21, ph21] = buildCone(x3(:, 1), R, 4096, 64, 64);
h33 = plot(180 * th11 / pi, 180 * ph11 / pi, '.', 'color', [0.25, 0.25, 1]);
h43 = plot(180 * th21 / pi, 180 * ph21 / pi, '.', 'color', [0.75, 0.75, 1]);

[~, ~, th11, ph11, th21, ph21] = buildCone(x4(:, 1), R, 4096, 64, 64);
h34 = plot(180 * th11 / pi, 180 * ph11 / pi, '.', 'color', [1, 1, 0.25]);
h44 = plot(180 * th21 / pi, 180 * ph21 / pi, '.', 'color', [1, 1, 0.75]);

plot(180 * theta1 / pi, 180 * phi1 / pi, 'k--')
plot(180 * theta2 / pi, 180 * phi2 / pi, 'k--')
plot(180 * theta3 / pi, 180 * phi3 / pi, 'k--')
plot(180 * theta4 / pi, 180 * phi4 / pi, 'k--')

h11 = plot(180 * theta1(1, 1) / pi, 180 * phi1(1, 1) / pi, 'ro');
h12 = plot(180 * theta2(1, 1) / pi, 180 * phi2(1, 1) / pi, 'go');
h13 = plot(180 * theta3(1, 1) / pi, 180 * phi3(1, 1) / pi, 'bo');
h14 = plot(180 * theta4(1, 1) / pi, 180 * phi4(1, 1) / pi, 'yo');
h21 = plot(180 * theta1(1, 1) / pi, 180 * phi1(1, 1) / pi, '*', 'color', [0.5, 0, 0]);
h22 = plot(180 * theta2(1, 1) / pi, 180 * phi2(1, 1) / pi, '*', 'color', [0, 0.5, 0]);
h23 = plot(180 * theta3(1, 1) / pi, 180 * phi3(1, 1) / pi, '*', 'color', [0, 0, 0.5]);
h24 = plot(180 * theta4(1, 1) / pi, 180 * phi4(1, 1) / pi, '*', 'color', [0.5, 0.5, 0]);

hold off

figure(4)
hold on
axis equal
if vnum >= 2019
    addToolbarExplorationButtons(gcf)
end
hsf = surf(MX, MY, MZ, zeros(size(MX)));
hsf.EdgeColor = 'none';
plot3(x1(1, :), x1(2, :), x1(3, :), 'r')
plot3(x2(1, :), x2(2, :), x2(3, :), 'g')
plot3(x3(1, :), x3(2, :), x3(3, :), 'b')
plot3(x4(1, :), x4(2, :), x4(3, :), 'y')

f11 = plot3(x1(1, 1), x1(2, 1), x1(3, 1), 'ro');
f12 = plot3(x2(1, 1), x2(2, 1), x2(3, 1), 'go');
f13 = plot3(x3(1, 1), x3(2, 1), x3(3, 1), 'bo');
f14 = plot3(x4(1, 1), x4(2, 1), x4(3, 1), 'yo');
f21 = plot3(x1(1, 1), x1(2, 1), x1(3, 1), 'r*');
f22 = plot3(x2(1, 1), x2(2, 1), x2(3, 1), 'g*');
f23 = plot3(x3(1, 1), x3(2, 1), x3(3, 1), 'b*');
f24 = plot3(x4(1, 1), x4(2, 1), x4(3, 1), 'y*');

[p1, p2] = buildCone(x1(:, 1), R, 4096, 64, 64);
f31 = plot3(p1(1, :), p1(2, :), p1(3, :), '.', 'color', [1, 0.25, 0.25]);
f41 = plot3(p2(1, :), p2(2, :), p2(3, :), '.', 'color', [1, 0.75, 0.75]);

[p1, p2] = buildCone(x2(:, 1), R, 4096, 64, 64);
f32 = plot3(p1(1, :), p1(2, :), p1(3, :), '.', 'color', [0.25, 1, 0.25]);
f42 = plot3(p2(1, :), p2(2, :), p2(3, :), '.', 'color', [0.75, 1, 0.75]);

[p1, p2] = buildCone(x3(:, 1), R, 4096, 64, 64);
f33 = plot3(p1(1, :), p1(2, :), p1(3, :), '.', 'color', [0.25, 0.25, 1]);
f43 = plot3(p2(1, :), p2(2, :), p2(3, :), '.', 'color', [0.75, 0.75, 1]);

[p1, p2] = buildCone(x4(:, 1), R, 4096, 64, 64);
f34 = plot3(p1(1, :), p1(2, :), p1(3, :), '.', 'color', [1, 1, 0.25]);
f44 = plot3(p2(1, :), p2(2, :), p2(3, :), '.', 'color', [1, 1, 0.75]);

hold off

%pause

if 1

    for idx_1 = [2:floor(size(th, 2) / 200):size(th, 2), 1]
        h11.XData = 180 * theta1(1, idx_1) / pi;
        h11.YData = 180 * phi1(1, idx_1) / pi;
        h12.XData = 180 * theta2(1, idx_1) / pi;
        h12.YData = 180 * phi2(1, idx_1) / pi;
        h13.XData = 180 * theta3(1, idx_1) / pi;
        h13.YData = 180 * phi3(1, idx_1) / pi;
        h14.XData = 180 * theta4(1, idx_1) / pi;
        h14.YData = 180 * phi4(1, idx_1) / pi;
        h21.XData = 180 * theta1(1, idx_1) / pi;
        h21.YData = 180 * phi1(1, idx_1) / pi;
        h22.XData = 180 * theta2(1, idx_1) / pi;
        h22.YData = 180 * phi2(1, idx_1) / pi;
        h23.XData = 180 * theta3(1, idx_1) / pi;
        h23.YData = 180 * phi3(1, idx_1) / pi;
        h24.XData = 180 * theta4(1, idx_1) / pi;
        h24.YData = 180 * phi4(1, idx_1) / pi;

        [p1, p2, th11, ph11, th21, ph21] = buildCone(x1(:, idx_1), R, 4096, 64, 64);
        h31.XData = 180 * th11 / pi;
        h31.YData = 180 * ph11 / pi;
        h41.XData = 180 * th21 / pi;
        h41.YData = 180 * ph21 / pi;
        f31.XData = p1(1, :);
        f31.YData = p1(2, :);
        f31.ZData = p1(3, :);
        f41.XData = p2(1, :);
        f41.YData = p2(2, :);
        f41.ZData = p2(3, :);

        [p1, p2, th11, ph11, th21, ph21] = buildCone(x2(:, idx_1), R, 4096, 64, 64);
        h32.XData = 180 * th11 / pi;
        h32.YData = 180 * ph11 / pi;
        h42.XData = 180 * th21 / pi;
        h42.YData = 180 * ph21 / pi;
        f32.XData = p1(1, :);
        f32.YData = p1(2, :);
        f32.ZData = p1(3, :);
        f42.XData = p2(1, :);
        f42.YData = p2(2, :);
        f42.ZData = p2(3, :);

        [p1, p2, th11, ph11, th21, ph21] = buildCone(x3(:, idx_1), R, 4096, 64, 64);
        h33.XData = 180 * th11 / pi;
        h33.YData = 180 * ph11 / pi;
        h43.XData = 180 * th21 / pi;
        h43.YData = 180 * ph21 / pi;
        f33.XData = p1(1, :);
        f33.YData = p1(2, :);
        f33.ZData = p1(3, :);
        f43.XData = p2(1, :);
        f43.YData = p2(2, :);
        f43.ZData = p2(3, :);

        [p1, p2, th11, ph11, th21, ph21] = buildCone(x4(:, idx_1), R, 4096, 64, 64);
        h34.XData = 180 * th11 / pi;
        h34.YData = 180 * ph11 / pi;
        h44.XData = 180 * th21 / pi;
        h44.YData = 180 * ph21 / pi;
        f34.XData = p1(1, :);
        f34.YData = p1(2, :);
        f34.ZData = p1(3, :);
        f44.XData = p2(1, :);
        f44.YData = p2(2, :);
        f44.ZData = p2(3, :);

        f11.XData = x1(1, idx_1);
        f11.YData = x1(2, idx_1);
        f11.ZData = x1(3, idx_1);
        f21.XData = x1(1, idx_1);
        f21.YData = x1(2, idx_1);
        f21.ZData = x1(3, idx_1);

        f12.XData = x2(1, idx_1);
        f12.YData = x2(2, idx_1);
        f12.ZData = x2(3, idx_1);
        f22.XData = x2(1, idx_1);
        f22.YData = x2(2, idx_1);
        f22.ZData = x2(3, idx_1);

        f13.XData = x3(1, idx_1);
        f13.YData = x3(2, idx_1);
        f13.ZData = x3(3, idx_1);
        f23.XData = x3(1, idx_1);
        f23.YData = x3(2, idx_1);
        f23.ZData = x3(3, idx_1);

        f14.XData = x4(1, idx_1);
        f14.YData = x4(2, idx_1);
        f14.ZData = x4(3, idx_1);
        f24.XData = x4(1, idx_1);
        f24.YData = x4(2, idx_1);
        f24.ZData = x4(3, idx_1);

        pause
    end
end
