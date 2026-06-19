%%  ----  Header  -----------------------------------------------------

fclose all;
close all
format long g
more off
clear
clc

%%  ----

% Number of test points
N = 501;

for rho2exp = log10(0.9) %linspace(-5, 5, ceil(sqrt(N)))
for dths = 0.25 * pi %linspace(0, 2 * pi, ceil(sqrt(N)))

% Points to be connected
rho1 = 1;
%rho2 = 1.2e-3;
rho2 = 10 ^ rho2exp;
%dth = pi / 4;
dth = dths;

% Parameters
xi = rho1 / rho2;
m1 = sqrt((xi - cos(dth)) ^ 2 + (sin(dth)) ^ 2);
m2 = sqrt((1  - cos(dth)) ^ 2 + (sin(dth)) ^ 2);
q1 = (1 - xi) / m1;
q2 = rho1 * m2 / m1;

% Angles
phi1 = -atan2((sin(dth)), (xi - cos(dth)));
phi2 = -atan2((sin(dth)), (1  - cos(dth)));

gamma1 = atan2(sin(phi1), (q1 + cos(phi1)));
gamma2 = atan2((q1 * sin(dth) + sin(phi1)), (q1 * cos(dth) + cos(phi1)));

%fprintf('DeltaGamma: %g\n', abs(gamma1 - gamma2))

%l3 = ([-1, 1] + 2 * (q1 < 0)) * pi * 0.5 - gamma1;
l3 = -phi1 + [-1, 1] * acos(abs(q1)) + pi * (q1 < 0);
if l3(1, 1) < -pi
    l3 = [-pi, l3(1, 2); (l3(1, 1) + 2 * pi), pi];
elseif l3(1, 2) > pi
    l3 = [-pi, (l3(1, 2) - 2 * pi); l3(1, 1), pi];
end

alpha = atan2((q1 + cos(phi1)), sin(phi1));
beta  = atan2((q1 * cos(dth) + cos(phi1)), (q1 * sin(dth) + sin(phi1))) ;


% Limits
if q1 ~= 0
    l1 = ([-1, 1] + 2 * (q1 < 0)) * pi * 0.5 - phi1;
    o1 = [false, false];
    l2 = ([-1, 1] + 2 * (q1 < 0)) * pi * 0.5 - phi2;
    o2 = [true, true];
else
    l1 = [-pi, pi];
    l2 = [-pi, pi];
end
if phi1 < phi2
    l = [l1(1, 1), l2(1, 2)];
    o = [o1(1, 1), o2(1, 2)];
elseif phi1 > phi2
    l = [l2(1, 1), l1(1, 2)];
    o = [o2(1, 1), o1(1, 2)];
else
    l = [-pi, pi];
end
if l(1, 1) < -pi
    l = [-pi, l(1, 2); (l(1, 1) + 2 * pi), pi];
elseif l(1, 2) > pi
    l = [-pi, (l(1, 2) - 2 * pi); l(1, 1), pi];
end

% Initial arguments
th0 = linspace(-pi, pi, N);

isv = zeros(3, N);
pp  = zeros(1, N);
ee  = zeros(1, N);

%peerr = zeros(2, N);

for idx_1 = 1:N
    %M  = [1, -rho1 * cos(th0(1, idx_1)); 1, -rho2 * cos(th0(1, idx_1) + dth)];
    %pe = M \ [rho1; rho2];

    d    = cos(th0(1, idx_1) + phi1);

    p    = q2 * cos(th0(1, idx_1) + phi2) / d;
    e    = q1 / d;

    %peerr(1, idx_1) = abs((pe(1, 1) - p) / p);
    %peerr(2, idx_1) = abs(pe(2, 1) - e);

    pp(1, idx_1) = p;
    ee(1, idx_1) = e;

    isv(1, idx_1) = NaN;
    isv(2, idx_1) = NaN;
    if ((p > 0) && (e >= 0))

        isv(1, idx_1) = e;
        isv(2, idx_1) = p;
    end
       isv(3, idx_1) = 1;
        if abs(e) >= 1
            isv(3, idx_1) = 2;
            thl = acos(-1 / e);
            if ((th0(1, idx_1) < (-thl)) || ((th0(1, idx_1) + dth) > thl))
                isv(3, idx_1) = 3;
            end
        end


end

tiv_i1 = find((isv(3, 2:end) ~= 3 & isv(3, 2:end) > 0 & isv(3, 1:(end - 1)) == 3) | ...
    (isv(3, 2:end) == 3 & isv(3, 1:(end - 1)) ~= 3 & isv(3, 1:(end - 1)) > 0));
tiv_i2 = tiv_i1 + 1;
th_tiv = 0.5 * (th0(1, tiv_i1) + th0(1, tiv_i2));
pp_tiv = 0.5 * (pp(1, tiv_i1)  + pp(1, tiv_i2));
ee_tiv = 0.5 * (ee(1, tiv_i1)  + ee(1, tiv_i2));

if size(l, 1) == 1
if all(th0(~isnan(isv(1, :))) >= l(1, 1)) && all(th0(~isnan(isv(1, :))) <= l(1, 2))
    %fprintf('Ok!\n')
else
    %fprintf('Porcodio: rho2 = %g, dth = %g*pi\n', rho2, (dth / pi))
    fprintf('Porcodio: ')
    if ~all(th0(~isnan(isv(1, :))) >= l(1, 1))
        fprintf('by %g radians on the left ', abs(min(th0(~isnan(isv(1, :)))) - l(1, 1)))
    end
    if  ~all(th0(~isnan(isv(1, :))) <= l(1, 2))
        fprintf('by %g radians on the right ', abs(max(th0(~isnan(isv(1, :)))) - l(1, 2)))
    end
    fprintf('\n')
end
else
if all((th0(~isnan(isv(1, :))) >= l(1, 1) & th0(~isnan(isv(1, :))) <= l(1, 2)) | ...
    ((th0(~isnan(isv(1, :))) >= l(2, 1)) & (th0(~isnan(isv(1, :))) <= l(2, 2))))
    %fprintf('Ok!\n')
else
    %fprintf('Porcodio: rho2 = %g, dth = %g*pi\n', rho2, (dth / pi))
    fprintf('Porcodio: ')
    if ~all(th0(~isnan(isv(1, :))) >= l(1, 1)) && ~all(th0(~isnan(isv(1, :))) >= l(2, 1))
        fprintf('by %g radians on the left ', abs(min(th0(~isnan(isv(1, :)))) - l(1, 1)))
    end
    if  ~all(th0(~isnan(isv(1, :))) <= l(1, 2)) &&  ~all(th0(~isnan(isv(1, :))) <= l(2, 2))
        fprintf('by %g radians on the right ', abs(max(th0(~isnan(isv(1, :)))) - l(1, 2)))
    end
    fprintf('\n')
end
end

end
end

hf(1, 1) = figure(1);
plot(th0, ee, '--', 'color', [1, 0.75, 0.75])
hold on
plot([-pi, -pi], [min(ee), max(ee)], '--', 'color', [0.75, 0.75, 0.75])
plot([pi, pi], [min(ee), max(ee)], '--', 'color', [0.75, 0.75, 0.75])
plot(th0, isv(1, :), 'r', 'LineWidth', 3)
plot(th0(~isnan(isv(1, :))), ones(size(th0(~isnan(isv(1, :))))), '--', 'color', [0.5, 0.5, 0.5])
if size(l, 1) == 1
plot([l(1, 1), l(1, 1)], [min(ee), max(ee)], 'm', 'LineWidth', 1.5)
plot(l(1, 1), 1, 'mo')
plot([l(1, 2), l(1, 2)], [min(ee), max(ee)], 'c', 'LineWidth', 1.5)
plot(l(1, 2), 1, 'co')
else
plot([l(2, 1), l(2, 1)], [min(ee), max(ee)], 'm', 'LineWidth', 1.5)
plot(l(2, 1), 1, 'mo')
plot([l(1, 2), l(1, 2)], [min(ee), max(ee)], 'c', 'LineWidth', 1.5)
plot(l(1, 2), 1, 'co')
end
plot(l3(1, 1), 0, 's', 'color', [0.5, 0.25, 0], 'MarkerSize', 9)
plot(l3(1, 2), 0, 's','color', [0.25, 0, 0.5], 'MarkerSize', 9)
plot(th_tiv, ee_tiv, 's','color', [0.75, 1, 0.25], 'MarkerSize', 9)
hold off
addToolbarExplorationButtons(hf(1, 1))

hf(1, 2) = figure(2);
plot(th0, pp, '--', 'color', [0.75, 0.75, 1])
hold on
plot([-pi, -pi], [min(pp), max(pp)], '--', 'color', [0.75, 0.75, 0.75])
plot([pi, pi], [min(pp), max(pp)], '--', 'color', [0.75, 0.75, 0.75])
plot(th0, isv(2, :), 'b', 'LineWidth', 3)
if size(l, 1) == 1
plot([l(1, 1), l(1, 1)], [min(pp), max(pp)], 'm', 'LineWidth', 1.5)
plot(l(1, 1), 0, 'mo')
plot([l(1, 2), l(1, 2)], [min(pp), max(pp)], 'c', 'LineWidth', 1.5)
plot(l(1, 2), 0, 'co')
else
plot([l(2, 1), l(2, 1)], [min(pp), max(pp)], 'm', 'LineWidth', 1.5)
plot(l(2, 1), 0, 'mo')
plot([l(1, 2), l(1, 2)], [min(pp), max(pp)], 'c', 'LineWidth', 1.5)
plot(l(1, 2), 0, 'co')
end
plot(l3(1, 1), 0, 's', 'color', [0.5, 0.25, 0], 'MarkerSize', 9)
plot(l3(1, 2), 0, 's','color', [0.25, 0, 0.5], 'MarkerSize', 9)
plot(th_tiv, pp_tiv, 's','color', [0.75, 1, 0.25], 'MarkerSize', 9)
hold off
addToolbarExplorationButtons(hf(1, 2))

pause
figure(2)
hold on
plot(th0(1, (isv(3, :) == 1)), pp(1, (isv(3, :) == 1)), 'g*')
plot(th0(1, (isv(3, :) == 2)), pp(1, (isv(3, :) == 2)), 'y*')
plot(th0(1, (isv(3, :) == 3)), pp(1, (isv(3, :) == 3)), 'k*')
hold off
figure(1)
hold on
plot(th0(1, (isv(3, :) == 1)), ee(1, (isv(3, :) == 1)), 'g*')
plot(th0(1, (isv(3, :) == 2)), ee(1, (isv(3, :) == 2)), 'y*')
plot(th0(1, (isv(3, :) == 3)), ee(1, (isv(3, :) == 3)), 'k*')
hold off