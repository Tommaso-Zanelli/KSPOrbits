%% --- Header -------------------------------------------------------------

fclose all;
close all
format long g
more off
clear
clc

%% --- Orbits and other parameters ----------------------------------------

% Starting (current) time
t0 = ydhmsTimeToSeconds([0 164 3 22 30]);

% Orbital parameters
mu_Sun    = 1.1723328e18;     % [m^3s^-2]
mu_Kerbin = 3.5316000e12;     % [m^3s^-2]
mu_Mun    = 6.5138398e10;     % [m^3s^-2]
mu_Duna   = 3.0136321e11;     % [m^3s^-2]
mu_Ike    = 1.8568369e10;     % [m^3s^-2]

% Radi... (Radii? Radiuses?)
R_Sun    = 261600000;     % [m]
R_Kerbin = 600000;        % [m]
R_Mun    = 200000;        % [m]
R_Duna   = 320000;        % [m]
R_Ike    = 130000;        % [m]

% Atmosphere altitudes
ha_Sun    = 600000;     % [m]
ha_Kerbin = 70000;      % [m]
ha_Mun    = 0;          % [m]
ha_Duna   = 50000;      % [m]
ha_Ike    = 0;          % [m]

% Initial mean anomalies
ime_Kerbin = 3.14;
ime_Mun    = 1.7;
ime_Duna   = 3.14;
ime_Ike    = 1.7;

% Kerbin orbit
[o_Kerbin] = rawOrbitToStructOrbit([[13599840256, 13599840256, 0, 0, 0, 0]; zeros(2, 6)]);
[o_Mun]    = rawOrbitToStructOrbit([[12000000, 12000000, 0, 0, 0, 0]; zeros(2, 6)]);
[o_Duna]   = rawOrbitToStructOrbit([[19669121365, 21783189163, 0.051, 135.5, 0.06, 0]; zeros(2, 6)]);
[o_Ike]    = rawOrbitToStructOrbit([[3104000, 3296000, 0.03, 0, 0.2, 0]; zeros(2, 6)]);

% Orbital periods
T_Kerbin = 2 * argStructToTime(pi, o_Kerbin, mu_Sun);        % [s]
T_Mun    = 2 * argStructToTime(pi, o_Mun,    mu_Kerbin);     % [s]
T_Duna   = 2 * argStructToTime(pi, o_Duna,   mu_Sun);        % [s]
T_Ike    = 2 * argStructToTime(pi, o_Ike,    mu_Duna);       % [s]

% Periapsis times
o_Kerbin.peT = -T_Kerbin * (ime_Kerbin / (2 * pi));     % [s]
o_Mun.peT    = -T_Mun    * (ime_Mun    / (2 * pi));     % [s]
o_Duna.peT   = -T_Duna   * (ime_Duna   / (2 * pi));     % [s]
o_Ike.peT    = -T_Ike    * (ime_Ike    / (2 * pi));     % [s]

% Starting orbit
[o_Start] = rawOrbitToStructOrbit([[67900, 680100, 0, 0, 0, 0]; zeros(2, 6)]);
o_Start.peT = t0;


%% --- Graphic output -----------------------------------------------------

% Static orbits
thFull = linspace(-pi, pi, 5001);
[x_Kerbin_Static, ~, ~] = structArgToPosVelTime(thFull, o_Kerbin, mu_Sun);
[x_Mun_Static, ~, ~]    = structArgToPosVelTime(thFull, o_Mun,    mu_Kerbin);
[x_Duna_Static, ~, ~]   = structArgToPosVelTime(thFull, o_Duna,   mu_Sun);
[x_Ike_Static, ~, ~]    = structArgToPosVelTime(thFull, o_Ike,    mu_Duna);

[x_Kerbin_Now, ~, ~] = structTimeToPosVelArg(t0, o_Kerbin, mu_Sun);
[x_Mun_Now, ~, ~]    = structTimeToPosVelArg(t0, o_Mun, mu_Kerbin);
[x_Duna_Now, ~, ~]   = structTimeToPosVelArg(t0, o_Duna,   mu_Sun);
[x_Ike_Now, ~, ~]    = structTimeToPosVelArg(t0, o_Ike,   mu_Duna);




figure(1)
hold on
axis equal
plot3(x_Kerbin_Static(1, :), x_Kerbin_Static(2, :), x_Kerbin_Static(3, :), 'b')
plot3(x_Duna_Static(1, :), x_Duna_Static(2, :), x_Duna_Static(3, :), 'r')
plot3(x_Kerbin_Now(1, :), x_Kerbin_Now(2, :), x_Kerbin_Now(3, :), 'o', 'color', [0, 0, 0.5], 'MarkerSize', 6, 'LineWidth', 3)
plot3(x_Duna_Now(1, :), x_Duna_Now(2, :), x_Duna_Now(3, :), 'o', 'color', [0.5, 0, 0], 'MarkerSize', 6, 'LineWidth', 3)
plot3(0, 0, 0, 'o', 'color', [0.5, 0.5, 0], 'MarkerSize', 6, 'LineWidth', 3)

plot3(x_Kerbin_Now(1, :) + x_Mun_Static(1, :), x_Kerbin_Now(2, :) + x_Mun_Static(2, :), x_Kerbin_Now(3, :) + x_Mun_Static(3, :), 'g')
plot3(x_Kerbin_Now(1, :) + x_Mun_Now(1, :), x_Kerbin_Now(2, :) + x_Mun_Now(2, :), x_Kerbin_Now(3, :) + x_Mun_Now(3, :), 'o', 'color', [0, 0.5, 0], 'MarkerSize', 6, 'LineWidth', 3)
plot3(x_Duna_Now(1, :) + x_Ike_Static(1, :), x_Duna_Now(2, :) + x_Ike_Static(2, :), x_Duna_Now(3, :) + x_Ike_Static(3, :), 'color', [0.5, 0.5, 0.5])
plot3(x_Duna_Now(1, :) + x_Ike_Now(1, :), x_Duna_Now(2, :) + x_Ike_Now(2, :), x_Duna_Now(3, :) + x_Ike_Now(3, :), 'o', 'color', [0.25, 0.25, 0.25], 'MarkerSize', 6, 'LineWidth', 3)
hold off