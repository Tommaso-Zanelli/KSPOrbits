function [acb] = KerbolSystemF()
% fclose all;
% close all
% format long g
% more off
% %clear
% clc

% Universal gravitational constant
G = 6.674e-11;

% Names
Names = {'Kerbol'; 'Moho'; 'Eve'; 'Kerbin'; 'Duna'; 'Dres'; 'Jool'; ...
         'Eeloo'; 'Gilly'; 'Mun'; 'Minmus'; 'Ike'; 'Laythe'; 'Vall'; ...
         'Tylo'; 'Bop'; 'Pol'};

% Mass, Radius, rotation period, semiMajorAxis, eccentricity, incl, lan,
% aop, meanAnomalyAt...0?
Kerbol = [1.756567e+28, 2.616e+08, 432000, 0, 0, 0, 0, 0, 0];

%Kerbol
Moho = [2.5263617e21, 250000, 1210000, 5263138304, 0.2, 7.0, 70.0, 15.0, 3.14];
Eve = [1.2244127e23, 700000, 80500, 9832684544, 0.01, 2.1, 15.0, 0, 3.14];
Kerbin = [5.2915793e22, 600000, 21600, 13599840256, 0.0, 0, 0, 0, 3.14];
Duna = [4.5154812e21, 320000, 65517.859, 20726155264, 0.051, 0.06, 135.5, 0, 3.14];
Dres = [3.2191322e20, 138000, 34800, 40839348203, 0.145, 5.0, 280.0, 90.0, 3.14];
Jool = [4.2332635e24, 6000000, 36000, 68773560320, 0.05, 1.304, 52.0, 0, 0.1];
Eeloo = [1.1149358e21, 210000, 19460, 90118820000, 0.26, 6.15, 50.0, 260.0, 3.14];

%Eve
Gilly = [1.2420512e17, 13000, 28255, 31500000, 0.55, 12.0, 80.0, 10.0, 0.9];

%Kerbin
Mun = [9.7600236e20, 200000, 138984.38, 12000000, 0.0, 0, 0, 0, 1.7];
Minmus = [2.6457897e19, 60000, 40400, 47000000, 0.0, 6.0, 78.0, 38.0, 0.9];
                                                   %incl,  lan, aop
%Duna
Ike = [2.7821949e20, 130000, 65517.862, 3200000, 0.03, 0.2, 0, 0, 1.7];

%Jool
Laythe = [2.9397663e22, 500000, 52980.879, 27184000, 0, 0, 0, 0, 3.14];
Vall = [3.1088028e21, 300000, 105962.09, 43152000, 0, 0, 0, 0, 0.9];
Tylo = [4.2332635e22, 600000, 211926.36, 68500000, 0, 0.025, 0, 0, 3.14];
Bop = [3.7261536e19, 65000, 544507.4, 128500000, 0.235, 15.0, 10.0, 25.0, 0.9];
Pol = [1.0813636e19, 44000, 901902.62, 179890000, 0.17085, 4.25, 2.0, 15.0, 0.9];

% All
acb = [0, Kerbol; ...
       1, Moho;   ...
       1, Eve;    ...
       1, Kerbin; ...
       1, Duna;   ...
       1, Dres;   ...
       1, Jool;   ...
       1, Eeloo;  ...
       3, Gilly;  ...
       4, Mun;    ...
       4, Minmus; ...
       5, Ike;    ...
       7, Laythe; ...
       7, Vall;   ...
       7, Tylo;   ...
       7, Bop;    ...
       7, Pol];

% Sign of the longitude of ascnding node
s_lan  = -1;

% Sign of the inclination
s_incl = 1;

% Sign of the argument of the periapsis
s_aop  = -1;

%
acb(:, 8) = acb(:, 8) * s_lan;
acb(:, 7) = acb(:, 7) * s_incl;
acb(:, 9) = acb(:, 9) * s_aop;

% Orbit angles converted to radians
acb(:, 7:9) = (acb(:, 7:9) * 2 * pi / 180);

% Adds a field containing the gravitational constant of each body
acb = [acb, G * acb(:, 2)];

% Adds a field for the semi-latus...pfft... rectum
acb = [acb, (acb(:, 5).* (1 - (acb(:, 6).^ 2)))];

% Adds a field for the initial anomaly
acb = [acb, zeros(size(acb, 1), 1)];

% Computes the initial anomaly
N = size(acb, 1);
for idx_1 = 1:N
    ec = acb(idx_1, 6);
    Mc = acb(idx_1, 10);
    q = @(E) (E - (ec * sin(E))) - Mc;
    Ec = fzero(q, Mc);
    acb(idx_1, 13) = 2 * atan(sqrt((1 + ec) / (1 - ec)) * tan(0.5 * Ec));
end

% Adds a field for the "periapsis time"
acb = [acb, zeros(size(acb, 1), 1)];

% Computes "periapsis time"
for idx_1 = 2:N
    acb(idx_1, 14) = -argToTime(pi, acb(idx_1, 12), acb(idx_1, 6), acb(acb(idx_1, 1), 11), 0) * acb(idx_1, 10) / pi;
end

% %positions of all bodies at instant t:
% t = 2916000;
% cc = zeros(N, 3);
% cc = [255, 255, 0; ...
%     229, 123, 53; ...
%     100, 60, 140; ...
%     118, 163, 89; ...
%     177, 9, 8; ...
%     56, 56, 56; ...
%     0, 255, 0; ...
%     163, 163, 163; ...
%     232, 48, 217; ...
%     107, 107, 107; ...
%     185, 203, 162; ...
%     225, 225, 225; ...
%     167, 193, 214; ...
%     101, 152, 193; ...
%     163, 163, 163; ...
%     24, 170, 220; ...
%     223, 9, 27] / 255;
% cc2 = 1 - 0.25 * (1 - cc);
% ths = linspace(-pi, pi, 5001);
% p = zeros(3, N);
% plot3(0, 0, 0, 'o', 'color', cc(1, :), 'LineWidth', 3)
% hold on
% for idx_1 = 2:N
%     th = timeToArg(t, acb(idx_1, 12), acb(idx_1, 6), acb(acb(idx_1, 1), 11), acb(idx_1, 14));
%     [p(:, idx_1), ~] = paramsToPosVel(th, acb(idx_1, 12), acb(idx_1, 6), acb(idx_1, 8), acb(idx_1, 7), acb(idx_1, 9), acb(acb(idx_1, 1), 11));
%     [xx, ~] = paramsToPosVel(ths, acb(idx_1, 12), acb(idx_1, 6), acb(idx_1, 8), acb(idx_1, 7), acb(idx_1, 9), acb(acb(idx_1, 1), 11));
%     p(:, idx_1) = p(:, idx_1) + p(:, acb(idx_1, 1));
%     xx = xx + repmat(p(:, acb(idx_1, 1)), 1, 5001);
%     if acb(idx_1, 1) == 1
%         lw = 2;
%     else
%         lw = 1;
%     end
%     plot3(p(1, idx_1), p(2, idx_1), p(3, idx_1),  'o', 'color', cc(idx_1, :), 'LineWidth', lw)
%     plot3(xx(1, :), xx(2, :), xx(3, :), 'color', cc2(idx_1, :))
% end
% axis equal
% hold off

% Q = [0, -1.999999999999993
% -2.3333333333333286, 0
% -1, -3.9999999999999893
% -0.3333333333333286, 7.999999999999984
% 7.999999999999979, -8.666666666666647
% 20.999999999999957, 6.666666666666654];
