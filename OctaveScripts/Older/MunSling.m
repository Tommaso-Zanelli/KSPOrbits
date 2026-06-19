%%  ----  Header  -----------------------------------------------------

fclose all;
close all
format long g
more off
clear
clc

%%  ----  Orbits and other parameters  --------------------------------

% Starting (current) time
t0 = ydhmsTimeToSeconds([0 172 0 0 0]);

% Orbital parameters
mu_Sun    = 1.1723328e18;     % [m^3s^-2]
mu_Kerbin = 3.5316000e12;     % [m^3s^-2]
mu_Mun    = 6.5138398e10;     % [m^3s^-2]

% Radi... (Radii? Radiuses?)
R_Sun    = 261600000;     % [m]
R_Kerbin = 600000;        % [m]
R_Mun    = 200000;        % [m]

% Atmosphere altitudes
ha_Sun    = 600000;     % [m]
ha_Kerbin = 70000;      % [m]
ha_Mun    = 0;          % [m]

% Initial mean anomalies
ime_Kerbin = 3.14;
ime_Mun    = 1.7;

% Kerbin orbit
[o_Kerbin] = rawOrbitToStructOrbit([[13599840256, 13599840256, 0, 0, 0, 0]; zeros(2, 6)]);
[o_Mun]    = rawOrbitToStructOrbit([[12000000, 12000000, 0, 0, 0, 0]; zeros(2, 6)]);

% Orbital periods
T_Kerbin = 2 * argStructToTime(pi, o_Kerbin, mu_Sun);        % [s]
T_Mun    = 2 * argStructToTime(pi, o_Mun,    mu_Kerbin);     % [s]

% Periapsis times
o_Kerbin.peT = -T_Kerbin * (ime_Kerbin / (2 * pi));     % [s]
o_Mun.peT    = -T_Mun    * (ime_Mun    / (2 * pi));     % [s]

% Spheres of influence
rSOI_Kerbin = ((mu_Kerbin / mu_Sun)    ^ 0.4)    * o_Kerbin.p ...
                                        / (1 - (o_Kerbin.e ^ 2));     % [m]
rSOI_Mun    = ((mu_Mun    / mu_Kerbin) ^ 0.4)    * o_Mun.p    ...
                                        / (1 - (o_Mun.e    ^ 2));     % [m]

%%  ----    ----

N = 5001;

pe  = 75000 + R_Kerbin;
ap1 = o_Mun.p - rSOI_Mun * 0.9;
ap2 = o_Mun.p + rSOI_Mun * 0.9;

ap_v = linspace(ap1, ap2, N);

%errArg = zeros(N, 2);

[xM, vM, aM] = structTimeToPosVelArg(t0, o_Mun, mu_Kerbin);
aM           = aM - 2 * pi * floor(aM / (2 * pi));

errPt = zeros(N, 3);

E0 = (mu_Kerbin / R_Kerbin);

minDist = zeros(N, 1);
VV_diff = zeros(N, 1);
VN_diff = zeros(N, 1);
EK_diff = zeros(N, 1);
ET_diff = zeros(N, 1);

for idx_1 = 1:N

    fprintf('%g of %g\n', idx_1, N)

    [o.p, o.e] = apoPeriToPe(ap_v(1, idx_1), pe);
    o.peT  = 0;
    o.lan  = 0;
    o.incl = 0;
    o.aop  = 0;
    T_2 = argStructToTime(pi, o, mu_Kerbin);
    o.peT = t0 - T_2;
    o.aop = aM - pi;

    [xo, vo, ao] = structTimeToPosVelArg(t0, o, mu_Kerbin);

    %errArg(idx_1, 1) = atan2(xo(2, 1), xo(1, 1)) - atan2(xM(2, 1), xM(1, 1));
    %errArg(idx_1, 2) = norm(xo - xM, 2) - abs(ap_v(1, idx_1) - o_Mun.p);

    [t_enc, d, th_c1, th_c2] = entryInSphere([o.peT, t0], rSOI_Mun, o, o_Mun, mu_Kerbin);

    errPt(idx_1, 1) = 1 - norm(structTimeToPosVelArg(t_enc, o, mu_Kerbin) - ...
        structTimeToPosVelArg(t_enc, o_Mun, mu_Kerbin), 2) / rSOI_Mun;

    [xMK, vMK, ~] = structTimeToPosVelArg(t_enc, o_Mun, mu_Kerbin);
    [xoK, voK, ~] = structTimeToPosVelArg(t_enc, o, mu_Kerbin);

    xoM = xoK - xMK;
    voM = voK - vMK;

    [th, oM, t2] = posVelToStruct(xoM, voM, mu_Mun, t_enc);

    th2 = acos(((oM.p / rSOI_Mun) - 1) / oM.e);

    %[th3] = timeStructToArg(t_enc, oM, mu_Mun);
    %errPt(idx_1, 2) = abs(th) - abs(th2);
    %errPt(idx_1, 3) = abs(th3) - abs(th2);
    %fprintf('%g - %g - %g - %g - %g\n', sign(th), oM.incl, errPt(idx_1, 2), ...
    %    (norm(structArgToPosVelTime(th, oM, mu_Mun), 2) - rSOI_Mun) / rSOI_Mun, ...
    %    norm(structArgToPosVelTime(th, oM, mu_Mun), 2) / norm(xoM, 2))
        %norm(structArgToPosVelTime(th, oM, mu_Mun) - xoM, 2) / rSOI_Mun)
        %(norm(structTimeToPosVelArg(t_enc, oM, mu_Mun), 2) - rSOI_Mun) / rSOI_Mun, ..
        %norm(structTimeToPosVelArg(t_enc, oM, mu_Mun) - xoM, 2) / rSOI_Mun)

    errPt(idx_1, 2) = (abs(norm(structArgToPosVelTime(th, oM, mu_Mun), 2) - rSOI_Mun) + ...
        abs(norm(structArgToPosVelTime(th2, oM, mu_Mun), 2) - rSOI_Mun)) / ...
        (2 * rSOI_Mun);

    %fprintf('%g\n', th2 - th)


    [xoM2, voM2, te] = structArgToPosVelTime(th2, oM, mu_Mun);
    [xMK2, vMK2, ~] = structTimeToPosVelArg(te, o_Mun, mu_Kerbin);

    xo2 = xoM2 + xMK2;
    vo2 = voM2 + vMK2;

    [the, o2, te_c] = posVelToStruct(xo2, vo2, mu_Kerbin, te);

    EK1 = 0.5 * (norm(voK, 2) ^ 2);
    EK2 = 0.5 * (norm(vo2, 2) ^ 2);

    E1 = EK1 - (mu_Kerbin / norm(xoK, 2)) + E0;
    E2 = EK2 - (mu_Kerbin / norm(xo2, 2)) + E0;
    %fprintf('%g\n', E2 - E1)

    minDist(idx_1, 1) = oM.p / (1 + oM.e);
    VV_diff(idx_1, 1) = norm((vo2 - voK), 2);
    VN_diff(idx_1, 1) = norm(vo2, 2) - norm(voK, 2);
    EK_diff(idx_1, 1) = EK2 - EK1;
    ET_diff(idx_1, 1) = E2 - E1;

end

figure(4)
plot(minDist, ET_diff, 'b')

figure(5)
plot(minDist, EK_diff, 'r')

figure(6)
plot(minDist, VV_diff, 'k')

figure(7)
plot(minDist, VN_diff, 'y')

%figure(1)
%plot(ap_v - o_Mun.p, errPt(:, 1), 'b')

%figure(2)
%plot(ap_v - o_Mun.p, errPt(:, 2), 'r')

%figure(3)
%plot(ap_v - o_Mun.p, errPt(:, 3), 'g')

%{
% Starting orbit details
h0K = 75000;
Tf = 0.5;

% Starting orbit structure
o0K.p    = (h0K + R_Kerbin);
o0K.e    = 0;
o0K.lan  = 0;
o0K.incl = 0;
o0K.aop  = 0;
o0K.peT  = 0;
%o0K.peT  = 2 * (Tf - 0.5) *  argStructToTime(pi, o0K, mu_Kerbin);

m1.t  = t0;

verbose = false;

N = 101;

maxV = 25000; %50000;
N1 = 1;%N;
N2 = N;
N3 = N;%N;
nm13 = 1 / (N1 * N2 * N3);
                        % Time              % Argument          % Max. Velocity
%[x1, x2, x3] = meshgrid(linspace(0, 1, N2), linspace(0, 1, N1), linspace(0, 1, N3));
[x1, x2, x3] = meshgrid(linspace(0, 1, N2), linspace(0.5, 0.5, N1), linspace(0, 1, N3));
pcc0 = -1;
if ~verbose
    fprintf('Computing orbits... (%7.4f%c)', pcc0, char(37))
end

ncp = 0;

for idx_1 = 1:N1
    for idx_2 = 1:N2
        for idx_3 = 1:N3

            pcc1 = floor((200000 * (idx_3 + N3 * (idx_2 - 1) + (N3 * N2) * (idx_1 - 1)) * nm13) + 0.5) * 0.0005;
            if pcc1 ~= pcc0

                if verbose
                    fprintf('Computing orbits... (%7.4f%c)\n\n', pcc1, char(37))
                else
                    fprintf('\b\b\b\b\b\b\b\b\b%7.4f%c)', pcc1, char(37))
                end
                pcc0 = pcc1;
            end

            Tf = x1(idx_1, idx_2, idx_3);
            th = 2 * pi * (x2(idx_1, idx_2, idx_3) - 0.5);
            vv = exp(log(maxV + 1) * x3(idx_1, idx_2, idx_3)) - 1;
            o0K.peT  = 0;
            o0K.peT  = 2 * (Tf - 0.5) *  argStructToTime(pi, o0K, mu_Kerbin);
            m1.v  = [vv * cos(th); 0; vv * sin(th)];
            [o1K] = applyManoeuvre(o0K, m1, mu_Kerbin);

            % Position, velocity and argument after manoeuvre
            [xx_m1, vv_m1, th_m1] = structTimeToPosVelArg(m1.t, o1K, mu_Kerbin);

            % Apoapsis and periapsis
            ap1 = o1K.p / (1 - o1K.e);
            pe1 = o1K.p / (1 + o1K.e);

            % Checks whether the orbit is periodic
            isPeriodic    = (o1K.e < 1) && (ap1 < rSOI_Kerbin);

            % Checks whether the object will plummet into Kerbin immediately
            doesItPlummet = (pe1 <= (R_Kerbin + ha_Kerbin)) && (th_m1 < 0);

            % If it does
            if doesItPlummet

                % Argument at the point of atmosphere entry
                th2 = -acos((1 / o1K.e) * ((o1K.p / (R_Kerbin + ha_Kerbin)) - 1));

                % Time of atmosphere entry
                t2  = argStructToTime(th2, o1K, mu_Kerbin);

                isPeriodic = false;

                if verbose
                    fprintf('\tPlummeted into Kerbin!\n')
                end

            % If it doesn't
            else

                % If it is not periodic (exits Kerbin's SOI)
                if ~isPeriodic

                    % Argument of SOI exit
                    th2 = acos((1 / o1K.e) * ((o1K.p / rSOI_Kerbin) - 1));

                    % Time of SOI exit
                    t2  = argStructToTime(th2, o1K, mu_Kerbin);

                % If it is periodic
                else

                    % If it will eventually plummet into Kerbin
                    if (pe1 <= (R_Kerbin + ha_Kerbin))

                        % Argument at the point of atmosphere entry
                        th2 = acos((1 / o1K.e) * ((o1K.p / (R_Kerbin + ha_Kerbin)) - 1));

                        % Time of atmosphere entry
                        t2  = argStructToTime(th2, o1K, mu_Kerbin);

                        isPeriodic = false;

                    % If it is just actually periodic
                    else

                        % Period
                        T = 2 * argStructToTime(pi, o1K, mu_Kerbin);

                        % Final time
                        t2 = t0 + T;

                    end

                end

            end

            % Defines search time interval
            ti = [t0, t2];

            % Searches closest period
            if isPeriodic
                if verbose
                    fprintf('\tPeriodic orbit! ')
                end

                [~, d0, ~, ~] = closestPass(ti, o1K, o_Mun, mu_Kerbin);
                n = 0;

                while true

                    [~, d1, ~, ~] = closestPass(ti + ((n + 1) * T), o1K, o_Mun, mu_Kerbin);

                    if d1 > d0
                        break;
                    else
                        n = n + 1;
                        d0 = d1;
                        if d1 <= rSOI_Mun
                            break;
                        end
                    end

                end

                ti = ti + T * n;
                if verbose
                    fprintf('Closest passage will be after %g rounds.\n', n)
                end
            end


            % Computes closest passage
            [t_cp, d_cp, th1, th2] = closestPass(ti, o1K, o_Mun, mu_Kerbin);

            % Checks whether
            if d_cp <= rSOI_Mun
                ncp = ncp + 1;
                man(ncp, :) = [x1(idx_1, idx_2, idx_3),  x2(idx_1, idx_2, idx_3),  x3(idx_1, idx_2, idx_3)];
            end

            if verbose
                fprintf('\tThe closest it comes to the Mun is: %g / %g\n\n', d_cp, rSOI_Mun)
            end

        end
    end
end

fprintf('\n')
%}
