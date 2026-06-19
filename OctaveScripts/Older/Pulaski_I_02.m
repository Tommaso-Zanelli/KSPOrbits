%% - Header ---------------------------------------------------------------

fclose all;
close all
format long g
more off
clear
clc

grph_out = true;

%% - Veichle capacity -----------------------------------------------------

% Delta-vs and corresponding durations for the three stages
dvt1 = [0 10 20 30 40 50 60 70 76.59; 0 25 51 77 103 129 155 180 197]';
dvt2 = [0 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200 210 220 230 240 250 260 270 280 290 300 310 320 330 340 350 360 370 380 390 400 410 420 430 440 450 460 470 480 490 500 510 520 530 540 550 560 570 580 590 600 610 620 630 630.54; 0 22 44 66 88 110 132 154 176 198 220 241 263 285 307 328 350 372 394 415 437 458 480 501 523 544 566 587 609 630 651 672 694 715 736 757 779 800 821 842 863 884 905 926 947 968 989 1009 1030 1051 1072 1093 1113 1134 1155 1175 1196 1217 1237 1258 1278 1299 1319 1340 1341]';
dvt3 = [0 50 100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300; 0 4 9 14 18 23 27 31 36 40 44 48 52 56 60 64 67 71 75 78 82 85 88 92 95 98 101]';

g0    = 9.81;               % [m / s^2] Reference acceleration
Isp   = 800;                % [s]       Specific impulse for the first two stages (nuclear engine)
Isp3  = 330;                % [s]       Specific impulse for the third stage
dmdt  = 1.53 * 5;           % [kg / s]  Mass flow for the first two stages (nuclear engine)
dmdt3 = 18.54 * (100 / 9);  % [kg / s]  Mass flow for the third stage

% Mass ratios for each possible delta-v considered: Dv = g0*Isp*ln(mr) -> mr = exp(Dv) / (go*Isp)
mr1 = exp(dvt1(2:end, 1) / (Isp * g0));
mr2 = exp(dvt2(2:end, 1) / (Isp * g0));
mr3 = exp(dvt3(2:end, 1) / (Isp3 * g0));

% Mass differences  for each possible delta-v considered: md = dt * dmdt
md1 = dvt1(2:end, 2) * dmdt;
md2 = dvt2(2:end, 2) * dmdt;
md3 = dvt3(2:end, 2) * dmdt3;

% Initial intervals for the masses of each initial stage
%{
M1_0 = [min(mr1.* md1./ (mr1 - 1)), max(mr1.* md1./ (mr1 - 1))];
M2_0 = [min(mr2.* md2./ (mr2 - 1)), max(mr2.* md2./ (mr2 - 1))];
M3_0 = [min(mr3.* md3./ (mr3 - 1)), max(mr3.* md3./ (mr3 - 1))];

% "cost functions" for the masses of each initial stage
cferrM1 = @(M) sum(((1 - 1./ mr1) * M / dmdt - dvt1(2:end, 2)).^ 2, 1);
cferrM2 = @(M) sum(((1 - 1./ mr2) * M / dmdt - dvt2(2:end, 2)).^ 2, 1);
cferrM3 = @(M) sum(((1 - 1./ mr3) * M / dmdt3 - dvt3(2:end, 2)).^ 2, 1);

% Initial masses of each stage [kg]
M1i = optiFun(M1_0, cferrM1, 1e-3, 16, 1e-14, 2500);
M2i = optiFun(M2_0, cferrM2, 1e-3, 16, 1e-14, 2500);
M3i = optiFun(M3_0, cferrM3, 1e-3, 16, 1e-14, 2500);
%}

% Previously computed values:
M1i = 155156.762505796222593635321;
M2i = 132848.541751198866473510861;
M3i =  63239.500950567175447940826;

% Corrected times [s]
tm1 = [0; (M1i * (1 - 1./ mr1) / dmdt)];
tm2 = [0; (M2i * (1 - 1./ mr2) / dmdt)];
tm3 = [0; (M3i * (1 - 1./ mr3) / dmdt3)];

% Time errors
%{
fprintf('Time error for the first stage  :  %g [s]\n', sum(abs(floor(tm1 + 0.5) - dvt1(:, 2)), 1))
fprintf('Time error for the second stage :  %g [s]\n', sum(abs(floor(tm2 + 0.5) - dvt2(:, 2)), 1))
fprintf('Time error for the third stage  :  %g [s]\n', sum(abs(floor(tm3 + 0.5) - dvt3(:, 2)), 1))
%}

% Total times [s]
t_stage1 = tm1(end, 1);
t_stage2 = tm2(end, 1);
t_stage3 = tm3(end, 1);

% Final masses of each stage [kg]
M1f = M1i - t_stage1 * dmdt;
M2f = M2i - t_stage2 * dmdt;
M3f = M3i - t_stage3 * dmdt3;

% Other masses [kg]
M1  = M1i;                  % [kg]      Initial mass
dM2 = M1f - M2i;            % [kg]      Dead weight discarded with first stage
dM3 = M2f - M3i;            % [kg]      Dead weight discarded with second stage


%% - Timings --------------------------------------------------------------

mtm = 91;       % [s] Maximum allowed manoeuvre duration
Dt  = 1.5;      % [s] Time interval between manoeuvres
Dt2 = 30;       % [s] Time interval when stage separation takes place

% Number of manoeuvres for each stage
nms1 = ceil(t_stage1 / mtm);
nms2 = ceil(t_stage2 / mtm);
nms3 = ceil(t_stage3 / mtm);

% Duration of individual manoeuvres for the first, second and third stages respectively [s]
tt1 = (t_stage1 / nms1) * ones(1, nms1);
tt2 = (t_stage2 / nms2) * ones(1, nms2);
tt3 = (t_stage3 / nms3) * ones(1, nms3);

%Number of manoeuvres
Nm1 = size(tt1, 2);         % Number of manoeuvres for the first stage
Nm2 = size(tt2, 2);         % Number of manoeuvres for the second stage
Nm3 = size(tt3, 2);         % Number of manoeuvres for the third stage
Nm  = (Nm1 + Nm2 + Nm3);    % Total number of manoeuvres

% Mass differences [kg]
dm1 = tt1 * dmdt;
dm2 = tt2 * dmdt;
dm3 = tt3 * dmdt3;

% Mass arrays [kg]
mm1 = M1i:-(t_stage1 * dmdt / nms1):M1f;
mm2 = M2i:-(t_stage2 * dmdt / nms2):M2f;
mm3 = M3i:-(t_stage3 * dmdt3 / nms3):M3f;

% Delta-vs for each manoeuvre [m / s]
dv1 = g0 * Isp * log(mm1(1, 1:(end - 1))./ mm1(1, 2:end));
dv2 = g0 * Isp * log(mm2(1, 1:(end - 1))./ mm2(1, 2:end));
dv3 = g0 * Isp3 * log(mm3(1, 1:(end - 1))./ mm3(1, 2:end));

% Initializes the manoeuvre time array and the delta-vs array
ttm = zeros(1, Nm);
dvm = zeros(1, Nm);

% Fills the manoeuvre time array
dvm(1, 1) = dv1(1, 1);
for idx_1 = 2:Nm1
    ttm(1, idx_1) = ttm(1, (idx_1 - 1)) + 0.5 * (tt1(1, idx_1) + tt1(1, (idx_1 - 1))) + Dt;
    dvm(1, idx_1) = dv1(1, idx_1);
end
ttm(1, (Nm1 + 1)) = ttm(1, Nm1) + 0.5 * (tt1(1, Nm1) + tt2(1, 1)) + Dt2;
dvm(1, (Nm1 + 1)) = dv2(1, 1);
for idx_1 = 2:Nm2
    ttm(1, (Nm1 + idx_1)) = ttm(1, (Nm1 + idx_1 - 1)) + 0.5 * (tt2(1, idx_1) + tt2(1, (idx_1 - 1))) + Dt;
    dvm(1, (Nm1 + idx_1)) = dv2(1, idx_1);
end
ttm(1, (Nm1 + Nm2 + 1)) = ttm(1, (Nm1 + Nm2)) + 0.5 * (tt2(1, Nm2) + tt3(1, 1)) + Dt2;
dvm(1, (Nm1 + Nm2 + 1)) = dv3(1, 1);
for idx_1 = 2:Nm3
    ttm(1, (Nm1 + Nm2 + idx_1)) = ttm(1, (Nm1 + Nm2 + idx_1 - 1)) + 0.5 * (tt3(1, idx_1) + tt3(1, (idx_1 - 1))) + Dt;
    dvm(1, (Nm1 + Nm2 + idx_1)) = dv3(1, idx_1);
end

% Total manoeuvre time [s]
mT = ttm(1, Nm);


%% - Orbit details --------------------------------------------------------

% Kerbin details
mu_K  = 3.5316e12;
R_K = 600000;

% Initial orbit details
ap0   = -300.710e6 + R_K;           % Apoapsis  [m]
pe0   = 74.9998e3  + R_K;           % Periapsis [m]
[p0, e0] = apoPeriToPe(ap0, pe0);   % Semi-latus rectum [m] and eccentricity
lan0  = 147.0 * (pi / 180);         % Longitude of ascending node
incl0 = 7.0   * (pi / 180);         % Orbital inclination
aop0  = 34.4  * (pi / 180);         % Argument of periapsis
pe0T  = 2856546.578;                % Periapsis absolute time [s]

% Outer limits for the initial orbit
[thl1] = timeToArg((pe0T - 1.5 * mT), p0, e0, mu_K, pe0T);
[thl2] = timeToArg((pe0T + 1.5 * mT), p0, e0, mu_K, pe0T);

% Defines a cost function corresponding to the third coordinate
cferrPlanar = @(th) (dot([0; 0; 1], paramsToPosVel(th, p0, e0, lan0, incl0, aop0, mu_K)) ^ 2);

% Finds the argiment of the equatorial point of the orbit
%{
th0fM   = optiFun([thl1, thl2], cferrPlanar, 1e-8, 16, 1e-15, 2500);
th0fM   = th0fM + eps(th0fM);
%}
% Pre-computed value
th0fM = -0.600393262686049357086;

% Time of the final manoeuvre (will change!)
Tmf   = argToTime(th0fM, p0, e0, mu_K, pe0T);

% Position and velocity in the original orbit at the final manoeuvre node
[x0fM, v0fM] = paramsToPosVel(th0fM, p0, e0, lan0, incl0, aop0, mu_K);

% Final orbit details
th_f  = atan2(x0fM(2, 1), x0fM(1, 1));              % Argument of the intersection
pf    = norm(x0fM, 2);                              % Semi-latus rectum [m]
ef    = 0.0;                                        % Eccentricity
lanf  = 0.0 * (pi / 180);                           % Longitude of ascending node
inclf = 0.0 * (pi / 180);                           % Orbital inclination
aopf  = 0.0 * (pi / 180);                           % Argument of periapsis
pefT  = Tmf - argToTime(th_f, pf, ef, mu_K, 0);     % Periapsis absolute time [s]

% Position and velocity in the final target orbit at the final manoeuvre node
[xffM, vffM] = paramsToPosVel(th_f, pf, ef, lanf, inclf, aopf, mu_K);

% Final dv (will change!)
dvf = norm((v0fM - vffM), 2);

% First manoeuvre node (fixed)
Tmi = Tmf - mT;                                 % Time of the manoeuvre [s]
thmi = timeToArg(Tmi, p0, e0, mu_K, pe0T);      % Argument of the manoeuvre
% Position and velocity at the node
[xmi, vmi] = paramsToPosVel(thmi, p0, e0, lan0, incl0, aop0, mu_K);

% Offsets the manoeuvre times
ttm_offs = ttm + Tmi;

%% - Manoeuvre design -----------------------------------------------------

%{
% Structure of the data array containing all intermetiade orbits
% Mdata(:, :, i) = [t1	t2	tf	p
%                   th1	th2	thf	e
%                   x1	x2	xf	lan
%                   y1	y2	yf	incl
%                   z1	z2	zf	aop
%                   u1	u2	uf	t_pe
%                   v1	v2	vf	av1
%                   w1	w2	wf	av2]

Mdata = zeros(8, 4, (Nm + 1));

% Sets the initial orbit
Mdata(:, :, 1) = [zeros(8, 1), ...
                  [ttm_offs(1, 1); thmi; xmi; vmi], ...
                  [Tmf; th0fM; x0fM; v0fM], ...
                  [p0; e0; lan0; incl0; aop0; pe0T; 0; 0]];

% Sets the final orbit
Mdata(:, :, (Nm + 1)) = [[0; th_f; xffM; vffM], ...
                         zeros(8, 1), ...
                         [0; th_f; xffM; vffM], ...
                         [pf; ef; lanf; inclf; aopf; pefT; 0; 0]];

for idx_1 = 1:2%(Nm1 + Nm2 - 1)

    xmi_c = Mdata(3:5, 2, idx_1); %xmi; % First manoeuvre node position
    vmi_c = Mdata(6:8, 2, idx_1);%vmi; % Current velocity at the manoeuvre node


    % Cost function
    cf = @(v1t, v2t) ((norm((vmi_c - v1t), 2) - dv1(1, 1)) ^ 2) + (1e308) * (norm((vffM - v2t), 2) > dvf);

    np = 16;        % The initial number of test points
    cc = 1e308;     % Initial "minimum"

    % Iterates until a decent minimum is found
    while (abs(cc) > 1e-10)

        % Finds the optimal transfer orbit for the given cost function
        [th_1t, th_2t, pt, et, lant, inclt, aopt] = optiTransfer2(xmi_c, x0fM, cf, mu_K, 4.7e-6, np);

        % Position and velocity at the extremes of the orbit
        [xx, vv] =  paramsToPosVel([th_1t, th_2t], pt, et, lant, inclt, aopt, mu_K);

        % Actual cost of the optimum
        cc = cf(vv(:, 1), vv(:, 2)) + norm((xx(:, 1) - xmi_c), 2) / norm(xmi, 2) + norm((xx(:, 2) - x0fM), 2) / norm(x0fM, 2);

        % Increses the number of test points for the next "round"
        np = np * 2;

        % Aborts if too many attempts are made
        if (np > 32768)
            fprintf('Error: transfer orbit n. %g not found!\n', 1)
            return;
        end

    end

    % Outputs the current deltav
    fprintf('%g, ', dvf)


    % Updates the velocity difference at the final manoeuvre node
    dvf = norm((vv(:, 2) - vffM), 2);

    % Computes the times of the two manoeuvres in the transfer orbit
    t_tmp = argToTime([th_1t, th_2t], pt, et, mu_K, 0);

    % Derives the (absolute) periapsis time of the transfer orbit from the initial manoeuvres time
    petT = ttm_offs(1, 1) - t_tmp(1, 1);

    % Finds the new (absolute) final manoeuvre time
    NFMT = petT + t_tmp(1, 2);

    (NFMT - ttm_offs(1, end))

    % Scales all subsequent manoeuvre times accordingly
    ttm_offs(1, (idx_1 + 1):end) = ttm_offs(1, (idx_1 + 1):end) + (NFMT - ttm_offs(1, end));

   % Outputs the time between the two manoeuvres
   fprintf('%g\n', (ttm_offs(1, (idx_1 + 1)) - ttm_offs(1, idx_1)))

    % Finds the argument of the subsequent manoeuvre
    thm2 = timeToArg(ttm_offs(1, (idx_1 + 1)), pt, et, mu_K, petT);

    % Finds position and velocity at the subsequent manoeuvre node
    [xx2, vv2] =  paramsToPosVel(thm2, pt, et, lant, inclt, aopt, mu_K);

    % Assigns the current transfer orbit data to the transfer orbit array
    Mdata(:, :, (idx_1 + 1)) = [[ttm_offs(1, idx_1); th_1t; xx(:, 1); vv(:, 1)], ...
                                [ttm_offs(1, (idx_1 + 1)); thm2; xx2; vv2], ...
                                [NFMT; th_2t; xx(:, 2); vv(:, 2)], ...
                                [pt; et; lant; inclt; aopt; petT; 0; 0]];

end

fprintf('%g\n', dvf)
%}

%% - Manoeuvre design, third attempt --------------------------------------

% Defines initial orbit object
o0 = paramsToStruct(p0, e0, lan0, incl0, aop0, pe0T);

% Base cost function 1 (altitude at the end)
cfb1 = @(A, T, da1, da2, h) (sqrt(dot([1; 1; 1; 0; 0; 0], ...
                             Pulaski_I_02_f01(o0, ttm_offs, dvm, Nm1, ...
                             Nm2, mu_K, A, T, da1, da2).^ 2)) - h) ^ 2;

% Base cost function 2 (z coordinate at the end)
cfb2 = @(T, da1, da2, h)    ((dot([0; 0; 1; 0; 0; 0], ...
                             Pulaski_I_02_f01(o0, ttm_offs, dvm, Nm1, ...
                             Nm2, mu_K, Pulaski_I_02_f02(cfb1, T, da1, ...
                             da2, h), T, da1, da2)))) ^ 2;

% Base cost function 3 (delta-v still necessary)
cfb3 = @(A, T, da1, da2, vf) dot([0; 0; 0; 1; 1; 1], ...
                             (Pulaski_I_02_f01(o0, ttm_offs, dvm, Nm1, ...
                             Nm2, mu_K, A, T, da1, da2) - [0; 0; 0; vf]).^ 2);

% Optifun is fun


% Applies all manoeuvres (using the nuclear engine)
A = Pulaski_I_02_f02(cfb1, 0, 0, 0, (75e3 + R_K));
[os, ths] = applyMultipleManoeuvres(o0, ttm_offs(1, 1:(Nm1 + Nm2)), ...
            dvm(1, 1:(Nm1 + Nm2)), A * ones(1, (Nm1 + Nm2)), ...
            zeros(1, (Nm1 + Nm2)), mu_K);
ths(2, (Nm1 + Nm2 + 1)) = timeToArg(ttm_offs(1, (Nm1 + Nm2 + 1)), ...
        os(1, (Nm1 + Nm2 + 1)).p, os(1, (Nm1 + Nm2 + 1)).e, mu_K, ...
        os(1, (Nm1 + Nm2 + 1)).peT);

% Builds all the transfers
Np    = 101;
x_tps = zeros(3, (Nm1 + Nm2 + 1));
x_tos = zeros(3, Np * (Nm1 + Nm2));
[x_tps(:, 1), ~] = paramsToPosVel(ths(2, 1), os(1, 1).p, os(1, 1).e, ...
                   os(1, 1).lan, os(1, 1).incl, os(1, 1).aop, mu_K);
   cth = linspace(ths(1, idx_1), ths(2, idx_1), Np);
for idx_1 = 2:(Nm1 + Nm2 + 1)
   idx_2 = idx_1 - 1;
   idx_3 = (idx_2 - 1) * Np + 1;
   idx_4 = idx_2 * Np;

   [x_tps(:, idx_1), ~] = paramsToPosVel(ths(2, idx_1), os(1, idx_1).p, ...
                          os(1, idx_1).e, os(1, idx_1).lan, ...
                          os(1, idx_1).incl, os(1, idx_1).aop, mu_K);
   cth = linspace(ths(1, idx_1), ths(2, idx_1), Np);
   [x_tos(:, idx_3:idx_4), ~] = paramsToPosVel(cth, os(1, idx_1).p, ...
                                os(1, idx_1).e, os(1, idx_1).lan, ...
                                os(1, idx_1).incl, os(1, idx_1).aop, mu_K);
end
%[x_tps(:, (Nm1 + Nm2 + 1)), ~] = paramsToPosVel(ths(2, (Nm1 + Nm2 + 1)), ...
%               os(1, (Nm1 + Nm2 + 1)).p, os(1, (Nm1 + Nm2 + 1)).e, ...
%               os(1, (Nm1 + Nm2 + 1)).lan, os(1, (Nm1 + Nm2 + 1)).incl, ...
%               os(1, (Nm1 + Nm2 + 1)).aop, mu_K);

%% - Graphic output -------------------------------------------------------

% Kerbin surface
[Ksph_X, Ksph_Y, Ksph_Z] = genSphere(R_K, 51);

% Arrays of arguments for initial and final orbits
th0s = linspace(thl1, thl2, 501);
thfs = linspace(-pi, pi, 501);

% Coordinates for initial and final orbits
[x0, ~] = paramsToPosVel(th0s, p0, e0, lan0, incl0, aop0, mu_K);
[xf, ~] = paramsToPosVel(thfs, pf, ef, lanf, inclf, aopf, mu_K);

% Creates a set of points in the original orbit corresponding to the
% initial guesses for the positions of the manoeuvres
th_man     =  timeToArg((ttm + Tmi), p0, e0, mu_K, pe0T);
[x_man, ~] = paramsToPosVel(th_man, p0, e0, lan0, incl0, aop0, mu_K);

% Graph
if grph_out

    figure(1)
    plot3(x0(1, :), x0(2, :), x0(3, :), 'b', 'LineWidth', 2)
    hold on
    plot3(xf(1, :), xf(2, :), xf(3, :), 'r', 'LineWidth', 1)

    % Initial manoeuvre node
    plot3(xmi(1, 1), xmi(2, 1), xmi(3, 1), 'bo')
    plot3(xmi(1, 1), xmi(2, 1), xmi(3, 1), 'k*')

    % Final manoeuvre node
    plot3(x0fM(1, 1), x0fM(2, 1), x0fM(3, 1), 'bo')
    plot3(x0fM(1, 1), x0fM(2, 1), x0fM(3, 1), 'rs')
    plot3(x0fM(1, 1), x0fM(2, 1), x0fM(3, 1), 'k*')

    % First guess positions of the manoeuvres along the initial orbit
    %plot3(x_man(1, :), x_man(2, :), x_man(3, :), 'ks')

    % Plots manouvres
    plot3(x_tps(1, :), x_tps(2, :), x_tps(3, :), 'ys')
    plot3(x_tps(1, :), x_tps(2, :), x_tps(3, :), 'g*')
    plot3(x_tos(1, :), x_tos(2, :), x_tos(3, :), 'm', 'LineWidth', 1.5)

    % Surface of Kerbin
    surf(Ksph_X, Ksph_Y, Ksph_Z, zeros(size(Ksph_X)));
    axis equal
    hold off

end

