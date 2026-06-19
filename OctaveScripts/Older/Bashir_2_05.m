%% --- Header -------------------------------------------------------------

fclose all;
close all
format long g
more off
clear
clc

optimize = [1, 1, 1, 1];
optimize2 = [1, 1, 1, 1];

rToD = 57.29577951308232286464772187173366546630859375;

%% --- Reference data -----------------------------------------------------

% Reference data for the first satellite
ref1 = [-5.67158e7, 1.62134e5, 1.0275, 297.0, 112.5, 013.9, 600000; ...
-1, 0, 0, 5, 16, NaN, 3.5316e12; ...
1, 0, 139, 3, 41, 16, 60000; ...
1, 0, 2, 0, 33, NaN, 1.7658e9; ...
1, 0, 139, 3, 42, 8, 2247428.4; ...
-1.18609e6, 1.05530e6, 207.6703, 088.1, 093.1, 268.0, 0; ...
1, 0, 2, 1, 30, NaN, 0; ...
1, 0, 139, 3, 41, 52, 0];

% Reference data for the second satellite
ref2 = [-4.65898e7, 4.90744e5, 1.0486, 015.1, 004.6, 302.9, 600000; ...
-1, 0, 0, 5, 19, 6, 3.5316e12; ...
1, 0, 139, 3, 42, 54, 60000; ...
1, 0, 2, 0, 57, NaN, 1.7658e9; ...
1, 0, 139, 3, 43, 16, 2247428.4; ...
-2.16502e6, 2.03449e6, 398.8652, 268.5, 150.1, 086.4, 0; ...
1, 0, 2, 1, 19, NaN, 0; ...
1, 0, 139, 3, 43, 59, 0];

% Reference data for the third satellite
ref3 = [-5.60756e7, 1.68758e5, 1.0281, 297.7, 018.8, 013.4, 600000; ...
-1, 0, 0, 5, 19, NaN, 3.5316e12; ...
1, 0, 139, 3, 44, 47, 60000; ...
1, 0, 2, 0, 27, NaN, 1.7658e9; ...
1, 0, 139, 3, 45, 21, 2247428.4; ...
-1.07817e6, 10.47406e5, 188.1318, 088.1, 096.8, 268.2, 0; ...
1, 0, 2, 1, 25, 0, 0; ...
1, 0, 139, 3, 45, 36, 0];

% Reference data for the fourth satellite
ref4 = [-6.89642e7, -2.59562e5, 1.0100, 037.1, 004.9, 265.2, 600000; ...
-1, 0, 0, 5, 16, NaN, 3.5316e12; ...
1, 0, 139, 3, 45, 53, 60000; ...
1, 0, 2, 1, 5, NaN, 1.7658e9; ...
1, 0, 139, 3, 46, 13, 2247428.4; ...
-2.15661e6, 2.02561e6, 380.0190, 267.8, 029.6, 086.4, 0; ...
1, 0, 2, 1, 29, NaN, 0; ...
1, 0, 139, 3, 46, 27, 0];

%% --- Minmus' orbit ------------------------------------------------------

R_K     = ref1(1, 7);
mu_K    = ref1(2, 7);
R_M     = ref1(3, 7);
mu_M    = ref1(4, 7);
rSOI_M  = ref1(5, 7);
[oM, ~] = rawOrbitToStruct([47000000, 47000000, 0, 78, 6, 38], 0);
TM_M    = argStructToTime(pi, oM, mu_K);
oM.peT  = -TM_M * 0.9 / pi;

%% --- Refines each encounter ---------------------------------------------

dv1 = [0; 0; 0];
if exist('Bashir_2_05_Sat1Data_02.mat', 'file')
    load('Bashir_2_05_Sat1Data_02.mat');
else
    a_Sat1 = zeros(12, 1);
end

if optimize(1, 1)
    cf1 = @(x) refineEncounter_cost(x, ref1, oM, true, 1);
    fprintf('Optimizing first satellite...\n')
    a_Sat1 = optiGrad(cf1, a_Sat1, 1e-9, 1e-7, 30);
    save('Bashir_2_05_Sat1Data_02.mat', 'a_Sat1');
end

dv2 = [0; 0; 0];
if exist('Bashir_2_05_Sat2Data_02.mat', 'file')
    load('Bashir_2_05_Sat2Data_02.mat');
else
    a_Sat2 = zeros(12, 1);
end

if optimize(1, 2)
    cf2 = @(x) refineEncounter_cost(x, ref2, oM, true, 1);
    fprintf('Optimizing second satellite...\n')
    a_Sat2 = optiGrad(cf2, a_Sat2, 1e-9, 1e-7, 30);
    save('Bashir_2_05_Sat2Data_02.mat', 'a_Sat2');
end

dv3 = [0; 0; 0];
if exist('Bashir_2_05_Sat3Data_02.mat', 'file')
    load('Bashir_2_05_Sat3Data_02.mat');
else
    a_Sat3 = zeros(12, 1);
end

if optimize(1, 3)
    cf3 = @(x) refineEncounter_cost(x, ref3, oM, true, 1);
    fprintf('Optimizing third satellite...\n')
    a_Sat3 = optiGrad(cf3, a_Sat3, 1e-9, 1e-7, 30);
    save('Bashir_2_05_Sat3Data_02.mat', 'a_Sat3');
end

dv4 = [0; 0; 0];
if exist('Bashir_2_05_Sat4Data_02.mat', 'file')
    load('Bashir_2_05_Sat4Data_02.mat');
else
    a_Sat4 = zeros(12, 1);
end

if optimize(1, 4)
    cf4 = @(x) refineEncounter_cost(x, ref4, oM, true, 1);
    fprintf('Optimizing fourth satellite...\n')
    a_Sat4 = optiGrad(cf4, a_Sat4, 1e-9, 1e-7, 30);
    save('Bashir_2_05_Sat4Data_02.mat', 'a_Sat4');
end

%% --- Computes the actual orbits -----------------------------------------

[cc1, o_Sat1_K, o_Sat1_M, t_enc_Sat1] = refineEncounter_cost(a_Sat1, ref1, oM, true, 0);
[cc2, o_Sat2_K, o_Sat2_M, t_enc_Sat2] = refineEncounter_cost(a_Sat2, ref2, oM, true, 0);
[cc3, o_Sat3_K, o_Sat3_M, t_enc_Sat3] = refineEncounter_cost(a_Sat3, ref3, oM, true, 0);
[cc4, o_Sat4_K, o_Sat4_M, t_enc_Sat4] = refineEncounter_cost(a_Sat4, ref4, oM, true, 0);

%% --- Tests closest passes -----------------------------------------------

load('Bashir_2_TargetsData.mat', 'oTM')
load('Bashir_2_04_TargetTrajectories.mat')

t0 = 3054600;

%{
[t1, d1, th1_1, th2_1] = closestPass([t0, t0 + 14400], o_Sat1_M, oTM(1), mu_M);
[t2, d2, th1_2, th2_2] = closestPass([t0, t0 + 14400], o_Sat2_M, oTM(2), mu_M);
[t3, d3, th1_3, th2_3] = closestPass([t0, t0 + 14400], o_Sat3_M, oTM(3), mu_M);
[t4, d4, th1_4, th2_4] = closestPass([t0, t0 + 14400], o_Sat4_M, oTM(4), mu_M);

[d1c, t_cp1, o1_pm1, o2_pm1, te_pm1] = encounterInSOICost([0; 0; 0], t0, ...
                            t0 + 14400, t_enc_Sat1, o_Sat1_K, o_Sat1_M, oM, oTM(1), mu_K, mu_M, rSOI_M);
[d2c, t_cp2, o1_pm2, o2_pm2, te_pm2] = encounterInSOICost([0; 0; 0], t0, ...
                            t0 + 14400, t_enc_Sat2, o_Sat2_K, o_Sat2_M, oM, oTM(2), mu_K, mu_M, rSOI_M);
[d3c, t_cp3, o1_pm3, o2_pm3, te_pm3] = encounterInSOICost([0; 0; 0], t0, ...
                            t0 + 14400, t_enc_Sat3, o_Sat3_K, o_Sat3_M, oM, oTM(3), mu_K, mu_M, rSOI_M);
[d4c, t_cp4, o1_pm4, o2_pm4, te_pm4] = encounterInSOICost([0; 0; 0], t0, ...
                            t0 + 14400, t_enc_Sat4, o_Sat4_K, o_Sat4_M, oM, oTM(4), mu_K, mu_M, rSOI_M);
%}

cf1 = @(v) encounterInSOICost(v, t0, ...
    t0 + 14400, t_enc_Sat1, o_Sat1_K, o_Sat1_M, oM, oTM(1), mu_K, mu_M, ...
    rSOI_M) ^ 2;

cf2 = @(v) encounterInSOICost(v, t0 + 300, ...
    t0 + 14400, t_enc_Sat2, o_Sat2_K, o_Sat2_M, oM, oTM(2), mu_K, mu_M, ...
    rSOI_M) ^ 2;

cf3 = @(v) encounterInSOICost(v, t0 + 600, ...
    t0 + 14400, t_enc_Sat3, o_Sat3_K, o_Sat3_M, oM, oTM(3), mu_K, mu_M, ...
    rSOI_M) ^ 2;

cf4 = @(v) encounterInSOICost(v, t0 + 900, ...
    t0 + 14400, t_enc_Sat4, o_Sat4_K, o_Sat4_M, oM, oTM(4), mu_K, mu_M, ...
    rSOI_M) ^ 2;

if optimize2(1, 1)
    fprintf('Designing manoeuvre for first satellite...\n')
    dv1 = optiGrad(cf1, dv1, 1e-9, 1e-5, 30);
end
if optimize2(1, 1) || ~exist('t_cp1', 'var') || ~exist('dv1_2', 'var') ...
        || ~exist('o_Sat1_K_pm', 'var') || ~exist('o_Sat1_M_pm', 'var')
    [~, t_cp1, o_Sat1_K_pm, o_Sat1_M_pm, te1_pm] = encounterInSOICost( ...
        dv1, t0, t0 + 14400, t_enc_Sat1, o_Sat1_K, o_Sat1_M, oM, ...
        oTM(1), mu_K, mu_M, rSOI_M);
    [xx1, v_cp1_1, ~] = structTimeToPosVelArg(t_cp1, o_Sat1_M_pm, mu_M);
    [~, v_cp1_2, ~] = structTimeToPosVelArg(t_cp1, oTM(1), mu_M);
    dv1_2 = refFrameLocal(xx1, v_cp1_1)' * (v_cp1_2 - v_cp1_1);
    save('Bashir_2_05_Sat1Data_02.mat', 'a_Sat1', 'dv1', 't_cp1', ...
        'dv1_2', 'o_Sat1_K_pm', 'o_Sat1_M_pm')
end

if optimize2(1, 2)
    fprintf('Designing manoeuvre for second satellite...\n')
    dv2 = optiGrad(cf2, dv2, 1e-9, 1e-6, 30);
end
if optimize2(1, 2) || ~exist('t_cp2', 'var') || ~exist('dv2_2', 'var') ...
        || ~exist('o_Sat2_K_pm', 'var') || ~exist('o_Sat2_M_pm', 'var')
    [~, t_cp2, o_Sat2_K_pm, o_Sat2_M_pm, te2_pm] = encounterInSOICost( ...
        dv2, t0 + 300, t0 + 14400, t_enc_Sat2, o_Sat2_K, o_Sat2_M, oM, ...
        oTM(2), mu_K, mu_M, rSOI_M);
    [xx2, v_cp2_1, ~] = structTimeToPosVelArg(t_cp2, o_Sat2_M_pm, mu_M);
    [~, v_cp2_2, ~] = structTimeToPosVelArg(t_cp2, oTM(2), mu_M);
    dv2_2 = refFrameLocal(xx2, v_cp2_1)' * (v_cp2_2 - v_cp2_1);
    save('Bashir_2_05_Sat2Data_02.mat', 'a_Sat2', 'dv2', 't_cp2', ...
        'dv2_2', 'o_Sat2_K_pm', 'o_Sat2_M_pm')
end

if optimize2(1, 3)
    fprintf('Designing manoeuvre for third satellite...\n')
    dv3 = optiGrad(cf3, dv3, 1e-9, 1e-6, 30);
end
if optimize2(1, 3) || ~exist('t_cp3', 'var') || ~exist('dv3_2', 'var') ...
        || ~exist('o_Sat3_K_pm', 'var') || ~exist('o_Sat3_M_pm', 'var')
    [~, t_cp3, o_Sat3_K_pm, o_Sat3_M_pm, te3_pm] = encounterInSOICost( ...
        dv3, t0 + 600, t0 + 14400, t_enc_Sat3, o_Sat3_K, o_Sat3_M, oM, ...
        oTM(3), mu_K, mu_M, rSOI_M);
    [xx3, v_cp3_1, ~] = structTimeToPosVelArg(t_cp3, o_Sat3_M_pm, mu_M);
    [~, v_cp3_2, ~] = structTimeToPosVelArg(t_cp3, oTM(3), mu_M);
    dv3_2 = refFrameLocal(xx3, v_cp3_1)' * (v_cp3_2 - v_cp3_1);
    save('Bashir_2_05_Sat3Data_02.mat', 'a_Sat3', 'dv3', 't_cp3', ...
        'dv3_2', 'o_Sat3_K_pm', 'o_Sat3_M_pm')
end

if optimize2(1, 4)
    fprintf('Designing manoeuvre for fourth satellite...\n')
    dv4 = optiGrad(cf4, dv4, 1e-9, 1e-6, 30);
end
if optimize2(1, 4) || ~exist('t_cp4', 'var') || ~exist('dv4_2', 'var') ...
        || ~exist('o_Sat4_K_pm', 'var') || ~exist('o_Sat4_M_pm', 'var')
    [~, t_cp4, o_Sat4_K_pm, o_Sat4_M_pm, te4_pm] = encounterInSOICost( ...
        dv4, t0 + 900, t0 + 14400, t_enc_Sat4, o_Sat4_K, o_Sat4_M, oM, ...
        oTM(4), mu_K, mu_M, rSOI_M);
    [xx4, v_cp4_1, ~] = structTimeToPosVelArg(t_cp4, o_Sat4_M_pm, mu_M);
    [~, v_cp4_2, ~] = structTimeToPosVelArg(t_cp4, oTM(4), mu_M);
    dv4_2 = refFrameLocal(xx4, v_cp4_1)' * (v_cp4_2 - v_cp4_1);
    save('Bashir_2_05_Sat4Data_02.mat', 'a_Sat4', 'dv4', 't_cp4', ...
        'dv4_2', 'o_Sat4_K_pm', 'o_Sat4_M_pm')
end

%% --- Output -------------------------------------------------------------

% ---- First satellite ----------------------------------------------------
fprintf('First satellite:\n')

fprintf('\n\tManoeuvre 1:\n')
fprintf('\t\tdv1 : %25.16f\n', dv1(1, 1))
fprintf('\t\tdv2 : %25.16f\n', dv1(2, 1))
fprintf('\t\tdv3 : %25.16f\n', dv1(3, 1))
fprintf('\t\t  t : %25.16f\n', t0)

[o_Sat1_M_pm_ap, o_Sat1_M_pm_pe] = peToApoPeri(o_Sat1_M_pm.p, o_Sat1_M_pm.e);
fprintf('\n\tPost-manoeuvre parameters (Minmus):\n')
fprintf('\t\t Apoapsis     :    %12.5e\n', o_Sat1_M_pm_ap - R_M)
fprintf('\t\t Periapsis    :    %12.5e\n', o_Sat1_M_pm_pe - R_M)
fprintf('\t\t Eccentricity : %10.4f\n', o_Sat1_M_pm.e)
fprintf('\t\t LAN          :   %5.1f\n', rToD * o_Sat1_M_pm.lan)
fprintf('\t\t Inclination  :   %5.1f\n', rToD * o_Sat1_M_pm.incl)
fprintf('\t\t Arg. pe.     :   %5.1f\n', rToD * o_Sat1_M_pm.aop)

fprintf('\n\tManoeuvre 2:\n')
fprintf('\t\tdv1 : %25.16f\n', dv1_2(1, 1))
fprintf('\t\tdv2 : %25.16f\n', dv1_2(2, 1))
fprintf('\t\tdv3 : %25.16f\n', dv1_2(3, 1))
fprintf('\t\t  t : %25.16f\n', t_cp1)

[oTM_1_ap, oTM_1_pe] = peToApoPeri(oTM(1).p, oTM(1).e);
fprintf('\n\tPost-manoeuvre parameters (Minmus):\n')
fprintf('\t\t Apoapsis     :    %12.5e\n', oTM_1_ap - R_M)
fprintf('\t\t Periapsis    :    %12.5e\n', oTM_1_pe - R_M)
fprintf('\t\t Eccentricity : %10.4f\n', oTM(1).e)
fprintf('\t\t LAN          :   %5.1f\n', rToD * oTM(1).lan)
fprintf('\t\t Inclination  :   %5.1f\n', rToD * oTM(1).incl)
fprintf('\t\t Arg. pe.     :   %5.1f\n', rToD * oTM(1).aop)


% ---- Second satellite ---------------------------------------------------
fprintf('\n\nSecond satellite:\n')

fprintf('\n\tManoeuvre 1:\n')
fprintf('\t\tdv1 : %25.16f\n', dv2(1, 1))
fprintf('\t\tdv2 : %25.16f\n', dv2(2, 1))
fprintf('\t\tdv3 : %25.16f\n', dv2(3, 1))
fprintf('\t\t  t : %25.16f\n', (t0 + 300))

[o_Sat2_M_pm_ap, o_Sat2_M_pm_pe] = peToApoPeri(o_Sat2_M_pm.p, o_Sat2_M_pm.e);
fprintf('\n\tPost-manoeuvre parameters (Minmus):\n')
fprintf('\t\t Apoapsis     :    %12.5e\n', o_Sat2_M_pm_ap - R_M)
fprintf('\t\t Periapsis    :    %12.5e\n', o_Sat2_M_pm_pe - R_M)
fprintf('\t\t Eccentricity : %10.4f\n', 0)
fprintf('\t\t Eccentricity : %10.4f\n', o_Sat2_M_pm.e)
fprintf('\t\t LAN          :   %5.1f\n', rToD * o_Sat2_M_pm.lan)
fprintf('\t\t Inclination  :   %5.1f\n', rToD * o_Sat2_M_pm.incl)
fprintf('\t\t Arg. pe.     :   %5.1f\n', rToD * o_Sat2_M_pm.aop)

fprintf('\n\tManoeuvre 2:\n')
fprintf('\t\tdv1 : %25.16f\n', dv2_2(1, 1))
fprintf('\t\tdv2 : %25.16f\n', dv2_2(2, 1))
fprintf('\t\tdv3 : %25.16f\n', dv2_2(3, 1))
fprintf('\t\t  t : %25.16f\n', t_cp2)

[oTM_2_ap, oTM_2_pe] = peToApoPeri(oTM(2).p, oTM(2).e);
fprintf('\n\tPost-manoeuvre parameters (Minmus):\n')
fprintf('\t\t Apoapsis     :    %12.5e\n', oTM_2_ap - R_M)
fprintf('\t\t Periapsis    :    %12.5e\n', oTM_2_pe - R_M)
fprintf('\t\t Eccentricity : %10.4f\n', oTM(2).e)
fprintf('\t\t LAN          :   %5.1f\n', rToD * oTM(2).lan)
fprintf('\t\t Inclination  :   %5.1f\n', rToD * oTM(2).incl)
fprintf('\t\t Arg. pe.     :   %5.1f\n', rToD * oTM(2).aop)


% ---- Third satellite ----------------------------------------------------
fprintf('\n\nThird satellite:\n')

fprintf('\n\tManoeuvre 1:\n')
fprintf('\t\tdv1 : %25.16f\n', dv3(1, 1))
fprintf('\t\tdv2 : %25.16f\n', dv3(2, 1))
fprintf('\t\tdv3 : %25.16f\n', dv3(3, 1))
fprintf('\t\t  t : %25.16f\n', (t0 + 600))

[o_Sat3_M_pm_ap, o_Sat3_M_pm_pe] = peToApoPeri(o_Sat3_M_pm.p, o_Sat3_M_pm.e);
fprintf('\n\tPost-manoeuvre parameters (Minmus):\n')
fprintf('\t\t Apoapsis     :    %12.5e\n', o_Sat3_M_pm_ap - R_M)
fprintf('\t\t Periapsis    :    %12.5e\n', o_Sat3_M_pm_pe - R_M)
fprintf('\t\t Eccentricity : %10.4f\n', o_Sat3_M_pm.e)
fprintf('\t\t LAN          :   %5.1f\n', rToD * o_Sat3_M_pm.lan)
fprintf('\t\t Inclination  :   %5.1f\n', rToD * o_Sat3_M_pm.incl)
fprintf('\t\t Arg. pe.     :   %5.1f\n', rToD * o_Sat3_M_pm.aop)

fprintf('\n\tManoeuvre 2:\n')
fprintf('\t\tdv1 : %25.16f\n', dv3_2(1, 1))
fprintf('\t\tdv2 : %25.16f\n', dv3_2(2, 1))
fprintf('\t\tdv3 : %25.16f\n', dv3_2(3, 1))
fprintf('\t\t  t : %25.16f\n', t_cp3)

[oTM_3_ap, oTM_3_pe] = peToApoPeri(oTM(3).p, oTM(3).e);
fprintf('\n\tPost-manoeuvre parameters (Minmus):\n')
fprintf('\t\t Apoapsis     :    %12.5e\n', oTM_3_ap - R_M)
fprintf('\t\t Periapsis    :    %12.5e\n', oTM_3_pe - R_M)
fprintf('\t\t Eccentricity : %10.4f\n', oTM(3).e)
fprintf('\t\t LAN          :   %5.1f\n', rToD * oTM(3).lan)
fprintf('\t\t Inclination  :   %5.1f\n', rToD * oTM(3).incl)
fprintf('\t\t Arg. pe.     :   %5.1f\n', rToD * oTM(3).aop)


% ---- Fourth satellite ---------------------------------------------------
fprintf('\n\nFourth satellite:\n')

fprintf('\n\tManoeuvre 1:\n')
fprintf('\t\tdv1 : %25.16f\n', dv4(1, 1))
fprintf('\t\tdv2 : %25.16f\n', dv4(2, 1))
fprintf('\t\tdv3 : %25.16f\n', dv4(3, 1))
fprintf('\t\t  t : %25.16f\n', (t0 + 900))

[o_Sat4_M_pm_ap, o_Sat4_M_pm_pe] = peToApoPeri(o_Sat4_M_pm.p, o_Sat4_M_pm.e);
fprintf('\n\tPost-manoeuvre parameters (Minmus):\n')
fprintf('\t\t Apoapsis     :    %12.5e\n', o_Sat4_M_pm_ap - R_M)
fprintf('\t\t Periapsis    :    %12.5e\n', o_Sat4_M_pm_pe - R_M)
fprintf('\t\t Eccentricity : %10.4f\n', o_Sat4_M_pm.e)
fprintf('\t\t LAN          :   %5.1f\n', rToD * o_Sat4_M_pm.lan)
fprintf('\t\t Inclination  :   %5.1f\n', rToD * o_Sat4_M_pm.incl)
fprintf('\t\t Arg. pe.     :   %5.1f\n', rToD * o_Sat4_M_pm.aop)

fprintf('\n\tManoeuvre 2:\n')
fprintf('\t\tdv1 : %25.16f\n', dv4_2(1, 1))
fprintf('\t\tdv2 : %25.16f\n', dv4_2(2, 1))
fprintf('\t\tdv3 : %25.16f\n', dv4_2(3, 1))
fprintf('\t\t  t : %25.16f\n', t_cp4)

[oTM_4_ap, oTM_4_pe] = peToApoPeri(oTM(4).p, oTM(4).e);
fprintf('\n\tPost-manoeuvre parameters (Minmus):\n')
fprintf('\t\t Apoapsis     :    %12.5e\n', oTM_4_ap - R_M)
fprintf('\t\t Periapsis    :    %12.5e\n', oTM_4_pe - R_M)
fprintf('\t\t Eccentricity : %10.4f\n', oTM(4).e)
fprintf('\t\t LAN          :   %5.1f\n', rToD * oTM(4).lan)
fprintf('\t\t Inclination  :   %5.1f\n', rToD * oTM(4).incl)
fprintf('\t\t Arg. pe.     :   %5.1f\n', rToD * oTM(4).aop)


% -------------------------------------------------------------------------