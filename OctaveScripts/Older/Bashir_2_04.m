fclose all;
close all
format long g
more off
clear
clc

vstr = version;
vnum = str2double(vstr(1, (end - 5):(end - 2)));

load('Bashir_2_OrbitData.mat');
load('Bashir_2_TargetsData.mat');

SOIr_K = 84159286;
N = 251;

% Initial time
t0 = 3016000;

% Initial argument
thiK = timeStructToArg(t0, oOK, mu_K);

% Exit argument
thOM = -thEM;

% Exit time
oT   = argStructToTime(thOM, oOM, mu_M);

% Minmus position and velocity at exit
[xxOMK, vvOMK] = structTimeToPosVelArg(oT, oMK, mu_K);

% Object position and velocity relative to Minmus at exit
[xxOOM, vvOOM] = structArgToPosVelTime(thOM, oOM, mu_M);

% Minmus crash argument
mKa = -acos((oOM.p / R_K - 1) / (oOM.e));
tC = argStructToTime(mKa, oOM, mu_M);

% Object position and velocity relative to Kerbin at exit
xxOOK = xxOOM + xxOMK;
vvOOK = vvOOM + vvOMK;

% New orbit around Kerbin after exit
[thOK, oOK2] = posVelToStruct(xxOOK, vvOOK, mu_K, oT);

% Kerbin exit argument for new orbit
thFK2 = acos((oOK2.p / SOIr_K - 1) / (oOK2.e));

% Kerbin exit time for new orbit
fT = argStructToTime(thFK2, oOK2, mu_K);

%
eTi = eT - 10 * 3600;
oTf = oT + 10 * 3600;

% Builds Bashir 2 default trajectory (through Minmus...)
thB1 = linspace(thiK, thEK,  N);
thB2 = linspace(thEM, thOM,  N);
thB3 = linspace(thOK, thFK2, N);
[XXB1, ~] = structArgToPosVelTime(thB1, oOK, mu_K);
[XXB2M, ~] = structArgToPosVelTime(thB2, oOM, mu_M);
ttT = linspace(eTi, oTf, N);
[XXM2K, ~, thM2B] = structTimeToPosVelArg(linspace(eT, oT, N), oMK, mu_K);
[XXM3K, ~, thM2M] = structTimeToPosVelArg(ttT, oMK, mu_K);
XXB2 = XXB2M + XXM2K;
[XXB3, ~] = structArgToPosVelTime(thB3, oOK2, mu_K);

% Target orbits around Minmus
[XXT1M, ~] = structTimeToPosVelArg(ttT, oTM(1), mu_M);
XXT1K = XXT1M + XXM3K;
[XXT2M, ~] = structTimeToPosVelArg(ttT, oTM(2), mu_M);
XXT2K = XXT2M + XXM3K;
[XXT3M, ~] = structTimeToPosVelArg(ttT, oTM(3), mu_M);
XXT3K = XXT3M + XXM3K;
[XXT4M, ~] = structTimeToPosVelArg(ttT, oTM(4), mu_M);
XXT4K = XXT4M + XXM3K;

m1.t = t0;
m1.v = [0; 0; 0];
m2.t = t0 + 2 * 180;
m2.v = [0; 0; 0];
m3.t = t0 + 4 * 180;
m3.v = [0; 0; 0];
m4.t = t0 + 6 * 180;
m4.v = [0; 0; 0];
x0   = [0; 0; 0];

if ~exist('Bashir2_03_OptData2.mat', 'file')
    fprintf('Target 1...\n')
    cf = @(a) Bashir_2_03_Cost01(a, m1, eT, mu_K, mu_M, ...
        oOK, tC, SOIr_M, oMK, R_M, oTM(1));
    [x1] = optiGrad(cf, x0, 1e-4, 1-6, 15);
    fprintf('Target 2...\n')
    cf = @(a) Bashir_2_03_Cost01(a, m2, eT, mu_K, mu_M, ...
        oOK, tC, SOIr_M, oMK, R_M, oTM(2));
    [x2] = optiGrad(cf, x0, 1e-6, 1-6, 15);
    fprintf('Target 3...\n')
    cf = @(a) Bashir_2_03_Cost01(a, m3, eT, mu_K, mu_M, ....
        oOK, tC, SOIr_M, oMK, R_M, oTM(3));
    [x3] = optiGrad(cf, x0, 1e-4, 1-6, 15);
    fprintf('Target 4...\n')
    cf = @(a) Bashir_2_03_Cost01(a, m4, eT, mu_K, mu_M, ...
        oOK, tC, SOIr_M, oMK, R_M, oTM(4));
    [x4] = optiGrad(cf, x0, 1e-4, 1-6, 15);
    save Bashir2_03_OptData2.mat x1 x2 x3 x4
else
    load Bashir2_03_OptData2.mat
    if ~exist('x1', 'var')
        fprintf('Target 1...\n')
        cf = @(a) Bashir_2_03_Cost01(a, m1, eT, mu_K, mu_M, ...
            oOK, tC, SOIr_M, oMK, R_M, oTM(1));
        [x1] = optiGrad(cf, x0, 1e-4, 1-6, 15);
    end
    if ~exist('x2', 'var')
        fprintf('Target 2...\n')
        x02 = [          15.2297160960409
          22.1490941236381
         -23.4070719830622];
        cf = @(a) Bashir_2_03_Cost01(a, m2, eT, mu_K, mu_M, ...
            oOK, tC, SOIr_M, oMK, R_M, oTM(2));
        [x2] = optiGrad(cf, x02, 1e-9, 1-6, 15);
    end
    if ~exist('x3', 'var')
            fprintf('Target 3...\n')
        cf = @(a) Bashir_2_03_Cost01(a, m3, eT, mu_K, mu_M, ....
            oOK, tC, SOIr_M, oMK, R_M, oTM(3));
        [x3] = optiGrad(cf, x0, 1e-4, 1-6, 15);
    end
    if ~exist('x4', 'var')
        cf = @(a) Bashir_2_03_Cost01(a, m4, eT, mu_K, mu_M, ...
            oOK, tC, SOIr_M, oMK, R_M, oTM(4));
        [x4] = optiGrad(cf, x0, 1e-4, 1-6, 15);
    end
    save Bashir2_03_OptData2.mat x1 x2 x3 x4
end

[cc1, tcp1, dv1] = Bashir_2_03_Cost01(x1, m1, eT, mu_K, mu_M, ...
    oOK, tC, SOIr_M, oMK, R_M, oTM(1));
[cc2, tcp2, dv2] = Bashir_2_03_Cost01(x2, m2, eT, mu_K, mu_M, ...
    oOK, tC, SOIr_M, oMK, R_M, oTM(2));
[cc3, tcp3, dv3] = Bashir_2_03_Cost01(x3, m3, eT, mu_K, mu_M, ...
    oOK, tC, SOIr_M, oMK, R_M, oTM(3));
[cc4, tcp4, dv4] = Bashir_2_03_Cost01(x4, m4, eT, mu_K, mu_M, ...
    oOK, tC, SOIr_M, oMK, R_M, oTM(4));

m1.v = x1;
[of1K] = applyManoeuvre(oOK, m1, mu_K);
[t1o, ~, ~, ~] = entryInSphere([m1.t, oT], SOIr_M, of1K, oMK, mu_K);
[x1OK, v1OK, ~] = structTimeToPosVelArg(t1o, of1K, mu_K);
[x1MK, v1MK, ~] = structTimeToPosVelArg(t1o, oMK, mu_K);
x1OM = x1OK - x1MK;
v1OM = v1OK - v1MK;
[~, of1M] = posVelToStruct(x1OM, v1OM, mu_M, t1o, 1);
[xx1OK, ~, ~] = structTimeToPosVelArg(linspace(m1.t, t1o, N), of1K, mu_K);
[xx1OKp, ~, ~] = structTimeToPosVelArg(linspace(m1.t, t1o, 2), of1K, mu_K);
tt1 = linspace(t1o, tcp1, N);
[xx1OM, ~, ~] = structTimeToPosVelArg(tt1, of1M, mu_M);
[xx1OMp, ~, ~] = structTimeToPosVelArg(tcp1, of1M, mu_M);
[xx1MK, ~, ~] = structTimeToPosVelArg(tt1, oMK, mu_K);
[xx1MKp, ~, ~] = structTimeToPosVelArg(tcp1, oMK, mu_K);
xx1OK = [xx1OK, (xx1OM + xx1MK)];
xx1OKp = [xx1OKp, (xx1OMp + xx1MKp)];
m1_2.t = tcp1;
m1_2.v = dv1;
[of1M2] = applyManoeuvre(of1M, m1_2, mu_M);
tt12 = linspace(tcp1, oTf, N);
[xx1OM2, ~, ~] = structTimeToPosVelArg(tt12, of1M2, mu_M);
[xx1OMp2, ~, ~] = structTimeToPosVelArg(oTf, of1M2, mu_M);
[xx1MK2, ~, ~] = structTimeToPosVelArg(tt12, oMK, mu_K);
[xx1MKp2, ~, ~] = structTimeToPosVelArg(oTf, oMK, mu_K);
xx1OK = [xx1OK, (xx1OM2 + xx1MK2)];
xx1OKp = [xx1OKp, (xx1OMp2 + xx1MKp2)];

m2.v = x2;
[of2K] = applyManoeuvre(oOK, m2, mu_K);
[t2o, ~, ~, ~] = entryInSphere([m2.t, oT], SOIr_M, of2K, oMK, mu_K);
[x2OK, v2OK, ~] = structTimeToPosVelArg(t2o, of2K, mu_K);
[x2MK, v2MK, ~] = structTimeToPosVelArg(t2o, oMK, mu_K);
x2OM = x2OK - x2MK;
v2OM = v2OK - v2MK;
[~, of2M] = posVelToStruct(x2OM, v2OM, mu_M, t2o, 1);
[xx2OK, ~, ~] = structTimeToPosVelArg(linspace(m2.t, t2o, N), of2K, mu_K);
[xx2OKp, ~, ~] = structTimeToPosVelArg(linspace(m2.t, t2o, 2), of2K, mu_K);
tt2 = linspace(t2o, tcp2, N);
[xx2OM, ~, ~] = structTimeToPosVelArg(tt2, of2M, mu_M);
[xx2OMp, ~, ~] = structTimeToPosVelArg(tcp2, of2M, mu_M);
[xx2MK, ~, ~] = structTimeToPosVelArg(tt2, oMK, mu_K);
[xx2MKp, ~, ~] = structTimeToPosVelArg(tcp2, oMK, mu_K);
xx2OK = [xx2OK, (xx2OM + xx2MK)];
xx2OKp = [xx2OKp, (xx2OMp + xx2MKp)];
m2_2.t = tcp2;
m2_2.v = dv2;
[of2M2] = applyManoeuvre(of2M, m2_2, mu_M);
tt22 = linspace(tcp2, oTf, N);
[xx2OM2, ~, ~] = structTimeToPosVelArg(tt22, of2M2, mu_M);
[xx2OMp2, ~, ~] = structTimeToPosVelArg(oTf, of2M2, mu_M);
[xx2MK2, ~, ~] = structTimeToPosVelArg(tt22, oMK, mu_K);
[xx2MKp2, ~, ~] = structTimeToPosVelArg(oTf, oMK, mu_K);
xx2OK = [xx2OK, (xx2OM2 + xx2MK2)];
xx2OKp = [xx2OKp, (xx2OMp2 + xx2MKp2)];

m3.v = x3;
[of3K] = applyManoeuvre(oOK, m3, mu_K);
[t3o, ~, ~, ~] = entryInSphere([m3.t, oT], SOIr_M, of3K, oMK, mu_K);
[x3OK, v3OK, ~] = structTimeToPosVelArg(t3o, of3K, mu_K);
[x3MK, v3MK, ~] = structTimeToPosVelArg(t3o, oMK, mu_K);
x3OM = x3OK - x3MK;
v3OM = v3OK - v3MK;
[~, of3M] = posVelToStruct(x3OM, v3OM, mu_M, t3o, 1);
[xx3OK, ~, ~] = structTimeToPosVelArg(linspace(m3.t, t3o, N), of3K, mu_K);
[xx3OKp, ~, ~] = structTimeToPosVelArg(linspace(m3.t, t3o, 2), of3K, mu_K);
tt3 = linspace(t3o, tcp3, N);
[xx3OM, ~, ~] = structTimeToPosVelArg(tt3, of3M, mu_M);
[xx3OMp, ~, ~] = structTimeToPosVelArg(tcp3, of3M, mu_M);
[xx3MK, ~, ~] = structTimeToPosVelArg(tt3, oMK, mu_K);
[xx3MKp, ~, ~] = structTimeToPosVelArg(tcp3, oMK, mu_K);
xx3OK = [xx3OK, (xx3OM + xx3MK)];
xx3OKp = [xx3OKp, (xx3OMp + xx3MKp)];
m3_2.t = tcp3;
m3_2.v = dv3;
[of3M2] = applyManoeuvre(of3M, m3_2, mu_M);
tt32 = linspace(tcp3, oTf, N);
[xx3OM2, ~, ~] = structTimeToPosVelArg(tt32, of3M2, mu_M);
[xx3OMp2, ~, ~] = structTimeToPosVelArg(oTf, of3M2, mu_M);
[xx3MK2, ~, ~] = structTimeToPosVelArg(tt32, oMK, mu_K);
[xx3MKp2, ~, ~] = structTimeToPosVelArg(oTf, oMK, mu_K);
xx3OK = [xx3OK, (xx3OM2 + xx3MK2)];
xx3OKp = [xx3OKp, (xx3OMp2 + xx3MKp2)];

m4.v = x4;
[of4K] = applyManoeuvre(oOK, m4, mu_K);
[t4o, ~, ~, ~] = entryInSphere([m4.t, oT], SOIr_M, of4K, oMK, mu_K);
[x4OK, v4OK, ~] = structTimeToPosVelArg(t4o, of4K, mu_K);
[x4MK, v4MK, ~] = structTimeToPosVelArg(t4o, oMK, mu_K);
x4OM = x4OK - x4MK;
v4OM = v4OK - v4MK;
[~, of4M] = posVelToStruct(x4OM, v4OM, mu_M, t4o, 1);
[xx4OK, ~, ~] = structTimeToPosVelArg(linspace(m4.t, t4o, N), of4K, mu_K);
[xx4OKp, ~, ~] = structTimeToPosVelArg(linspace(m4.t, t4o, 2), of4K, mu_K);
tt4 = linspace(t4o, tcp4, N);
[xx4OM, ~, ~] = structTimeToPosVelArg(tt4, of4M, mu_M);
[xx4OMp, ~, ~] = structTimeToPosVelArg(tcp4, of4M, mu_M);
[xx4MK, ~, ~] = structTimeToPosVelArg(tt4, oMK, mu_K);
[xx4MKp, ~, ~] = structTimeToPosVelArg(tcp4, oMK, mu_K);
xx4OK = [xx4OK, (xx4OM + xx4MK)];
xx4OKp = [xx4OKp, (xx4OMp + xx4MKp)];
m4_2.t = tcp4;
m4_2.v = dv4;
[of4M2] = applyManoeuvre(of4M, m4_2, mu_M);
tt42 = linspace(tcp4, oTf, N);
[xx4OM2, ~, ~] = structTimeToPosVelArg(tt42, of4M2, mu_M);
[xx4OMp2, ~, ~] = structTimeToPosVelArg(oTf, of4M2, mu_M);
[xx4MK2, ~, ~] = structTimeToPosVelArg(tt42, oMK, mu_K);
[xx4MKp2, ~, ~] = structTimeToPosVelArg(oTf, oMK, mu_K);
xx4OK = [xx4OK, (xx4OM2 + xx4MK2)];
xx4OKp = [xx4OKp, (xx4OMp2 + xx4MKp2)];

fprintf('First relay:\n')
fprintf('\tManoeuvre 1:\n')
fprintf('\t\tv1:\t%11.3f\n', m1.v(1, 1))
fprintf('\t\tv2:\t%11.3f\n', m1.v(2, 1))
fprintf('\t\tv3:\t%11.3f\n', m1.v(3, 1))
fprintf('\t\tt: \t%11.3f\n', m1.t)
fprintf('\tManoeuvre 2:\n')
fprintf('\t\tv1:\t%11.3f\n', m1_2.v(1, 1))
fprintf('\t\tv2:\t%11.3f\n', m1_2.v(2, 1))
fprintf('\t\tv3:\t%11.3f\n', m1_2.v(3, 1))
fprintf('\t\tt: \t%11.3f\n', m1_2.t)
fprintf('\nSecond relay:\n')
fprintf('\tManoeuvre 1:\n')
fprintf('\t\tv1:\t%11.3f\n', m2.v(1, 1))
fprintf('\t\tv2:\t%11.3f\n', m2.v(2, 1))
fprintf('\t\tv3:\t%11.3f\n', m2.v(3, 1))
fprintf('\t\tt: \t%11.3f\n', m2.t)
fprintf('\tManoeuvre 2:\n')
fprintf('\t\tv1:\t%11.3f\n', m2_2.v(1, 1))
fprintf('\t\tv2:\t%11.3f\n', m2_2.v(2, 1))
fprintf('\t\tv3:\t%11.3f\n', m2_2.v(3, 1))
fprintf('\t\tt: \t%11.3f\n', m2_2.t)
fprintf('\nThird relay:\n')
fprintf('\tManoeuvre 1:\n')
fprintf('\t\tv1:\t%11.3f\n', m3.v(1, 1))
fprintf('\t\tv2:\t%11.3f\n', m3.v(2, 1))
fprintf('\t\tv3:\t%11.3f\n', m3.v(3, 1))
fprintf('\t\tt: \t%11.3f\n', m3.t)
fprintf('\tManoeuvre 2:\n')
fprintf('\t\tv1:\t%11.3f\n', m3_2.v(1, 1))
fprintf('\t\tv2:\t%11.3f\n', m3_2.v(2, 1))
fprintf('\t\tv3:\t%11.3f\n', m3_2.v(3, 1))
fprintf('\t\tt: \t%11.3f\n', m3_2.t)
fprintf('\nFourth relay:\n')
fprintf('\tManoeuvre 1:\n')
fprintf('\t\tv1:\t%11.3f\n', m4.v(1, 1))
fprintf('\t\tv2:\t%11.3f\n', m4.v(2, 1))
fprintf('\t\tv3:\t%11.3f\n', m4.v(3, 1))
fprintf('\t\tt: \t%11.3f\n', m4.t)
fprintf('\tManoeuvre 2:\n')
fprintf('\t\tv1:\t%11.3f\n', m4_2.v(1, 1))
fprintf('\t\tv2:\t%11.3f\n', m4_2.v(2, 1))
fprintf('\t\tv3:\t%11.3f\n', m4_2.v(3, 1))
fprintf('\t\tt: \t%11.3f\n', m4_2.t)

%% ------------------------------------------------------------------------
%  Graphic output

figure(1)
hold on
if vnum >= 2019
    addToolbarExplorationButtons(gcf)
end
axis equal
plot3(0, 0, 0, 'bo', 'LineWidth', 2)
plot3(XXM3K(1, :), XXM3K(2, :), XXM3K(3, :), 'color', [0, 1, 0], 'LineWidth', 2)
plot3(XXB1(1, :), XXB1(2, :), XXB1(3, :), 'color', [0.75, 0, 0])
plot3(XXB1(1, 1), XXB1(2, 1), XXB1(3, 1), '*', 'color', [0.75, 0, 0])
plot3(XXB1(1, end), XXB1(2, end), XXB1(3, end), '*', 'color', [0.75, 0, 0])
plot3(XXB2(1, :), XXB2(2, :), XXB2(3, :), 'color', [1, 0, 0])
plot3(XXB2(1, 1), XXB2(2, 1), XXB2(3, 1), 'o', 'color', [1, 0, 0])
plot3(XXB2(1, end), XXB2(2, end), XXB2(3, end), 'o', 'color', [1, 0, 0])
plot3(XXB3(1, :), XXB3(2, :), XXB3(3, :), 'color', [1, 0.25, 0.25])
plot3(XXB3(1, 1), XXB3(2, 1), XXB3(3, 1), '*', 'color', [1, 0.25, 0.25])
plot3(XXB3(1, end), XXB3(2, end), XXB3(3, end), '*', 'color', [1, 0.25, 0.25])

plot3(XXT1K(1, :), XXT1K(2, :), XXT1K(3, :), 'color', [1, 1, 0])
plot3(XXT2K(1, :), XXT2K(2, :), XXT2K(3, :), 'color', [1, 0, 1])
plot3(XXT3K(1, :), XXT3K(2, :), XXT3K(3, :), 'color', [0, 1, 1])
plot3(XXT4K(1, :), XXT4K(2, :), XXT4K(3, :), 'color', [1, 0.5, 0])

plot3(xx1OK(1, :), xx1OK(2, :), xx1OK(3, :), 'color', [0.75, 0.75, 0])
plot3(xx1OKp(1, :), xx1OKp(2, :), xx1OKp(3, :), 'o', 'color', [0.5, 0.5, 0])
plot3(xx2OK(1, :), xx2OK(2, :), xx2OK(3, :), 'color', [0.75, 0, 0.75])
plot3(xx2OKp(1, :), xx2OKp(2, :), xx2OKp(3, :), 'o', 'color', [0.5, 0, 0.5])
plot3(xx3OK(1, :), xx3OK(2, :), xx3OK(3, :), 'color', [0, 0.75, 0.75])
plot3(xx3OKp(1, :), xx3OKp(2, :), xx3OKp(3, :), 'o', 'color', [0, 0.5, 0.5])
plot3(xx4OK(1, :), xx4OK(2, :), xx4OK(3, :), 'color', [0.75, 0.375, 0])
plot3(xx4OKp(1, :), xx4OKp(2, :), xx4OKp(3, :), 'o', 'color', [0.5, 0.25, 0])
hold off

%% ---

save('Bashir_2_04_TargetTrajectories.mat', 'of1K', 'of1M', 'of2K', 'of2M', 'of3K', 'of3M', 'of4K', 'of4M')
