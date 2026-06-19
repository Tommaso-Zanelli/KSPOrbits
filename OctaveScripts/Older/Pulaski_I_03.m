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

%% - Manoeuvre design, fourth (?) attempt ---------------------------------

T = 293.8;
A = 2.79675370863456417476982;

%ttm_offs = ttm_offs + T;

% Defines initial orbit object
o0 = paramsToStruct(p0, e0, lan0, incl0, aop0, pe0T);

a = [A * ones(1, (Nm1 + Nm2)), zeros(1, (Nm1 + Nm2))];
cf = @(a, T) Pulaski_I_03_f01(o0, ttm_offs + T, dvm, Nm1, Nm2, mu_K, 0, 0, a(1, 1:(Nm1 + Nm2)), a(1, (Nm1 + Nm2 + 1):(2 * (Nm1 + Nm2))), 675000);

% Gradm

da = 2^-25;

%a0 = [0.76910994764397909529, 2.79515399943993260000, 2.79532975356887410000, 2.79550369961216380000, 2.79530149346537900000, 2.79568858863589710000, 2.79608809116235160000, 2.79648434996469010000, 2.79686054056215870000, 2.79723597297779620000, 2.79757957043331950000, 2.79788006531175170000, 2.79814658557853810000, 2.79832724567565230000, 2.79842787229498220000, 2.79841235491330000000, 2.79827049378983620000, 2.79797098270853790000, 2.79749215430761080000, -0.03266484261183293100, -0.03192395389637653300, -0.03116981321578796500, -0.04704933335051102900, -0.04542530120473203800, -0.04373107222260826800, -0.04194460479927680800, -0.04004342314700173900, -0.03799972209575697800, -0.03578508652451663200, -0.03336091328757122300, -0.03068939230380004400, -0.02772061934249517300, -0.02440004091094999400, -0.02066725045802538300, -0.01645657642837961400, -0.01170960468758313800, -0.00640001914888779780];
%a0 = [7.6912634931603141e-01, 2.7951442274140526e+00, 2.7953201159360765e+00, 2.7954942056131244e+00, 2.7952869672494418e+00, 2.7956743756614166e+00, 2.7960742298318517e+00, 2.7964708832531628e+00, 2.7968475170251792e+00, 2.7972234481007310e+00, 2.7975676069352375e+00, 2.7978687346586582e+00, 2.7981359698739841e+00, 2.7983174381768472e+00, 2.7984189797754411e+00, 2.7984044988401315e+00, 2.7982638112578644e+00, 2.7979656240119857e+00, 2.7974882738037952e+00, -3.2670316001790038e-02, -3.1929340788156256e-02, -3.1175111732548215e-02, -4.7057507327045521e-02, -4.5433282491463750e-02, -4.3738849539934330e-02, -4.1952163071533773e-02, -4.0050743326413660e-02, -3.8006779969311698e-02, -3.5791852207260305e-02, -3.3367349992186775e-02, -3.0695455253033112e-02, -2.7726254017821871e-02, -2.4405181372140764e-02, -2.0671816887786774e-02, -1.6460472783953722e-02, -1.1712717481771629e-02, -6.4022182292352032e-03];
%a0 = [7.8095946427693219e-01, 2.7862473325531267e+00, 2.7866899987730487e+00, 2.7871497420392228e+00, 2.7829567974815648e+00, 2.7839751810469116e+00, 2.7850383197398911e+00, 2.7861431463812738e+00, 2.7872612559936862e+00, 2.7884203448102558e+00, 2.7895889858576695e+00, 2.7907520168320876e+00, 2.7919196755401270e+00, 2.7930184702497209e+00, 2.7940668593056097e+00, 2.7950112277150767e+00, 2.7958154635907255e+00, 2.7964270315810982e+00, 2.7967701361963653e+00, -3.6687824307318311e-02, -3.5866067163835010e-02, -3.5026770733526086e-02, -5.2849236319596163e-02, -5.1043055766062884e-02, -4.9160887308501378e-02, -4.7170012056366700e-02, -4.5045850318214134e-02, -4.2755616909519532e-02, -4.0276609899012980e-02, -3.7556674877427448e-02, -3.4558990911250295e-02, -3.1221145017892714e-02, -2.7490617981556675e-02, -2.3288148147900905e-02, -1.8545463627954336e-02, -1.3202157117743382e-02, -7.2161272705254867e-03];
%a0 = [7.9325291525392694e-01, 2.7772524251857624e+00, 2.7779623499944339e+00, 2.7787107736595904e+00, 2.7704815100631355e+00, 2.7721316639822264e+00, 2.7738594130649012e+00, 2.7756725691442021e+00, 2.7775409488145040e+00, 2.7794904175353556e+00, 2.7814892669139071e+00, 2.7835257816145962e+00, 2.7855973575503019e+00, 2.7876341293936910e+00, 2.7896426895957789e+00, 2.7915602922536706e+00, 2.7933248530631389e+00, 2.7948615197320295e+00, 2.7960392007499362e+00, -4.0921982913155401e-02, -4.0011903493825642e-02, -3.9081190053555814e-02, -5.8947549524176024e-02, -5.6948985601348902e-02, -5.4866587885049292e-02, -5.2664329601807822e-02, -5.0309651317087610e-02, -4.7765134791400862e-02, -4.5006450512841328e-02, -4.1981193133771133e-02, -3.8642796222087980e-02, -3.4921153314602126e-02, -3.0755781111153561e-02, -2.6058631182868212e-02, -2.0755096018116214e-02, -1.4778924401635644e-02, -8.0807903123033646e-03];
%a0 = [8.0605365169853072e-01, 2.7681510591597793e+00, 2.7691359096145023e+00, 2.7701732132174999e+00, 2.7578666180408495e+00, 2.7601442974977193e+00, 2.7625358366624493e+00, 2.7650595472566453e+00, 2.7676803699424402e+00, 2.7704272419786724e+00, 2.7732640701664515e+00, 2.7761840097873818e+00, 2.7791662190093711e+00, 2.7821561222380278e+00, 2.7851445501220544e+00, 2.7880481840900098e+00, 2.7907900728090231e+00, 2.7932689210323849e+00, 2.7952935586185590e+00, -4.5395104224253358e-02, -4.4390044957780039e-02, -4.3365571817751809e-02, -6.5378663596694295e-02, -6.3183869530467673e-02, -6.0900239500799591e-02, -5.8470522112462980e-02, -5.5875086028086655e-02, -5.3068510937281807e-02, -5.0019355444069019e-02, -4.6667797112708236e-02, -4.2970731835631083e-02, -3.8841134942980908e-02, -3.4217614051834246e-02, -2.9000000417039196e-02, -2.3102571256812678e-02, -1.6456492372948210e-02, -9.0003870535238176e-03];
%a0 = [8.9005069850136820e-01,  2.7151753771058447e+00,  2.7176809890190659e+00,  2.7203342170336455e+00,  2.6846181513415948e+00,  2.6903579092298520e+00,  2.6964799326999978e+00,  2.7030437243448646e+00,  2.7099863937541508e+00,  2.7173375638753461e+00,  2.7250626401110010e+00,  2.7331290566712534e+00,  2.7414894252526536e+00,  2.7500649738508867e+00,  2.7587924336143539e+00,  2.7675257224797760e+00,  2.7760212670984030e+00,  2.7840035097436742e+00,  2.7909736831685343e+00, -7.6396624213930689e-02, -7.4803949869812189e-02, -7.3177070725171101e-02, -1.0993630831364457e-01, -1.0645522642834421e-01, -1.0282796281804532e-01, -9.8935829574330447e-02, -9.4763848343298471e-02, -9.0221326603961260e-02, -8.5243739270560670e-02, -7.9733424499890043e-02, -7.3590363812967144e-02, -6.6686140792760376e-02, -5.8880956484678927e-02, -5.0018054617874931e-02, -3.9942483696825806e-02, -2.8495621181030679e-02, -1.5610469649286928e-02];
%a0 = [ 1.0030801900684871e+00,  2.6579092494151357e+00,  2.6624821989284806e+00,  2.6672563164780971e+00,  2.6096396844050740e+00,  2.6198211048725430e+00,  2.6307347959396283e+00,  2.6423935619333183e+00,  2.6547492200668934e+00,  2.6677613396586768e+00,  2.6813302618644168e+00,  2.6953704738756810e+00,  2.7097820291187822e+00,  2.7242854619644064e+00,  2.7387509717027334e+00,  2.7528195872842232e+00,  2.7661816460612383e+00,  2.7783568522168998e+00,  2.7886752740939778e+00, -1.0438262760647268e-01, -1.0406404237416024e-01, -1.0367456004362459e-01, -1.5880785157978167e-01, -1.5820338471425122e-01, -1.5739574353731073e-01, -1.5618957398273603e-01, -1.5454731200015828e-01, -1.5221216638235871e-01, -1.4901223598651595e-01, -1.4466659061504406e-01, -1.3880925600526503e-01, -1.3092588758515494e-01, -1.2033286265916802e-01, -1.0618569450120631e-01, -8.7603823523659641e-02, -6.3940714197028709e-02, -3.5441238385665444e-02];
%a0 = [1.1935780385124302e+00, 2.6032810833041156e+00, 2.6111657346498922e+00, 2.6195120791273285e+00, 2.5573640504706163e+00, 2.5738204693121332e+00, 2.5915548033337417e+00, 2.6105208463284413e+00, 2.6304978053043055e+00, 2.6512863517850542e+00, 2.6725203497317174e+00, 2.6938807712445327e+00, 2.7149178317145686e+00, 2.7349072313883385e+00, 2.7532751146834809e+00, 2.7692775363594988e+00, 2.7822400456392691e+00, 2.7916586582185285e+00, 2.7969416251915726e+00, -1.7228044178985968e-01, -1.7489655042898283e-01, -1.7739983071487817e-01, -2.6765375789578233e-01, -2.7365300344503202e-01, -2.7936501898591265e-01, -2.8444307891919263e-01, -2.8863388692441067e-01, -2.9145709303482875e-01, -2.9254515488472621e-01, -2.9105521269417195e-01, -2.8618423155949174e-01, -2.7657873213483791e-01, -2.6044534346934506e-01, -2.3543750979707212e-01, -1.9871990208235915e-01, -1.4792362304676313e-01, -8.3086755914948171e-02];
a0 = [1.1827594641422032e+00, 2.6028043835324182e+00, 2.6106869684818865e+00, 2.6190302872996680e+00, 2.5560458896093023e+00, 2.5725846740321656e+00, 2.5903646624985202e+00, 2.6093659771512416e+00, 2.6293564841039281e+00, 2.6501416083459737e+00, 2.6713600872547070e+00, 2.6927013992181057e+00, 2.7137370738418181e+00, 2.7337408708989686e+00, 2.7521446367167575e+00, 2.7682282962016354e+00, 2.7813301841139588e+00, 2.7909589081873589e+00, 2.7965317145635229e+00, -1.6426539238850163e-01, -1.6700627926657657e-01, -1.6962948564648189e-01, -2.5655820114821737e-01, -2.6281823492276885e-01, -2.6878943036158198e-01, -2.7412621510560264e-01, -2.7859283994503875e-01, -2.8171724561258316e-01, -2.8314356544942187e-01, -2.8204362255678495e-01, -2.7763860555611097e-01, -2.6859468523822499e-01, -2.5315874222858242e-01, -2.2901648307616207e-01, -1.9339626338720795e-01, -1.4398734515438805e-01, -8.0870313921172193e-02];

% Applies all manoeuvres (using the nuclear engine)
[os, ths, dvs] = applyMultipleManoeuvres(o0, ttm_offs(1, 1:(Nm1 + Nm2)) + 382 * a0(1, 1), ...
            dvm(1, 1:(Nm1 + Nm2)), a0(1, 2:(Nm1 + Nm2 + 1)), ...
            a0(1, (Nm1 + Nm2 + 2):(2 * (Nm1 + Nm2) + 1)), mu_K);
ths(2, (Nm1 + Nm2 + 1)) = timeToArg(ttm_offs(1, (Nm1 + Nm2 + 1)) + 382 * a0(1, 1), ...
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

    % Last delta-v
    [~, vvf1] = paramsToPosVel(ths(2, (Nm1 + Nm2 + 1)), os(1, (Nm1 + Nm2 + 1)).p, ...
                      os(1, (Nm1 + Nm2 + 1)).e, os(1, (Nm1 + Nm2 + 1)).lan, ...
                      os(1, (Nm1 + Nm2 + 1)).incl, os(1, (Nm1 + Nm2 + 1)).aop, mu_K);
    [lan, incl, aop_th] = rotAngles(x_tps(:, end), [-x_tps(2, end); x_tps(1, end); 0]);
    [~, vvf2] = paramsToPosVel(0, norm(x_tps(:, end), 2), 0, lan, incl, aop_th, mu_K);
    [T] = refFrameLocal(x_tps(:, end), vvf1);
    ldv = T' * (vvf2 - vvf1);


%% - Numerical output -----------------------------------------------------

for idx_1 = 1:(Nm1 + Nm2)
    fprintf('dv1: %23.16f\ndv2: %23.16f\ndv3: %23.16f\nat : %23.16f\n\n', ...
        dvs(1, idx_1), dvs(2, idx_1), dvs(3, idx_1), ...
        (ttm_offs(1, idx_1) + 382 * a0(1, 1)))
end
fprintf('dv1: %23.16f\ndv2: %23.16f\ndv3: %23.16f\nat : %23.16f\n\n', ...
        ldv(1, 1), ldv(2, 1), ldv(3, 1), ...
        (ttm_offs(1, (Nm1 + Nm2 + 1)) + 382 * a0(1, 1))) 