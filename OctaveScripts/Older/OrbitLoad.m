%%

fclose all;
close all
format long g
more off
clear
clc

%%

N = 72;

cellForm = csvToCells('TwoOrbits3.csv');

optStep = false;

% Minmus data
mu_Ms  = str2double(cellForm{14, 2});
R_Ms   = str2double(cellForm{15, 2});

% Initial orbit data
ap_1   = str2double(cellForm{8, 2}) + R_Ms;
pe_1   = str2double(cellForm{9, 2}) + R_Ms;
pe_1_T = str2double(cellForm{10, 2});
incl_1 = zeroToMinusPi(degToGrads(str2double(cellForm{11, 2})));
lan_1  = zeroToMinusPi(degToGrads(str2double(cellForm{12, 2})));
aop_1  = zeroToMinusPi(degToGrads(str2double(cellForm{13, 2})));

% Final orbit data
ap_2   = str2double(cellForm{2, 2}) + R_Ms;
pe_2   = str2double(cellForm{3, 2}) + R_Ms;
incl_2 = zeroToMinusPi(degToGrads(str2double(cellForm{4, 2})));
lan_2  = zeroToMinusPi(degToGrads(str2double(cellForm{5, 2})));
aop_2  = zeroToMinusPi(degToGrads(str2double(cellForm{6, 2})));

% Minmus data
%mu_Ms  = str2double(cellForm{14, 2});
%R_Ms   = str2double(cellForm{15, 2});

% Forbidden times
FTimes = [];
for idx_1 = 17:size(cellForm, 1)
  FTimes = [FTimes; str2double(cellForm{idx_1, 1})];
end

%%

% Sphere of influence radius
a    = 4.7e7;
mu_K = 3.5316e12;
R_soi_Ms = a * ((mu_Ms / mu_K) ^ (2 / 5));

% Orbital computations
p1 = 2 * ap_1 * pe_1 / (ap_1 + pe_1);
e1 = (ap_1 - pe_1) / (ap_1 + pe_1);

p2 = 2 * ap_2 * pe_2 / (ap_2 + pe_2);
e2 = (ap_2 - pe_2) / (ap_2 + pe_2);

R1 = rOrb(lan_1, incl_1, aop_1);
R2 = rOrb(lan_2, incl_2, aop_2);

if e1 >= 1
    th1_0 = acos((1 / e1) * ((p1 / R_soi_Ms) - 1));
else
    th1_0 = pi;
end
th1 = linspace(-th1_0, th1_0, N);
if e2 >= 1
    th2_0 = acos((1 / e2) * ((p2 / R_soi_Ms) - 1));
else
    th2_0 = pi;
end
th2 = linspace(-th2_0, th2_0, N);

t_1 = argToTime(th1, p1, e1, mu_Ms, pe_1_T);

xx1 = repmat((p1./ (1 + e1 * cos(th1))), 3, 1).* [cos(th1); sin(th1); zeros(1, N)];
xx2 = repmat((p2./ (1 + e2 * cos(th2))), 3, 1).* [cos(th2); sin(th2); zeros(1, N)];

xx1_e = (p1 / (1 + e1 * cos(-th1_0))) * [cos(-th1_0); -sin(th1_0); 0];
xx1_e = [xx1_e, ((p1 / (1 + e1 * cos(th1_0))) * [cos(th1_0); sin(th1_0); 0])];

xx1 = R1 * xx1;
xx2 = R2 * xx2;
xx1_e = R1 * xx1_e;

%% Wrong time intervals

dmn = 600;
WTints = [(FTimes - dmn), (FTimes + dmn)];

t_1v = t_1;

for idx_1 = 1:size(WTints, 1)
   t_1v(t_1v >= WTints(idx_1, 1) & t_1v <= WTints(idx_1, 2)) = NaN;
end

t_1v = t_1v(~isnan(t_1v));
th1v = timeToArg(t_1v, p1, e1, mu_Ms, pe_1_T);

xx1v = repmat((p1./ (1 + e1 * cos(th1v))), 3, 1).* [cos(th1v); sin(th1v); zeros(size(th1v))];
xx1v = R1 * xx1v;

% MSph
[Msph_X, Msph_Y, Msph_Z] = genSphere(R_Ms, 51);

th1_1 = 1.28318573174794;
th2_2 = 0.0872664625997164;



%[th1t, th2t, otp, dv1, dv2, er] = getTransferOrbitStatic([th1_1, th2_2, pi], paramsToVector(p1, e1, lan_1, incl_1, aop_1), paramsToVector(p2, e2, lan_2, incl_2, aop_2), mu_Ms);
%[cf, th1t, th2t, otp, dv1, dv2, er] = transferCostStaticTimes([th1_1, th2_2, pi], paramsToVector(p1, e1, lan_1, incl_1, aop_1), paramsToVector(p2, e2, lan_2, incl_2, aop_2), pe_1_T, FTimes, mu_Ms);
tic
%[cf, th1t, th2t, th0_3, otp, dv1, dv2, er] = transferCostStatOpt2p(th1_1, th2_2, paramsToVector(p1, e1, lan_1, incl_1, aop_1), paramsToVector(p2, e2, lan_2, incl_2, aop_2), mu_Ms);
[cf, th1t, th2t, th0_3, otp, dv1, dv2, er] = transferCostStaticTimes(th1_1, th2_2, paramsToVector(p1, e1, lan_1, incl_1, aop_1), paramsToVector(p2, e2, lan_2, incl_2, aop_2), pe_1_T, FTimes, mu_Ms);
toc
[p3, e3, lan_3, incl_3, aop_3] = vectorToParams(otp);

[x1s, v1s] = paramsToPosVel(th1_1, p1, e1, lan_1, incl_1, aop_1, mu_Ms);
[x2f, v2f] = paramsToPosVel(th2_2, p2, e2, lan_2, incl_2, aop_2, mu_Ms);

t1  = argToTime(th1_1, p1, e1, mu_Ms, pe_1_T);
while th1t > th2t
    th1t = th1t - 2 * pi;
end
tt1 = argToTime(th1t, p3, e3, mu_Ms, 0);
tt2 = argToTime(th2t, p3, e3, mu_Ms, 0);
t2  = t1 + tt2 - tt1;

if th1t > th2t
    th1t = th1t - 2 * pi;
end
th3 = linspace(th1t, th2t, N);

[xx3, v3] = paramsToPosVel(th3, p3, e3, lan_3, incl_3, aop_3, mu_Ms);

%{
ir1_t  = x1s / norm(x1s, 2);
k1_t   = cross(ir1_t, v1s);
k1_t   = k1_t / norm(k1_t, 2);
ith1_t = cross(k1_t, ir1_t);

ir2_t  = xx3(:, end) / norm(xx3(:, end), 2);
k2_t   = cross(ir2_t, v3(:, end));
k2_t   = k2_t / norm(k2_t, 2);
ith2_t = cross(k2_t, ir2_t);
%}

ith1_t = v1s / norm(v1s, 2);
k1_t   = cross(x1s, ith1_t);
k1_t   = k1_t / norm(k1_t, 2);
ir1_t = cross(ith1_t, k1_t);

ith2_t = v3(:, end) / norm(v3(:, end), 2);
k2_t   = cross(xx3(:, end), ith2_t);
k2_t   = k2_t / norm(k2_t, 2);
ir2_t = cross(ith2_t, k2_t);


round([t1, dot(dv1, ir1_t), dot(dv1, ith1_t), dot(dv1, k1_t)]' * 1000) / 1000
round([t2, dot(dv2, ir2_t), dot(dv2, ith2_t), dot(dv2, k2_t)]' * 1000) / 1000

plot3(xx2(1, :), xx2(2, :), xx2(3, :), 'b', 'LineWidth', 1)
hold on
plot3(xx3(1, :), xx3(2, :), xx3(3, :), 'k', 'LineWidth', 1)
plot3(xx1(1, :), xx1(2, :), xx1(3, :), 'r', 'LineWidth', 1)

plot3(xx3(1, 1), xx3(2, 1), xx3(3, 1), 'co', 'LineWidth', 1)
plot3(x1s(1, 1), x1s(2, 1), x1s(3, 1), 'cs', 'LineWidth', 1)

plot3(xx3(1, end), xx3(2, end), xx3(3, end), 'mo', 'LineWidth', 1)
plot3(x2f(1, 1), x2f(2, 1), x2f(3, 1), 'ms', 'LineWidth', 1)

%plot3(xx1v(1, :), xx1v(2, :), xx1v(3, :), 'gs')
%plot3(xx1_e(1, :), xx1_e(2, :), xx1_e(3, :), 'ro')
surf(Msph_X, Msph_Y, Msph_Z, zeros(size(Msph_X)));
axis equal
hold off

%%


t1 = argToTime(th1, p1, e1, mu_Ms, pe_1_T);
t1 = t1 - pe_1_T;
t2 = FTimes - pe_1_T;
t3 = WTints - pe_1_T;

c = zeros(size(t1));
g = @(x, a) exp(-((x / a).^ 2));

for idx_1 = 1:size(t2, 1)
    c = c + 1000 * g((t1 - t2(idx_1, 1)), 509.593080172811);
end

figure(2)
plot(t1, c, 'r')
hold on
plot(t2(t2 >= min(t1) & t2 <= max(t1)), zeros(sum(t2 >= min(t1) & t2 <= max(t1)), 1), 'gs')

plot(t3(t3(:, 1) >= min(t1) & t3(:, 1) <= max(t1), 1), zeros(sum(t3(:, 1) >= min(t1) & t3(:, 1) <= max(t1), 1), 1), 'bo')
plot(t3(t3(:, 2) >= min(t1) & t3(:, 2) <= max(t1), 2), zeros(sum(t3(:, 2) >= min(t1) & t3(:, 2) <= max(t1), 1), 1), 'bo')
hold off

th1_1s = -th1_0;
th2_2s = -pi;
th1_1f = th1_0; %th1_1s + 71 * (2 * pi) / 72;
th2_2f = th2_2s + 71 * (2 * pi) / 72;
mnv    = Inf;
ff1 = 0;
ff2 = 0;

if optStep

    for idx_1 = 1:5

        fprintf('Optimization level %g of %g, current optimum: %g\n', idx_1, 5, mnv)

        dth1 = (th1_1f - th1_1s) / 71;
        dth2 = (th2_2f - th2_2s) / 71;

        [th1_1M, th2_2M] = meshgrid(linspace(th1_1s, th1_1f, 72), linspace(th2_2s, th2_2f, 72));

        for idx_2 = 1:72

            fprintf('\t\tStep %g of %g, current optimum: %g\n', idx_2 * 72, 72 * 72, mnv)

            for idx_3 = 1:72

                [cf] = transferCostStaticTimes(th1_1M(idx_2, idx_3), th2_2M(idx_2, idx_3), paramsToVector(p1, e1, lan_1, incl_1, aop_1), paramsToVector(p2, e2, lan_2, incl_2, aop_2), pe_1_T, FTimes, mu_Ms);

                if cf < mnv
                    mnv = cf;
                    ff1 = idx_2;
                    ff2 = idx_3;
                end

                th1_1s = th1_1M(idx_2, idx_3) - dth1;
                th1_1f = th1_1M(idx_2, idx_3) + dth1;
                th2_2s = th2_2M(idx_2, idx_3) - dth2;
                th2_2f = th2_2M(idx_2, idx_3) + dth2;

            end
        end

    end

    [cf, th1t, th2t, th0_3, otp, dv1, dv2, er] = transferCostStaticTimes(th1_1M(ff1, ff2), th2_2M(ff1, ff2), paramsToVector(p1, e1, lan_1, incl_1, aop_1), paramsToVector(p2, e2, lan_2, incl_2, aop_2), pe_1_T, FTimes, mu_Ms);

    %tic
    %th0_3 = linspace(-pi, pi, 3601);
    %cst = zeros(1, 3600);
    %for idx_1 = 1:3600
    %[cst(1, idx_1), th1t, th2t, otp, dv1, dv2, er] = transferCostStaticTimes([th1_1, th2_2, th0_3(1, idx_1)], paramsToVector(p1, e1, lan_1, incl_1, aop_1), paramsToVector(p2, e2, lan_2, incl_2, aop_2), pe_1_T, FTimes, mu_Ms);
    %end
    %toc

end


%{
dvr = -(185:2.5:195)';
aps = [45.4913;10.7465;5.55822;3.47593;2.35945] * 1000000 + R_Ms;
pes = ones(5, 1) * 1879380 + R_Ms;
[ps, es] = apoPeriToPe(pes, aps);
[x1, v1] = paramsToPosVel(0, p1, e1, lan_1, incl_1, aop_1, mu_Ms);
ir = x1 / norm(x1);
k = cross(ir, v1);
k = k / norm(k);
ith = cross(k, ir);
xc = zeros(5, 3);
vc = zeros(5, 3);
for idx_1 = 1:5
[xct, vct] = paramsToPosVel(0, ps(idx_1, 1), es(idx_1, 1), lan_1, incl_1, aop_1, mu_Ms);
xc(idx_1, :) = xct';
vc(idx_1, :) = vct';
end
dvs = vc - repmat(v1', 5, 1);
dvf = dvs * ith;
figure(3)
plot(dvr, dvf, 'b')
%}