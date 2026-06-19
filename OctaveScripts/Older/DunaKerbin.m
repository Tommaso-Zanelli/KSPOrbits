%% --- Header -------------------------------------------------------------

fclose all;
close all
format long g
more off
clear
clc

vstr = version;
vnum = str2double(vstr(1, (end - 5):(end - 2)));

spheresActive = false;

%% --- Sun data -----------------------------------------------------------

mu_S   = 1.1723328e18;      % Standard gravitational parameter
R_S    = 261600000;         % Radius

%% --- Kerbin orbit / data ------------------------------------------------

mu_K   = 3.5316000e12;      % Standard gravitational parameter
R_K    = 600000;            % Radius
rSOI_K = 84159286;          % Sphere of influence radius
ime_K  = 3.14;              % Initial mean anomaly

% Orbit
[oK] = rawOrbitToStructOrbit([[13599840256, 13599840256, 0 0, 0, 0]; ...
            zeros(2, 6)]);

% Orbital period
TK = 2 * argStructToTime(pi, oK, mu_S);

% Periapsis time
oK.peT = -TK * (ime_K / (2 * pi));

%% --- Duna orbit / data --------------------------------------------------

mu_D   = 3.0136321e11;      % Standard gravitational parameter
R_D    = 320000;            % Radius
rSOI_D = 47921949;          % Sphere of influence radius
ime_D  = 3.14;              % Initial mean anomaly

% Orbit
[oD] = rawOrbitToStructOrbit([[19669121365, 21783189163, 0.051 135.5, ...
            0.06, 0]; zeros(2, 6)]);

% Orbital period
TD = 2 * argStructToTime(pi, oD, mu_S);

% Periapsis time
oD.peT = -TD * (ime_D / (2 * pi));

%% --- Other data ---------------------------------------------------------

5085789


ct = 3054600;
ti = 3 * 6 * 3600;
tf = 84 * 6 * 3600;
%(5085789-ct-ti-tf)/3600

ri = 6 * 6 * 3600;
rf = 42 * 6 * 3600;

N = 5001;

%% --- Transfer orbit -----------------------------------------------------

th0 = linspace(-pi, pi, N);
%th0 = -0.5510298098286258;
%th0 = -0.5734760746443260;
th0 = -0.6515450820811063;

if size(th0, 2) == 1
    transfer_exists = true;
    mdt = 0;
else
    dt = zeros(1, N);
    transfer_exists = false;
    mdt = Inf;
end

    n_iter = 0;

while (mdt > 1e-16 && n_iter < 100 && ~transfer_exists) || (n_iter < 1 && transfer_exists)

    [xxi_K, v1, ~] = structTimeToPosVelArg(ct + ti, oK, mu_S);
    [xxf_D, v2, ~] = structTimeToPosVelArg(ct + ti + tf, oD, mu_S);
    [oT.lan, oT.incl, aop_th] = rotAngles(xxi_K, (xxf_D-xxi_K));

    for idx_1 = 1:size(th0, 2)
        d_th = atan2(norm(cross(xxi_K, xxf_D), 2), dot(xxi_K, xxf_D));
        rho_1 = norm(xxi_K, 2);
        rho_2 = norm(xxf_D, 2);
        %[th_1, th_2, oT.p, oT.e, oT.lan, oT.incl, oT.aop] = transferOrbit(xxi_K, xxf_D, th0(1, idx_1), 0)
        oT.p = (rho_1*(rho_2*cos(d_th + th0(1, idx_1)) - rho_2*cos(th0(1, idx_1))))/(rho_2*cos(d_th + th0(1, idx_1)) - rho_1*cos(th0(1, idx_1)));
        oT.e = (rho_1 - rho_2)/(rho_2*cos(d_th + th0(1, idx_1)) - rho_1*cos(th0(1, idx_1)));
        th_1 = th0(1, idx_1);
        th_2 = th0(1, idx_1) + d_th;
        oT.aop = aop_th - th_1;
        oT.peT = 0;
        dt(1, idx_1) = argStructToTime(th_2, oT, mu_S) - argStructToTime(th_1, oT, mu_S);
        oT.peT = ct + ti - argStructToTime(th_1, oT, mu_S);
    end

    if ~transfer_exists
        dt(isnan(dt)) = Inf;
        dt = abs((dt - tf).^ 2);

        [mdt, imdt] = min(dt);
        newl = [th0(1, max([1, (imdt - 1)])); th0(1, min([N, (imdt + 1)]))];
        fprintf('Minimum (seconds): %g\n', sqrt(mdt))
        fprintf('New interval: %24.16f - %24.16f\n', newl(1, 1), newl(2, 1))
        th0 = linspace(newl(1, 1), newl(2, 1), N);
    end

    n_iter = n_iter + 1;

end

%% --- Return transfer ----------------------------------------------------

th0r = linspace(-pi, pi, N);
%th0r = 1.4539773959531785;
%th0r = 1.8233331333118599;
th0r = -2.9607660442793375;

if size(th0r, 2) == 1
    return_exists = true;
    mdt = 0;
else
    dt_r = zeros(1, N);
    return_exists = false;
    mdt = Inf;
end

    n_iter = 0;

while (mdt > 1e-16 && n_iter < 100 && ~return_exists) || (n_iter < 1 && return_exists)

    [xxri_D, v3, ~] = structTimeToPosVelArg(ct + ti + tf + ri, oD, mu_S);
    [xxrf_K, v4, ~] = structTimeToPosVelArg(ct + ti + tf + ri + rf, oK, mu_S);
    [oR.lan, oR.incl, aop_th_R] = rotAngles(xxri_D, (xxrf_K-xxri_D));

    for idx_1 = 1:size(th0r, 2)
        d_thr = atan2(norm(cross(xxri_D, xxrf_K), 2), dot(xxri_D, xxrf_K));
        rho_1r = norm(xxri_D, 2);
        rho_2r = norm(xxrf_K, 2);
        %[th_1, th_2, oT.p, oT.e, oT.lan, oT.incl, oT.aop] = transferOrbit(xxi_K, xxf_D, th0(1, idx_1), 0)
        oR.p = (rho_1r*(rho_2r*cos(d_thr + th0r(1, idx_1)) - rho_2r*cos(th0r(1, idx_1))))/(rho_2r*cos(d_thr + th0r(1, idx_1)) - rho_1r*cos(th0r(1, idx_1)));
        oR.e = (rho_1r - rho_2r)/(rho_2r*cos(d_thr + th0r(1, idx_1)) - rho_1r*cos(th0r(1, idx_1)));
        th_1r = th0r(1, idx_1);
        th_2r = th0r(1, idx_1) + d_thr;
        oR.aop = aop_th_R - th_1r;
        oR.peT = 0;
        dt_r(1, idx_1) = argStructToTime(th_2r, oR, mu_S) - argStructToTime(th_1r, oR, mu_S);
        oR.peT = ct + ti + + tf + ri - argStructToTime(th_1r, oR, mu_S);
    end

    if ~return_exists
        dt_r(isnan(dt_r)) = Inf;
        dt_r = abs((dt_r - rf).^ 2);

        [mdt, imdt] = min(dt_r);
        newl = [th0r(1, max([1, (imdt - 1)])); th0r(1, min([N, (imdt + 1)]))];
        fprintf('Minimum (seconds): %g\n', sqrt(mdt))
        fprintf('New interval: %24.16f - %24.16f\n', newl(1, 1), newl(2, 1))
        th0r = linspace(newl(1, 1), newl(2, 1), N);
    end

    n_iter = n_iter + 1;

end

%% --- Output data --------------------------------------------------------

% Current positions
xxc_S = [0; 0; 0];
[xxc_K, ~, ~] = structTimeToPosVelArg(ct, oK, mu_S);
[xxc_D, ~, ~] = structTimeToPosVelArg(ct, oD, mu_S);

if spheresActive

    % Spheres
    [XS_S, YS_S, ZS_S] = genSphere(R_S, 26);
    [XS_K, YS_K, ZS_K] = genSphere(R_K, 26);
    [XS_D, YS_D, ZS_D] = genSphere(R_D, 26);

    % Corrects spheres
    XS_K = XS_K + xxc_K(1, 1);
    YS_K = YS_K + xxc_K(2, 1);
    ZS_K = ZS_K + xxc_K(3, 1);
    XS_D = XS_D + xxc_D(1, 1);
    YS_D = YS_D + xxc_D(2, 1);
    ZS_D = ZS_D + xxc_D(3, 1);

end

% Orbits
th_o = linspace(-pi, pi, 5001);
[xxo_K, ~, ~] = structArgToPosVelTime(th_o, oK, mu_S);
[xxo_D, ~, ~] = structArgToPosVelTime(th_o, oD, mu_S);

if spheresActive

    % Sphere colours
    CS(:, :, 1) = ones(size(ZS_S));
    CS(:, :, 2) = ones(size(ZS_S));
    CS(:, :, 3) = zeros(size(ZS_S));
    CK(:, :, 1) = zeros(size(ZS_K));
    CK(:, :, 2) = zeros(size(ZS_K));
    CK(:, :, 3) = ones(size(ZS_K));
    CD(:, :, 1) = ones(size(ZS_D));
    CD(:, :, 2) = zeros(size(ZS_D));
    CD(:, :, 3) = zeros(size(ZS_D));

end

% Transfer
if transfer_exists
    th_t = linspace(th_1, th_2, 5001);
    [xxo_T, vvo_T, ~] = structArgToPosVelTime(th_t, oT, mu_S);

    v1t = vvo_T(:, 1);
    v2t = vvo_T(:, end);
    dv1 = v1t - v1;
    dv2 = v2 - v2t;
    fprintf('First  burn dv: %8.1f\n', norm(dv1, 2))
    fprintf('Second burn dv: %8.1f\n', norm(dv2, 2))
end

% Return
if return_exists
    th_t = linspace(th_1r, th_2r, 5001);
    [xxo_R, vvo_R, ~] = structArgToPosVelTime(th_t, oR, mu_S);

    v1r = vvo_R(:, 1);
    v2r = vvo_R(:, end);
    dv1r = v1r - v3;
    dv2r = v4 - v2r;
    fprintf('Third  burn dv: %8.1f\n', norm(dv1r, 2))
    fprintf('Fourth burn dv: %8.1f\n', norm(dv2r, 2))
end

% Delta v
if transfer_exists && return_exists

    dv_tot = norm(dv1, 2) + norm(dv2, 2) + norm(dv1r, 2) + norm(dv2r, 2);
    fprintf('\nTotal dv: %g\n', dv_tot)

end

%% --- Graphic output -----------------------------------------------------

figure(1)
m = 25000000000;
hold on
if vnum >= 2019
    addToolbarExplorationButtons(gcf)
end
plot3(xxo_K(1, :), xxo_K(2, :), xxo_K(3, :), 'color', [0, 0, 0.5])
plot3(xxo_D(1, :), xxo_D(2, :), xxo_D(3, :), 'color', [0.5, 0, 0])

plot3(xxc_K(1, :), xxc_K(2, :), xxc_K(3, :), 'o', 'color', [0, 0, 1])
plot3(xxc_D(1, :), xxc_D(2, :), xxc_D(3, :), 'o', 'color', [1, 0, 0])
plot3(xxc_S(1, :), xxc_S(2, :), xxc_S(3, :), 'o', 'color', [1, 1, 0])

if transfer_exists
    plot3(xxi_K(1, :), xxi_K(2, :), xxi_K(3, :), '*', 'color', [0, 0, 0])
    plot3(xxf_D(1, :), xxf_D(2, :), xxf_D(3, :), '*', 'color', [0.25, 0.25, 0.25])

    plot3(xxo_T(1, :), xxo_T(2, :), xxo_T(3, :), 'c')
end

if return_exists
    plot3(xxri_D(1, :), xxri_D(2, :), xxri_D(3, :), '*', 'color', [0.75, 0.75, 0.75])
    plot3(xxrf_K(1, :), xxrf_K(2, :), xxrf_K(3, :), '*', 'color', [0.375, 0.375, 0.375])

    plot3(xxo_R(1, :), xxo_R(2, :), xxo_R(3, :), 'm')
end

if spheresActive

    hS_S = surf(XS_S, YS_S, ZS_S, CS);
    hS_K = surf(XS_K, YS_K, ZS_K, CK);
    hS_D = surf(XS_D, YS_D, ZS_D, CD);
    hS_S.EdgeColor = 'none';
    hS_K.EdgeColor = 'none';
    hS_D.EdgeColor = 'none';

end

axis ([-m, m, -m, m, -m, m])
axis equal
hold off