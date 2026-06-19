%% --- Header -------------------------------------------------------------

fclose all;
close all
format long g
more off
clear
clc

vstr = version;
vnum = str2double(vstr(1, (end - 5):(end - 2)));

N = 5001;
n = 51;

%% --- Initial data -------------------------------------------------------

% Duna data
R_D     =	320000;
mu_D    =   3.0136321e11;
rSOI_D  =   47921949;

% Ike data
R_I     =   130000;
mu_I    =   1.8568369e10;
rSOI_I  =   1049598.9;

% Ike orbit
[oI]	=	rawOrbitToStruct([3104000, 3296000, 0.03, 0, 0.2, 0]);
oI.peT  =   argStructToTime(-pi, oI,mu_D) * (1.7 / pi);

% Duna relays distance
dRD     =   22000000;

%% --- Ike's orbit --------------------------------------------------------

th = linspace(-pi, pi, N);
[xI, ~, ~] = structArgToPosVelTime(th, oI, mu_D);

%% --- Ike's torus of doom ------------------------------------------------

nIT = ceil(n * 0.5);

th_IT = linspace(-pi, pi, nIT);
[xx_IT, ~, ~] = structArgToPosVelTime(th_IT, oI, mu_D);

ii_r = xx_IT./ repmat(sqrt(sum((xx_IT.^ 2), 1)), 3, 1);
ii_k = repmat([0; 0; 1], 1, nIT);
ii_k = ii_k - ii_r.* repmat(sum((ii_k.* ii_r), 1), 3, 1);
ii_k = ii_k./ repmat(sqrt(sum((ii_k.^ 2), 1)), 3, 1);

XX_IT = repmat(reshape(xx_IT', nIT, 1, 3), 1, nIT, 1);
II_r  = repmat(reshape(ii_r', nIT, 1, 3), 1, nIT, 1);
II_k  = repmat(reshape(ii_k', nIT, 1, 3), 1, nIT, 1);
TH_IT = repmat(th_IT, nIT, 1, 3);
XX_T  = XX_IT + rSOI_I * (II_r.* cos(TH_IT) + II_k.* sin(TH_IT));

%% --- Spheres ------------------------------------------------------------

% Duna
[XS_D, YS_D, ZS_D] = genSphere(R_D, n);
CS_D(:, :, 1) =  ones(size(XS_D));
CS_D(:, :, 2) = zeros(size(XS_D));
CS_D(:, :, 3) = zeros(size(XS_D));

% Duna's sphere of influence
[XS_rSOI_D, YS_rSOI_D, ZS_rSOI_D] = genSphere(rSOI_D, ceil(n * 0.5));


%% --- Duna relays --------------------------------------------------------

% Maximum relay distance
dMR     = sqrt((dRD ^ 2) - (R_D ^ 2));
dMR2    = dMR ^ 2;

% Relay orbits
o_DR1 = rawOrbitToStruct([dRD, dRD, 0, 0, 0,0]);%0.875, 0]);
o_DR2 = rawOrbitToStruct([dRD, dRD, 0, 120, 0,0]);%0.875, 0])
o_DR3 = rawOrbitToStruct([dRD, dRD, 0, 240, 0,0]);%0.875, 0])

[xx_DR1, ~, ~] = structArgToPosVelTime(th, o_DR1, mu_D);
[xx_DR2, ~, ~] = structArgToPosVelTime(th, o_DR2, mu_D);
[xx_DR3, ~, ~] = structArgToPosVelTime(th, o_DR3, mu_D);

fprintf('Loading sphere mesh...\n')
load('sphere.mat');
%load('~/Documents/sphereHuge.mat')
fprintf('Sphere mesh loaded.\n')
Sphere = R_D * Sphere(2:end, :);
%Sphere = R_D * Sphere([2:7, 9:end], :);

NS	= size(Sphere, 1);
B	= true(NS, 1);
Q   = false(NS, 1);
pt  = zeros(N, 1);

% For each point
fprintf('\n\nCoverage checking... 000.00')
for idx_1 = 1:N

    fprintf('\b\b\b\b\b\b%5.1f%c', 100 * (idx_1 / N), 37)

    q = ((sum(((Sphere - repmat((xx_DR1(:, idx_1)'), NS, 1)).^ 2), 2) <= (dMR2)) | ...
         (sum(((Sphere - repmat((xx_DR2(:, idx_1)'), NS, 1)).^ 2), 2) <= (dMR2)) | ...
         (sum(((Sphere - repmat((xx_DR3(:, idx_1)'), NS, 1)).^ 2), 2) <= (dMR2)));

     pt(idx_1, 1) = 100 * (sum(q) / NS);

     B = B & q;
     Q = Q | q;

end
fprintf('\n\n\n%10.6f%c of the surface is always covered.\n', 100 * (sum(B) / NS), 37)
fprintf('\n%10.6f%c of the surface is covered at some time.\n', 100 * (sum(Q) / NS), 37)
fprintf('\n%10.6f%c of the surface is covered on average.\n', mean(pt), 37)
fprintf('\nMaximum surface covered: %10.6f%c.\n', max(pt), 37)
fprintf('\nMinimum surface covered: %10.6f%c.\n', min(pt), 37)

%% --- Graphic output -----------------------------------------------------

figure(1)
if vnum >= 2019
	addToolbarExplorationButtons(gcf)
end
hold on
axis equal
plot3(0,0,0,'ro')
plot3(0,0,0,'r*')
h_D = surf(XS_D, YS_D, ZS_D, CS_D);
h_D.EdgeColor = 'none';

plot3(xI(1, :), xI(2, :), xI(3, :), 'b', 'LineWidth', 2)

plot3(xx_DR1(1, :), xx_DR1(2, :), xx_DR1(3, :), 'color', [0.25, 1, 0.25], 'LineWidth', 2)
plot3(xx_DR2(1, :), xx_DR2(2, :), xx_DR2(3, :), 'color', [0, 1, 0], 'LineWidth', 2)
plot3(xx_DR3(1, :), xx_DR3(2, :), xx_DR3(3, :), 'color', [0, 0.75, 0], 'LineWidth', 2)

plot3(Sphere(~B, 1), Sphere(~B, 2), Sphere(~B, 3), 'bo')
plot3(Sphere(~Q, 1), Sphere(~Q, 2), Sphere(~Q, 3), 'y*', 'MarkerSize', 2)

h_rSOI_I = surf(reshape(XX_T(:, :, 1), nIT, nIT), ...
                reshape(XX_T(:, :, 2), nIT, nIT), ...
                reshape(XX_T(:, :, 3), nIT, nIT));
h_rSOI_I.FaceColor = 'none';
h_rSOI_I.EdgeColor = [0.625, 0.125, 0.125];

%h_rSOI_D = surf(XS_rSOI_D, YS_rSOI_D, ZS_rSOI_D);
%h_rSOI_D.FaceColor = 'none';
%h_rSOI_D.EdgeColor = [0.75, 0.75, 0.75];
hold off

figure(2)
if vnum >= 2019
	addToolbarExplorationButtons(gcf)
end
hold on
plot(th, pt, 'b')
hold off
