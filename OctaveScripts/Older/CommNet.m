%%  ====  Header  =====================================================

fclose all;
close all
format long g
more off
clear
clc

% Number of points in a curve
N = 5001;

%%  ====  Planetary data  =============================================

% Initializes planetary index
ip = 0;

% Moho
ip = ip + 1;
pName{ip, 1} = 'Moho';
mu(ip, 1)    = 1.6860938e11;    % Gravitational parameter
R(ip, 1)     = 250000;          % Radius
r_SOI(ip, 1) = 9646663.0;       % Sphere of influence radius

% Eve
ip = ip + 1;
pName{ip, 1} = 'Eve';
mu(ip, 1)    = 8.1717302e12;    % Gravitational parameter
R(ip, 1)     = 700000;          % Radius
r_SOI(ip, 1) = 85109365;        % Sphere of influence radius

% Kerbin
ip = ip + 1;
pName{ip, 1} = 'Kerbin';
mu(ip, 1)    = 3.5316e12;       % Gravitational parameter
R(ip, 1)     = 600000;          % Radius
r_SOI(ip, 1) = 84159286;        % Sphere of influence radius

% Duna
ip = ip + 1;
pName{ip, 1} = 'Duna';
mu(ip, 1)    = 3.0136321e11;    % Gravitational parameter
R(ip, 1)     = 320000;          % Radius
r_SOI(ip, 1) = 47921949;        % Sphere of influence radius

% Dres
ip = ip + 1;
pName{ip, 1} = 'Dres';
mu(ip, 1)    = 2.1484489e10;    % Gravitational parameter
R(ip, 1)     = 138000;          % Radius
r_SOI(ip, 1) = 32832840;        % Sphere of influence radius

% Jool
ip = ip + 1;
pName{ip, 1} = 'Jool';
mu(ip, 1)    = 2.8252800e14;    % Gravitational parameter
R(ip, 1)     = 6000000;          % Radius
r_SOI(ip, 1) = 2455985200;        % Sphere of influence radius

% Eeloo
ip = ip + 1;
pName{ip, 1} = 'Eeloo';
mu(ip, 1)    = 7.4410815e10;    % Gravitational parameter
R(ip, 1)     = 210000;          % Radius
r_SOI(ip, 1) = 1.1908294e8;     % Sphere of influence radius

%%  ====  Relay data  =================================================

% Antenna powers
pa_1 = 100e9;
pa_2 = 5e3;
pa_3 = 500e3;

% Antenna powers (Jool)
pa_4 = pa_1 * (5^0.75);
pa_5 = 2e9; %5e6;

% Antenna ranges
ra1 = ones(ip, 1) * sqrt(pa_1 * pa_2);
ra2 = ones(ip, 1) * sqrt(pa_1 * pa_3);

% Antenna range (Jool)
ra1(6, 1) = sqrt(pa_4 * pa_2);
ra2(6, 1) = sqrt(pa_4 * pa_5);

% Relay altitude
hr       = (2.2e7) * ones(ip, 1); %ra1 - 4e5;
hr(1, 1) = 0.95 * r_SOI(1, 1);
hr(2, 1) = 0.5 * hr(3, 1);

hr(6, 1) = 1.2e9;

% Action oblate spheroid major axes
fsp = @(ra, hr) 0.5 * hr + sqrt((ra.^ 2) - 0.75 * (hr.^ 2));

% Action oblate spheroid minor axis
fma = @(ra, hr) sqrt((ra.^ 2) - (hr.^ 2));

% Relay positions
thr = [0; 2; 4] * (pi / 3);

for ii = 1:ip
    for jj = 1:3
        xr(ii, jj, :) = hr(ii, 1) * [cos(thr(jj, 1)), sin(thr(jj, 1))];
    end
end

%%  ====  Natural satellite clearance

% Satellite index
is = 0;

% Gilly
is = is + 1;
sName{is, 1}  = 'Gilly';
planet(is, 1) = 2;
ap(is, 1)     = 48825000;
pe(is, 1)     = 14175000;
r_SOIs(is, 1) = 126123.27;	% Sphere of influence radius

% Mun
is = is + 1;
sName{is, 1}  = 'Mun';
planet(is, 1) = 3;
ap(is, 1)     = 12000000;
pe(is, 1)     = 12000000;
r_SOIs(is, 1) = 2429559.1;	% Sphere of influence radius

% Minmus
is = is + 1;
sName{is, 1}  = 'Minmus';
planet(is, 1) = 3;
ap(is, 1)     = 47000000;
pe(is, 1)     = 47000000;
r_SOIs(is, 1) = 2247428.4;	% Sphere of influence radius

% Ike
is = is + 1;
sName{is, 1}  = 'Ike';
planet(is, 1) = 4;
ap(is, 1)     = 3296000;
pe(is, 1)     = 3104000;
r_SOIs(is, 1) = 11049598.9;	% Sphere of influence radius

GP = 388587;
heq_G = (((GP / (2 * pi)) ^ 2) * mu(2, 1)) ^ (1 / 3);

%%  ====  Circle coordinates  =========================================

xx = @(R, x0) x0(1, 1) + R * cos(linspace(-pi, pi, N));
yy = @(R, x0) x0(1, 2) + R * sin(linspace(-pi, pi, N));

%%  ====  Graphic output  =============================================

% For each planet aside from Moho
for ii=6:6 %ii = 1:ip
    hl(ii, 1) = figure(ii);
    title(pName{ii, 1})
    hold on

    % Planet
    plot(xx(R(ii, 1), [0, 0]), yy(R(ii, 1), [0, 0]), 'b')

    % Sphere of influence
    plot(xx(r_SOI(ii, 1), [0, 0]), yy(r_SOI(ii, 1), [0, 0]), 'c')

    % Relay positions and action spheres
    for jj = 1:3
        plot(xr(ii, jj, 1), xr(ii, jj, 2), 'k*')
        plot(xx(ra1(ii, 1), xr(ii, jj, :)), yy(ra1(ii, 1), xr(ii, jj, :)), 'g')
        plot(xx(ra2(ii, 1), xr(ii, jj, :)), yy(ra2(ii, 1), xr(ii, jj, :)), 'r')
    end

    % Collective action spheres
    plot(xx(fsp(ra1(ii, 1), hr(ii, 1)), [0, 0]), yy(fsp(ra1(ii, 1), hr(ii, 1)), [0, 0]), 'y')
    plot(xx(fsp(ra2(ii, 1), hr(ii, 1)), [0, 0]), yy(fsp(ra2(ii, 1), hr(ii, 1)), [0, 0]), 'm')

    axis equal
    hold off

    fprintf('For %s:\n', pName{ii, 1})
    fprintf('Ellipsoid 1 height vs planet radius       : %g / %g = %g\n', fma(ra1(ii, 1), hr(ii, 1)), R(ii, 1), (fma(ra1(ii, 1), hr(ii, 1)) / R(ii, 1)))
    fprintf('Ellipsoid 2 height vs Sphere of influence : %g / %g = %g\n', fma(ra2(ii, 1), hr(ii, 1)), r_SOI(ii, 1), (fma(ra2(ii, 1), hr(ii, 1)) / r_SOI(ii, 1)))
    fprintf('T = %g [s]\n\n', 2 * pi / sqrt(mu(ii, 1) / ((hr(ii, 1)) ^ 3)))

end

for ii = 1:is

    in = [(pe(ii, 1) - r_SOIs(ii, 1)), (ap(ii, 1) + r_SOIs(ii, 1))];
    fprintf('%s\n', sName{ii, 1})
    fprintf('[%g - %g - %g]\n\n', in(1, 1), hr(planet(ii, 1), 1), in(1, 2))

end