function [t, d, th1, th2] = entryInSphere(ti, rS, o1, o2, mu, of, tta)
%
%   function [t, d, th1, th2] = entryInSphere(ti, rS, o1, o2, mu, of, tta)
%
%   Given a time interval, finds the entry time of object 01 in a sphere
%   of radius rS around object o2
%

    % If no additional arguments for "optiFun" were specified
    if nargin < 6
        of = [];
    end

    % Maximum number of iterations, if unspecified
    if ~isfield(of, 'nIterMax') || numel(of.nIterMax)
        of.nIterMax = 100;
    end

    % First derivative tolerance, if unspecified
    if ~isfield(of, 'd1Tol') || numel(of.d1Tol)
        of.d1Tol = eps(mean(ti)) * (2 ^ 26);
    end

    % Number of section points, if unspecified
    if ~isfield(of, 'nPoints') || numel(of.nPoints)
        of.nPoints = 12;
    end

    % Variable increment, if unspecified
    if ~isfield(of, 'dx') || numel(of.dx)
        of.dx = eps(mean(ti)) * (2 ^ 26);
    end

    % If no additional arguments for "structTimeToPosVelArg" were specified
    if nargin < 7
        tta = [];
    end

    % Minimum absolute distance
    [tm, dm] = closestPass(ti, o1, o2, mu, of, tta);

    % Time interval
    tn = [ti(1, 1), tm];

    % Cost function
    cf = @(t) (sqrt(sum(((structTimeToPosVelArg(t, o1, mu, tta) - ...
        structTimeToPosVelArg(t, o2, mu, tta)).^ 2), 1)) - rS).^ 2 + ...
        (t > tm).* ((dm * (t - tm)).^ 2);

    % Finds the instant at which distance is minimum
    [t, ~, ~] = optiFun(tn, cf, of.dx, of.nPoints, of.d1Tol, of.nIterMax);

    %tt = (t-10):of.dx:(t+10);
    %{
    tt = linspace(ti(1, 1), ti(1, 2), 5001);
    cc = zeros(size(tt));
    for idx_1 = 1:size(tt, 2)
        cc(1, idx_1) = cf(tt(1, idx_1));
    end
    plot(tt, cc, 'b')
    hold on
    plot(t, cf(t), 'ro')
    plot(tm, cf(tm), 'ms')
    hold off
    %}

    % If the actual distance was requested
    if nargout >= 2
        d = sqrt(cf(t));
    end

    % If the argument of the first object was requsted
    if nargout >= 3
        th1 = timeStructToArg(t, o1, mu, tta);
    end

    % If the argument of the second object was requested
    if nargout >= 4
        th2 = timeStructToArg(t, o2, mu, tta);
    end

end