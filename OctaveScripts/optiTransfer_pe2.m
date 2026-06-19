function [p, e, th0] = optiTransfer_pe2(r1, r2, dth, costf, delta_th, n_s, nm_nt, ep_nt)
%
%   function [p, e] = optiTransfer_pe2(r1, r2, dth, costf)
%
%   Given two co-planar points, their angular distance, and the velocities
%   required at them, computes the optimal transfer orbit between them.
%
%   costf must be in the form: c = costf(p, e, th0, dth)
%

    % If unspecified, sets the tolerance for the Newton's method iterations
    if nargin < 8
        ep_nt = 1e-9;
    end

    % If unspecified, sets th number of Newton's method iterations
    if nargin < 7
        nm_nt  = 1000;
    end

    % If unspecified, sets the number of section intervals
    if nargin < 6
        n_s = 12;
    end

    % If unspecified, sets the delta_th
    if nargin < 5
        delta_th = 4e-6;
    end

    % Degenerate cases
    if dth == 0

        % Same points
        if r1 == r2
            p = r1;
            e = 0;
            return;

        % Aligned points
        else

        end
    end

    % Radius ratio
    xi = r1 / r2;

    % Determinant phase
    phi_d = atan2(sin(dth), (cos(dth) - xi));

    % p numerator's phase
    phi_p = atan2(sin(dth), (cos(dth) - 1));

    % Finds search intervals
    [pp1, pn1, pp2, pn2] = posNegAngles(phi_d, phi_p);
    fprintf('%g %g %g %g\n', pp1, pp2, pn1, pn2)

    % Chooses the interval for which "e" is positive
    if xi > 1
        thi = [pp1 + (2 * delta_th), pp2 - (2 * delta_th)];
    else
        thi = [pn1 + (2 * delta_th), pn2 - (2 * delta_th)];
    end

    th_s = linspace(thi(1, 1), thi(1, 2), 501);
    figure(3)
    hold on
        plot(thi, [-0.3, -0.3], 's-', 'color', [0.5, 0, 0])
        plot(th_s, cos(th_s + phi_d), '--')
        plot(th_s, cos(th_s + phi_p), '--')
    hold off
    cc_s = zeros(size(th_s));
    dd_s = zeros(size(th_s));
    pp_s = zeros(size(th_s));
    for idx_1 = 1:size(th_s, 2)
        cc_s(1, idx_1) = cost_pe(r1, xi, dth, th_s(1, idx_1), costf);
        dd_s(1, idx_1) = det_pe(xi, dth, th_s(1, idx_1));
        pp_s(1, idx_1) = ppr_pe(r1, dth, th_s(1, idx_1), dd_s(1, idx_1));
    end
    figure(2)
    hold on
    plot(th_s, cc_s, 'r')
    plot(th_s, dd_s, 'k')
    plot(th_s, pp_s, 'b')
    hold off

    % Interal spanning portion
    if n_s > 2
        th_span = linspace(thi(1, 1), thi(1, 2), n_s);
        cs_span = zeros(1, n_s);
        for idx_1 = 1:n_s
            cs_span(1, idx_1) = cost_pe(r1, xi, dth, th_span(1, idx_1), costf);
        end
        [cc, im] = min(cs_span);
        th0 = th_span(1, im);
        figure(2)
        hold on
            plot(th_span, cs_span, 'mo')
            plot(th0, cc, 'k*')
        hold off
    % If two or less test points were requested
    else

        % Chooses the midpoint of the interval as the initial guess
        th0 = 0.5 * (thi(1, 1) + thi(1, 2));
        cc  = cost_pe(r1, xi, dth, th0, costf);

    end

    % Initialization
    idx_1 = 0;

    % In(De)creased values
    cm1 = cost_pe(r1, xi, dth, th0 - delta_th, costf);
    cp1 = cost_pe(r1, xi, dth, th0 + delta_th, costf);

    % Finite differences
    d1 = cp1 - cm1;
    d2 = cp1 + cm1 - 2 * cc;

    % Newton's optimization steps
    while (idx_1 < nm_nt && abs(d1) > ep_nt)

        % Current increase
        if d2 > 0
            dth0 = 0.5 * delta_th * (-d1) / d2;
        else
            dth0 = -sign(d1) * sqrt(delta_th);
        end

        % New value
        th0 = max([min([th0 + dth0, thi(1, 2)]), thi(1, 1)]);

        % Re-computes in(de)creased values
        cc = cost_pe(r1, xi, dth, th0, costf);
        cm1 = cost_pe(r1, xi, dth, th0 - delta_th, costf);
        cp1 = cost_pe(r1, xi, dth, th0 + delta_th, costf);

        % Re-computes finite differences
        d1 = cp1 - cm1;
        d2 = cp1 + cm1 - 2 * cc;

        % Increases steps counter
        idx_1 = idx_1 + 1;

        figure(2)
        hold on
            plot(th0, cc, 'gs')
        hold off

    end

    % Output values
    d = det_pe(xi, dth, th0);
    p = ppr_pe(r1, dth, th0, d);
    e = e_pe(xi, d);

end

% Determinant
function d = det_pe(xi, dth, th0)
    d = cos(dth + th0) - xi * cos(th0);
end

% Semi-latus rectum
function p = ppr_pe(r1, dth, th0, d)
    p = r1 * (cos(dth + th0) - cos(th0))./ d;
end

% Eccentricity
function e = e_pe(xi, d)
    e = (xi - 1)./ d;
end

% Cost function
function c = cost_pe(r1, xi, dth, th0, costf)
    d = det_pe(xi, dth, th0);
    c = costf(ppr_pe(r1, dth, th0, d), e_pe(xi, d), th0, dth);
end
