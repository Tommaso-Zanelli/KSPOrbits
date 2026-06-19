function [p, e, th0] = optiTransfer_pe(r1, r2, dth, costf)
%
%   function [p, e] = optiTransfer_pe(r1, r2, dth, costf)
%
%   Given two co-planar points, their angular distance, and the velocities
%   required at them, computes the optimal transfer orbit between them.
%
%   costf must be in the form: c = costf(p, e, th0, dth)
%

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

    % Finds possible initial guesses
    [pp1, pn1, pp2, pn2] = posNegAngles(phi_d, phi_p);
    pp = 0.5 * (pp1 + pp2);
    pn = 0.5 * (pn1 + pn2);

    % Chooses the one for which "e" is positive
    if xi > 1
        th0 = pp;
    else
        th0 = pn;
    end

    % Initialization
    c1 = 1e308;
    c2 = cost_pe(r1, xi, dth, th0, costf);
    idx_1 = 1;
    delta_th = 1e-6;

    plot(linspace(-pi, pi, 501), cost_pe(r1, xi, dth, linspace(-pi, pi, 501), costf))
    hold on
        plot(th0, cost_pe(r1, xi, dth, th0, costf), 'ro')


    % Newton iterations (to be tweaked/made tweakable)
    while (((abs(c1 - c2) / abs(c1)) > 1e-3) && (idx_1 < 1000))

        % Surrounding values
        cm1 = cost_pe(r1, xi, dth, (th0 - delta_th), costf);
        cp1 = cost_pe(r1, xi, dth, (th0 + delta_th), costf);

        % Finite differences
        d1 = cp1 - cm1;
        d2 = cp1 + cm1 - 2 * c2;

        % If the second derivative is null, a small value is chosen
        if (d2 <= 0)
            d2 = max([delta_th ^ 2, abs(d2)]);
        end

        % Newton iteration
        th0 = th0 - 0.5 * delta_th * d1 / d2;
        plot(th0, cost_pe(r1, xi, dth, th0, costf), 'k*')
        pause

        % New cost
        c1 = c2;
        c2 = cost_pe(r1, xi, dth, th0, costf);

        % Increases counter
        idx_1 = idx_1 + 1;

    end
    hold off

    % Output values
    d = det_pe(xi, dth, th0);
    p = ppr_pe(r1, dth, th0, d);
    e = e_pe(xi, d);

end

function d = det_pe(xi, dth, th0)
    d = cos(dth + th0) - xi * cos(th0);
end

function p = ppr_pe(r1, dth, th0, d)
    p = r1 * (cos(dth + th0) - cos(th0))./ d;
end

function e = e_pe(xi, d)
    e = (xi - 1)./ d;
end

function c = cost_pe(r1, xi, dth, th0, costf)
    d = det_pe(xi, dth, th0);
    c = costf(ppr_pe(r1, dth, th0, d), e_pe(xi, d), th0, dth);
end
