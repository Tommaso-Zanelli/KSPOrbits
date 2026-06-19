function [th_t1, th_t2, p_t, e_t, lan_t, incl_t, aop_t, er] = optiTransfer(x1, x2, costf, mu, delta_th, n_s, nm_nt, ep_ntS)
%
%   function [dv1, dv2, th_t1, th_t2, p_t, e_t, lan_t, incl_t, aop_t] = optiTransfer(x1, x2, v1, v2, costf, mu)
%
%   Given two points around a body, and a cost function, finds the optimal
%   transfer orbit. The cost function must be in the form costf(vt1, vt2)
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
        n_s = 36;
    end

    % If unspecified, sets the delta_th
    if nargin < 5
        delta_th = 4e-6;
    end

    % Moves on the plane of the two points selected
    i1 = x1 / norm(x1, 2);
    i3 = cross(x1, x2);
    if norm(i3, 2) == 0
        i3 = [0; 0; 1];
        i3 = i3 - i1 * dot(i3, i1);
    end
    while norm(i3, 2) == 0
        i3 = rand(3, 1);
        i3 = i3 - i1 * dot(i3, i1);
    end
    i3 = i3 / norm(i3, 2);
    i2 = cross(i3, i1);

    % Frame of reference
    T = [i1, i2, i3]';

    % The two points in the new frame of reference
    x1p = T * x1;
    x2p = T * x2;

    % The angle between the two points
    dth = atan2(norm(cross(x1p, x2p), 2), dot(x1p, x2p));

    %{
    % Coplanar velocities
    vx = @(mu, p, th)   -sqrt(mu / p) * sin(th);
    vy = @(mu, p, e, th) sqrt(mu / p) * (e + cos(th));

    %{
    vx2 = @(mu, p, e, th0, dth) vx(mu, p, (th0 + dth)) * cos(th0) + vy(mu, p, e, (th0 + dth)) * sin(th0);
    vy2 = @(mu, p, e, th0, dth) -vx(mu, p, (th0 + dth)) * sin(th0) + vy(mu, p, e, (th0 + dth)) * cos(th0);
    %}

    vx2 = @(mu, p, e, th) vx(mu, p, th);
    vy2 = @(mu, p, e, th) vy(mu, p, e, th);
    %}

    % Velocities, n-th attempt
    vr = @(mu, p, e, th) sqrt(mu / p) * e * sin(th);
    vth = @(mu, p, e, th) sqrt(mu / p) * (1 + e * cos(th));

    vx1 = @(mu, p, e, th) vr(mu, p, e, th);
    vy1 = @(mu, p, e, th) vth(mu, p, e, th);

    vx2 = @(mu, p, e, th) vr(mu, p, e, th) * cos(dth) - vth(mu, p, e, th) * sin(dth);
    vy2 = @(mu, p, e, th) vr(mu, p, e, th) * sin(dth) + vth(mu, p, e, th) * cos(dth);

    % Re-defines the velocities
    v1f = @(mu, p, e, th) ((T') * [vx1(mu, p, e, th); vy1(mu, p, e, th); 0]);
    v2f = @(mu, p, e, th) ((T') * [vx2(mu, p, e, th); vy2(mu, p, e, th); 0]);

    % Re-defines the cost function
    cost_n = @(p, e, th0, dth) costf(v1f(mu, p, e, th0), v2f(mu, p, e, (th0 + dth)));

    % Finds the optimal transfer orbit
    [p, e, th0] = optiTransfer_pe2(norm(x1p, 2), norm(x2p, 2), dth, cost_n, delta_th, n_s, nm_nt, ep_nt);

    % Finds the velocities at the two points
    v1t = v1f(mu, p, e, th0);
    v2t = v2f(mu, p, e, (th0 + dth));


    % If the error is requested
    if nargout > 8

	% Finds all the orbital parameters
    	[th_t1, p_t, e_t, lan_t, incl_t, aop_t, er1] = posVelToParams(x1, v1t, mu, 1);
    	[th_t2, p_tc, e_tc, lan_tc, incl_tc, aop_tc, er2] = posVelToParams(x2, v2t, mu, 1);

    	% Re-evaluates the error
	er = max([er1, er2, abs((p_tc - p_t) / p_t), ...
                            abs((p - p_t) / p_t), ...
                            abs((e_tc - e_t)), ...
                            abs((e - e_t)), ...
                            abs((lan_tc - lan_t) / (2 * pi)), ...
                            abs((incl_tc - incl_t) / (2 * pi)), ...
                            abs((aop_tc - aop_t) / (2 * pi))]);
    else

	% Finds all the orbital parameters

    	[th_t1, p_t, e_t, lan_t, incl_t, aop_t] = posVelToParams(x1, v1t, mu, 1);
    	[th_t2, p_tc, e_tc, lan_tc, incl_tc, aop_tc] = posVelToParams(x2, v2t, mu, 1);
        if th_t1 > th_t2
            th_t1 = th_t1 - 2 * pi;
        end
        fprintf('p   : %g, %g, %g\n', p, p_t, p_tc)
        fprintf('e   : %g, %g, %g\n', e, e_t, e_tc)
        fprintf('lan : %g, %g\n', lan_t, lan_tc)
        fprintf('incl: %g, %g\n', incl_t, incl_tc)
        fprintf('aop : %g, %g\n', aop_t, aop_tc)
        fprintf('cc  : %g\n', cost_n(p, e, th0, dth))
    %{
    th_t1 = 0;
    th_t2 = dth;
    p_t = p;
    e_t = e;
    lan_t  = 0;
    incl_t = -pi;
    aop_t  = th0;
    %}
    end

end
