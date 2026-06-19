function [cf, th1t, th2t, th0_3, otp, dv1, dv2, er] = transferCostStaticTimes(th1_1, th2_2, o1p, o2p, tAp1, tf, mu)
%
% [cf, th1t, th2t, otp, dv1, dv2, er] = transferCostStatic(prm, o1p, o2p, tAp1, tf, mu_Ms)
%
% Computes the cost of a tranfer in terms of required deltaV and penalizes
% times too close to a set of "forbidden" ones

    [cf, th1t, th2t, th0_3, otp, dv1, dv2, er] = transferCostStatOpt2p(th1_1, th2_2, o1p, o2p, mu);

    % Gets the time of the first manoeuvre
    [p1, e1, ~, ~, ~] = vectorToParams(o1p);
    t1 = argToTime(th1_1, p1, e1, mu, 0);

    % Gets the transfer time interval
    [p3, e3, ~, ~, ~] = vectorToParams(otp);
    th1tc = th1t;
    while th1tc > th2t
        th1tc = th1t - 2 * pi;
    end
    t1_3 = argToTime(th1tc, p3, e3, mu, 0);
    t2_3 = argToTime(th2t, p3, e3, mu, 0);
    Dt = t2_3 - t1_3;

    % Gets the time of the second manoeuvre
    t2 = t1 * Dt;

    % Scales the "forbidden times"
    tfs = tf - tAp1;

    % Adds the penalty for the forbidden times
    for idx_1 = 1:size(tfs, 1)
            cf = cf + 1000 * (exp(-(((t1 - tfs(idx_1, 1)) /  509.593080172811).^ 2)) + exp(-(((t2 - tfs(idx_1, 1)) /  509.593080172811).^ 2)));
    end

end

