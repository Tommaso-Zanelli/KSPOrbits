function [cf, th1t, th2t, th0_3, otp, dv1, dv2, er] = transferCostStatOpt2p(th1_1, th2_2, o1p, o2p, mu_Ms)
%
% function [cf, th1t, th2t, th0_3, otp, dv1, dv2, er] = transferCostStatOpt2p(th1_1, th2_2, o1p, o2p, tAp1, tf, mu_Ms)
%
% Uses the transferCostStatic function to compute, given two points, the
% optimal transfer orbit among all possible porbits connecting them

    % Initializes the "first guesses"
    th0_3s = -pi;
    th0_3f = th0_3s + (2 * pi) * 71 / 72;
    cst    = zeros(72, 1);

    for idx_1 = 1:10

        % (Re-)initializes the research interval
        th0_3v = linspace(th0_3s, th0_3f, 72)';

        % delta angle
        dth0_3 = (th0_3f - th0_3s) / 71;

        for idx_2 = 1:72

            % Computes the cost
            [cst(idx_2, 1)] = transferCostStatic([th1_1, th2_2, th0_3v(idx_2, 1)], o1p, o2p, mu_Ms);

        end

        % Finds the minimum cost
        ff = find(cst == min(cst));
        ff = ff(1, 1);

        % finds the extremes of the minimum interval
        th0_3s = th0_3v(ff, 1) - dth0_3;
        th0_3f = th0_3v(ff, 1) + dth0_3;

        % Gets the optimum
        th0_3 = th0_3v(ff, 1);

    end

    % Returns all parameters as requested
    [cf, th1t, th2t, otp, dv1, dv2, er] = transferCostStatic([th1_1, th2_2, th0_3], o1p, o2p, mu_Ms);

end