function [cf, th1t, th2t, otp, dv1, dv2, er] = transferCostStatic(prm, o1p, o2p, mu)
%
% function transferCostStatic(prm, o1p, o2p, mu_Ms)
%
% Computes the cost of a tranfer in terms of required deltaV

    [th1t, th2t, otp, dv1, dv2, er] = getTransferOrbitStatic(prm, o1p, o2p, mu);

    if otp(2, 1) < 0 || otp(1, 1) < 0
        cf = Inf;
        return;
    end

    cf = norm(dv1, 2) + norm(dv2, 2);

end

