function [cc] = Bashir_2_03_Cost02(a, eT, mu_K, mu_M, oOK, tC, oT, SOIr_M, oMK, R_M, oTM, t0)
%
%   [cc] = Bashir_2_03_Cost02(a, eT, mu_K, mu_M, oOK, tC, oT, SOIr_M, oMK, R_M, oTM, t0)
%
%   Whatever...
%

    c   = zeros(4, 1);
    tcp = zeros(4, 1);
    dv  = zeros(4, 1);
    for idx_1 = 1:4
        mm(idx_1).t  = t0 + (tC - t0) * a(idx_1, 1);
        idx_2 = 3 * idx_1 + 2;
        idx_3 = 3 * idx_1 + 4;
        mm(idx_1).v  = 1800 * a(idx_2:idx_3, 1) - 900;
        [c(idx_1, 1), tcp(idx_1, 1), dv(idx_1, 1)] = Bashir_2_03_Cost01(...
            [0; 0; 0], mm(idx_1), eT, mu_K, mu_M, oOK, tC, SOIr_M, oMK, ...
            R_M, oTM(idx_1));
        dv(idx_1, 1) = dv(idx_1, 1) + norm(mm(idx_1).v, 2);
        c(idx_1, 1)  = c(idx_1, 1) + costVel(dv(idx_1, 1));
    end

    cc = sum(c, 1) + timeCost([mm(1).t, mm(2).t, mm(3).t, mm(4).t, ...
                                tcp'], 180, oT);

end


% Costvel
function [cv] = costVel(dv)

    cv = (dv > 1800) * exp(dv - 1800);

end

% Timecost
function [ct] = timeCost(tt, dt, tm)

    ct = 0;

    for idx_1 = 1:size(tt, 2)
        for idx_2 = 1:size(tt, 2)
            if idx_1 ~= idx_2

                dtt = tt(1, idx_2) - tt(1, idx_1);
                if abs(dtt) < 1.01 * dt
                    ct = (1e9) * cc(ddt / dt);
                end

            end
        end
    end

    for idx_1 = 5:size(tt, 2)
        tl = tt(1, idx_1);
        ct = ct + (tl > tm) * ((tl - tm) ^ 4);
    end
end

function [c] = cc(x)

    if abs(x) < 1
        c = 1;
    else

    end

end