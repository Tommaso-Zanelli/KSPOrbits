function [pp, pn, pp2, pn2] = posNegAngles(th1, th2)
%
%   function [pp, pn] = posNegAngles(th1, th2)
%
%   Given two angles th1 and th2 in [-pi, pi] finds pp and pn such that:
%       -   cos(th1 + pp) > 0
%       -   cos(th2 + pp) > 0
%       -   cos(th1 + pn) < 0
%       -   cos(th2 + pn) < 0
%


    % Assigns the positive and negative intervals for (th + phi)
    c1 = (-th1 > -0.5 * pi) + (-th1 > 0.5 * pi);
    c2 = (-th2 > -0.5 * pi) + (-th2 > 0.5 * pi);
    intvt = zeros(3, 8);
    intvt([1 2], 1:4) = [pi * ((0:3) - 0.5 * (1 + 2 * c1)) - th1;
                         ([-1, 1, -1, 1]) * ((-1) ^ (c1 + 1))];
    intvt([1 3], 5:8) = [pi * ((0:3) - 0.5 * (1 + 2 * c2)) - th2 ;
                         ([-1, 1, -1, 1]) * ((-1) ^ (c2 + 1))];


    figure(3)
    hold on
        for idx_2 = 1:3
            idx_3 = idx_2 + 4;
            if intvt(2, idx_2) == 1
                c = 'bo-';
            else
                c = 'c';
            end
            plot(intvt(1, idx_2:(idx_2 + 1)), [0.5, 0.5], c)
            plot(intvt(1, idx_2:(idx_2 + 1)), [0.25, 0.25], 'go-')
            if intvt(3, idx_3) == 1
                c = 'ro-';
            else
                c = 'm';
            end
            plot(intvt(1, idx_3:(idx_3 + 1)), [0, 0], c)
            plot(intvt(1, idx_3:(idx_3 + 1)), [0.25, 0.25], 'go-')
            plot([intvt(1, idx_2 + 1), intvt(1, idx_2 + 1)], [-3, 3], 'k--')
            plot([intvt(1, idx_3 + 1), intvt(1, idx_3 + 1)], [-3, 3], 'k--')
            plot([intvt(1, idx_2), intvt(1, idx_2)], [-3, 3], 'k--')
            plot([intvt(1, idx_3), intvt(1, idx_3)], [-3, 3], 'k--')
        end
        plot([-pi, -pi], [-3, 3], 'y')
        plot([pi, pi], [-3, 3], 'y')
        ttot = linspace(min(intvt(1, :)), max(intvt(1, :)), 501);
        plot(ttot(cos(ttot + th1) > 0), cos(ttot(cos(ttot + th1) > 0) + th1)+0.5)
        plot(ttot(cos(ttot + th2) > 0), cos(ttot(cos(ttot + th2) > 0) + th2)+0.5)
        plot(ttot(cos(ttot + th1) < 0), cos(ttot(cos(ttot + th1) < 0) + th1))
        plot(ttot(cos(ttot + th2) < 0), cos(ttot(cos(ttot + th2) < 0) + th2))
    hold off

    % Reorders the boundaries of the intervals
    [~, ii] = sort(intvt(1, :), 2);
    intvt = intvt(:, ii);

    % Assigns signs for missing intervals
    for idx_1 = 2:8
        if intvt(2, idx_1) == 0
            intvt(2, idx_1) = intvt(2, (idx_1 - 1));
        end
        if intvt(3, idx_1) == 0
            intvt(3, idx_1) = intvt(3, (idx_1 - 1));
        end
    end

    % Finds positive and negative intervals
    ip = find(intvt(2, :) == 1 & intvt(3, :) == 1);
    in = find(intvt(2, :) == -1 & intvt(3, :) == -1);
    ip = ip(1, 1);
    in = in(1, 1);

    % Finds the midpoints of confirmed positive and negative intervals
    if nargout < 3
        pp = 0.5 * sum(intvt(1, ip:(ip + 1)), 2);
        pn = 0.5 * sum(intvt(1, in:(in + 1)), 2);
    else
        pp = intvt(1, ip);
        pn = intvt(1, in);
        pp2 = intvt(1, ip + 1);
        pn2 = intvt(1, in + 1);
    end

    figure(3)
    hold on
        plot([pp, pp2], [0.75, 0.75], 'ks-')
        plot([pn, pn2], [-0.25, -0.25], 's-', 'color', [0.5, 0.5, 0.5])
        posv_th = linspace(pp, pp2, 501);
        negv_th = linspace(pn, pn2, 501);

    hold off

end