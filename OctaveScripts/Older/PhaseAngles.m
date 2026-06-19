%% ------------------------------------------------------------------------

fclose all;
close all
format long g
more off
clc

%% ------------------------------------------------------------------------

th1s = -2:2;
th2s = th1s;

%th1s = -1
%th2s = 1

s = size(th1s, 2);

for idx_a = 1:s
    for idx_b = 1:s

        th1 = th1s(1, idx_a);

        th2 = th2s(1, idx_b);

        c = @(th) (-th > -0.5 * pi) + (-th > 0.5 * pi);
        intv = @(th) [-th + pi * ((0:3) - 0.5 * (1 + 2 * (c(th)))); ([-1 1 -1 1]) * ((-1) ^ (c(th) - 1))];

        c1 = c(th1);
        c2 = c(th2);

        intvt = zeros(3, 8);

        intvt([1 2], 1:4) = intv(th1);
        intvt([1 3], 5:8) = intv(th2);

        q = s * (idx_a - 1) + idx_b;
        figure(q)
        addToolbarExplorationButtons(gcf)
        hold on
        for idx_1 = 1:3
            x1 = intvt(1, (idx_1):(idx_1 + 1));
            y1 = [0, 0];
            if intvt(2, idx_1) == 1
                c = 'bo-';
            else
                c = 'c';
            end
            plot(-th1, 0, 'k*')
            plot(x1, y1, c)
            x2 = intvt(1, (idx_1 + 4):(idx_1 + 5));
            y2 = [-0.5, -0.5];
            if intvt(3, idx_1 + 4) == 1
                c = 'ro-';
            else
                c = 'm';
            end
            plot(x2, y2, c)
            plot(-th2, -0.5, 'ko')
            axis equal
            plot(intvt(1, (idx_1):(idx_1 + 1)), [-0.25, -0.25], 'go-')
            plot(intvt(1, (idx_1 + 4):(idx_1 + 5)), [-0.25, -0.25], 'go-')
            for idx_2 = [idx_1, (idx_1 + 1), (idx_1 + 4), (idx_1 + 5)]
                plot(repmat(intvt(1, idx_2), 1, 2), [-3, 3], 'k--')
            end
        end
        plot([-pi, -pi], [-3, 3], 'y')
        plot([pi, pi], [-3, 3], 'y')


        [~, ii] = sort(intvt(1, :), 2);

        intvt = intvt(:, ii);

        for idx_1 = 2:8
            if intvt(2, idx_1) == 0
                intvt(2, idx_1) = intvt(2, (idx_1 - 1));
            end
            if intvt(3, idx_1) == 0
                intvt(3, idx_1) = intvt(3, (idx_1 - 1));
            end
        end

        intvt(2, :) = (0.5 * sum(intvt(2:3, :), 1));
        %intvt(2, intvt(1, :) < - pi | intvt(1, :) > pi) = 0;

        ip = find(intvt(2, :) == 1);
        in = find(intvt(2, :) == -1);
        ip = ip(1, 1);
        in = in(1, 1);

        pp = 0.5 * sum(intvt(1, ip:(ip + 1)), 2);
        pn = 0.5 * sum(intvt(1, in:(in + 1)), 2);


        if pp > pi
            pp = pp - 2 * pi;
        end
        if pn > pi
            pn = pn - 2 * pi;
        end
        if pp < -pi
            pp = pp + 2 * pi;
        end
        if pn < -pi
            pn = pn + 2 * pi;
        end


        plot(pp, -0.25, 'ks')
        plot(pn, -0.25, 's', 'color', [0.5, 0.5, 0.5])
        hold off

        %{
        fprintf('%g:\n', s * (idx_a - 1) + idx_b)
        fprintf('cos(%g - pp) = %g\n', th1, cos(th1 - pp))
        fprintf('cos(%g - pp) = %g\n', th2, cos(th2 - pp))
        fprintf('cos(%g - pn) = %g\n', th1, cos(th1 - pn))
        fprintf('cos(%g - pn) = %g\n', th2, cos(th2 - pn))
        fprintf('\n')
        %}

        fprintf('%g:\n', s * (idx_a - 1) + idx_b)
        fprintf('cos(%g + pp) = %g\n', th1, cos(th1 + pp))
        fprintf('cos(%g + pp) = %g\n', th2, cos(th2 + pp))
        fprintf('cos(%g + pn) = %g\n', th1, cos(th1 + pn))
        fprintf('cos(%g + pn) = %g\n', th2, cos(th2 + pn))
        fprintf('\n')


    end
end