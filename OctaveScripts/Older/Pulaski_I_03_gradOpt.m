ar = a0;%[0.76910994764397909529, A * ones(1, (Nm1 + Nm2)), zeros(1, (Nm1 + Nm2))];
crs = cf(a0(1, 2:end), 382 * a0(1, 1));
dt = -0.125;
fprintf('\nSTARTING MINIMUM: %13.16f\n\n', crs)
prevred = false;

while abs(dt) >= 1e-3

    fprintf('New cycle...\n\n')
    tic
    a = ar;
    a(1, 1) = a(1, 1) + dt;
    cus = 0;
    cnt = true;
    while cnt

        cc = cf(a(1, 2:end), 382 * a(1, 1));
        ccr = cc;
        fprintf('Current cost is: %g\nComputing gradient...\n', cc)
        tic; g = Pulaski_I_03_grad2(cf, a, da, 1); toc
        dq = da;
        %q = -dq * g;

        %cn = cf(a(1, 2:end) + q(1, 2:end), 382 * (a(1, 1) + q(1, 1)));

        while dq >= da
            q = -dq * g;
            cn = cf(a(1, 2:end) + q(1, 2:end), 382 * (a(1, 1) + q(1, 1)));
            if cn < cc
               a = a + q;
               cc = cn;
               dq = 2 * dq;
            else
               dq = 0.5 * dq;
            end
        end

        %while cn < cc
        %dq = 2 * dq;
        %q = -dq * g;
        %cn = cf(a(1, 2:end) + q(1, 2:end), 382 * (a(1, 1) + q(1, 1)));
        %end
        %dq = 0.5 * dq;
        %q = -dq * g;
        %cn = cf(a(1, 2:end) + q(1, 2:end), 382 * (a(1, 1) + q(1, 1)));
        %if cn < cc
        %a = a + q;
        %else
        %cnt = false;
        %end

        if abs(cc - ccr) < (1e-6) * abs(ccr)
            cus = cus + 1;
            fprintf('New cycle in %g...\n', 10 - cus + 1)
        else
           cus = 0;
        end
        if cus > 10
            cnt = false;
        end
        if cc < crs
            fprintf('Alright!\n')
        end

    end

    if cc < crs
        ar = a;
        crs = cc;
        fprintf('\nNEW MINIMUM FOUND: %23.16f\n\n', crs);
        if ~prevred
            dt = 2 * dt;
        end
    else
        dt = 0.5 * dt;
        prevred = true;
    end
    toc
    fprintf('For this cycle.\n\n')

end

fprintf('a = [')
for idx_1 = 1:size(a, 2)
    if idx_1 < size(a, 2)
        sp = ', ';
    else
        sp = ']\n';
    end
    fprintf('%23.16e%s', ar(1, idx_1), sp)
end
