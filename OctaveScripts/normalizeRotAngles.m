function [lan, incl, aop] = normalizeRotAngles(lan, incl, aop)

    % Size of the (smaller) input
    N = min([size(lan, 2), size(incl, 2), size(aop, 2)]);

    % No need to perform any action if only one input was given
    if N == 1
        return
    end

    % Reference values
    lan_ref  = lan(1, 1);
    incl_ref = incl(1, 1);
    aop_ref  = aop(1, 1);

    % Original values
    lan_1  = lan(1, 2:N);
    incl_1 = incl(1, 2:N);
    aop_1  = aop(1, 2:N);

    % Alternative equivalent values
    [lan_2, incl_2, aop_2] = rotEquivalent(lan_1, incl_1, aop_1);

    % Error values
    er_1 = abs(lan_1 - lan_ref) + abs(incl_1 - incl_ref) + abs(aop_1 - aop_ref);
    er_2 = abs(lan_2 - lan_ref) + abs(incl_2 - incl_ref) + abs(aop_2 - aop_ref);

    % Indexes
    id_1 = er_1 <= er_2;
    id_2 = er_2 <= er_1;

    % Re-initializes the output
    lan  = zeros(1, (N - 1));
    incl = zeros(1, (N - 1));
    aop  = zeros(1, (N - 1));

    % Assigns the values
    lan(id_1)  = lan_1(id_1);
    lan(id_2)  = lan_2(id_2);
    incl(id_1) = incl_1(id_1);
    incl(id_2) = incl_2(id_2);
    aop(id_1)  = aop_1(id_1);
    aop(id_2)  = aop_2(id_2);

    % Adds the first values
    lan  = [lan_ref, lan];
    incl = [incl_ref, incl];
    aop  = [aop_ref, aop];

end