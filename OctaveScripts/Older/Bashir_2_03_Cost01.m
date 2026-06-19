function [cc, tcp, dv] = Bashir_2_03_Cost01(a, m0, eT, mu_K, mu_M, oOK, oT, SOIr_M, oMK, R_M, oTM)
%
%   [cc, tcp, dv] = Bashir_2_03_Cost01(a, m0, eT, mu_K, mu_M, oOK, oT, SOIr_M, oMK, R_M, oTM)
%
%
%

    % Epsilon value
    ep = 4e-6;

    % Initializes the cost
    cc = 0;

    % Checks whether the manoeuvre is already in the sphere of influence
    if m0.t < eT

        % Adds the current value to the initial guess
        m = m0;
        m.v = m.v + a;

        % Applies the manoeuvre
        opM = applyManoeuvre(oOK, m, mu_K);

        % Checks how close to Minmus it actually gets
        [~, d, ~, ~] = closestPass([m.t, oT], opM, oMK, mu_K);

        % If no Minmus encounter happens, cost is infinite
        if d > SOIr_M
            cc = Inf;
            return;
        end

        % Penalty if too close to Minmus
        if d > (1 - ep) * SOIr_M
            cc = cc + (((1.158e77) * ((d - (1 - ep) * SOIr_M)) ...
                / ((ep) * SOIr_M)) ^ 4);
        end

        % Finds Minmus encounter
        [teM, ~, ~, ~] = entryInSphere([m.t, oT], SOIr_M, opM, oMK, mu_K);

        % Then, the manoeuvre around Minmus is:
        [xOK1, vOK1, ~] = structTimeToPosVelArg(teM, opM, mu_K);
        [xMK1, vMK1, ~] = structTimeToPosVelArg(teM, oMK, mu_K);
        [~, oM] = posVelToStruct((xOK1 - xMK1), (vOK1 - vMK1), mu_M, teM, 1);

    else

        % Applies the manoeuvre
        oM = applyManoeuvre(oOM, m, mu_K);

        % Initial time is manoeuvre time
        teM = m.t;

    end

    % Checks whether the object will plummet into Minmus or not
    peM = oM.p / (1 + oM.e);
    apM = oM.p / (1 - oM.e);
    if peM < R_M

        % Finds Minmus encounter, this is the final time
        thMe = -acos((oM.p / R_M - 1) / oM.e);
        thMe = [thMe, -thMe];
        tMe   = argStructToTime(thMe, oM, mu_M);
        toM   = min(tMe(tMe > teM));

    else

        % If the orbit is closed and inside the SOI:
        if oM.e < 1 && apM < SOIr_M

            % The final time
            toM = 10 * argStructToTime(pi, oM, mu_M);

        else

            % The final time: exit from Minmus' SOI
            thMe = acos((oM.p / SOIr_M - 1) / oM.e);
            toM = 10 * argStructToTime(thMe, oM, mu_M);

        end

    end

    % Finds the closest it gets to the target object
    [tcp, d, ~, ~] = closestPass([teM, toM], oM, oTM, mu_M);

    % Adds the distance to the cost
    cc = cc + d;

    % Computes the velocity difference at closest pass
    [x1tmp, vo1, ~] = structTimeToPosVelArg(tcp, oM, mu_M);
    T = refFrameLocal(x1tmp, vo1);
    [~, vo2, ~] = structTimeToPosVelArg(tcp, oTM, mu_M);
    dv = (T') * (vo2 - vo1);

end

