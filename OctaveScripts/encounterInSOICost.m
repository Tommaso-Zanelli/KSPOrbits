function [d, t_cp, o1_pm, o2_pm, te_pm] = encounterInSOICost(mv, tm, ...
                            tf, te, o1, o2, oA, oT, mu_1, mu_2, rSOI)
%
%   function [d, t_cp, o1_pm, o2_pm, te_pm] = encounterInSOICost(mv, tm, ...
%                               tf, te, o1, o2, oA, oT, mu_1, mu_2, rSOI)
%
%   An attractor of standard gravitational parameter mu_2 is orbiting an
%   attractor of standard gravitational parameter mu_1 with orbit oA, its
%   sphere of influence has size rSOI. We'll denote them as attractor 2 and
%   attractor 1 respectively.
%   A target object is orbiting attractor 2 with orbit oT.
%   An object is orbiting attractor 1 with orbit o1, which will encounter
%   attractor 2's sphere of influence at time te, and proceed around it
%   with orbit o2.
%   A manoeuvre with delta_v described by mv is applied at time tm. If tm
%   is less than te, the following orbit around attractor 1, o1_pm, is
%   computed, its encounter with attractor 2 found (if no encounter
%   happens, d is set to infinite) and its orbit around attactor 2, o2_pm,
%   is computed.
%   If tm is larger than te, the subsequent orbit around attractor 2,
%   o2_pm, is computed (and o1_pm is set to null/invalid).
%   The closest passage along orbit 02_pm to the object of orbit oT is
%   found, the minimum distance is assigned to d and the time of passage to
%   t_cp.
%   A maximum time tf is also provided.
%

    % Builds the manoeuvre as a stucture
    m.t  = tm;
    m.v  = mv;

    % ---------------------------------------------------------------------
    % If the manoeuvre takes place before encounter with attractor 2:
    if tm < te

        % Applies the manoeuvre
        [o1_pm] = applyManoeuvre(o1, m, mu_1);

        % Finds the encounter
        [te_pm, ~, th_e1, th_eA] = entryInSphere([tm, tf], rSOI, ...
                                                        o1_pm, oA, mu_1);

        % Finds position and velocity of the object and of attractor 2 at
        % the encounter:
        [xx_e1, vv_e1, ~] = structArgToPosVelTime(th_e1, o1_pm, mu_1);
        [xx_eA, vv_eA, ~] = structArgToPosVelTime(th_eA, oA, mu_1);

        % Computes the post-manoeuvre orbit around attractor 2
        [~, o2_pm] = posVelToStruct((xx_e1 - xx_eA), (vv_e1 - vv_eA), ...
            mu_2, te_pm);


    % ---------------------------------------------------------------------
    % If the manoeuvre takes place after encounter with attractor 2:
    else

        % Encounter time is the one provided
        te_pm = te;

        % Orbit 1 is set as null/invalid
        [o1_pm] = rawOrbitToStruct([1,0,0,0,0],0);

        % Applies the manoeuvre
        [o2_pm] = applyManoeuvre(o2, m, mu_2);

    end
    % ---------------------------------------------------------------------

    % Computes the closest passage to the target object
    [t_cp, d, ~, ~] = closestPass([te_pm, tf], o2_pm, oT, mu_2);


end