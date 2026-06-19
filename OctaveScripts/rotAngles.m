function [lan, incl, aop_th] = rotAngles(x, v)
%
% function [lan, incl, aoppth] = rotAngles(x, v)
%
% Computes the rotational angles for the orbital plane, given a point  and
% the local velocity

    % Checks wether the input is oriented in the right way
    if size(x, 2) == 3 && size(x, 1) ~= 3
        x = x';
    end
    if size(v, 2) == 3 && size(v, 1) ~= 3
        v = v';
    end

    % Checks wether the input has the correct size, returns errors otherwise
    if size(x, 2) ~= size(v, 2) || size(x, 1) ~= 3
        fprintf('Runtime error: function rotAngles received input with inconsistent sizes.\n')
        lan    = NaN;
        incl   = NaN;
        aop_th = NaN;
        return;
    end

    if size(x, 2) > 1

        % Initializes the output containers
        lan    = zeros(1, size(x, 2));
        incl   = zeros(1, size(x, 2));
        aop_th = zeros(1, size(x, 2));

        % For each input element provided
        for idx_1 = 1:size(x, 2)

            % Recursive call to the function
            [lan(1, idx_1), incl(1, idx_1), aop_th(1, idx_1)] = rotAngles(x(:, idx_1), v(:, idx_1));

        end

    end

    ir     = x / norm(x, 2);                                                % radial unit vector
    iv     = v / norm(v, 2);                                                % velocity direction
    k      = cross(ir, iv);
    k      = k / norm(k, 2);                                                % unit vector normal to the orbital plane
    ith    = cross(k, ir);
    ith    = ith / norm(ith, 2);                                            % tangential unit vector

    if dot(iv, ith) < 0
        %fprintf('Interesting situation\n');
        k = -k;
        ith = -ith;
    end

    M      = [ir, ith, k];                                                  % reconstructed rotation matrix

    incl   = atan2(sqrt((M(3, 1) ^ 2) + (M(3, 2) ^ 2)), M(3, 3));           % orbital inclination
    %incl   = atan2(sqrt((M(1, 3) ^ 2) + (M(2, 3) ^ 2)), M(3, 3));

    % Whether inclination is null or not:
    if abs(incl) < eps(4)
        incl = 0;
        lan = 0;
        aop_th = atan2(M(2, 1), M(2, 2));
        %aop_th = atan2(-M(1, 2), M(1, 1));
    elseif abs(abs(incl) - pi) < eps(4)
        incl = pi;
        lan = 0;
        aop_th = atan2(M(2, 1), M(2, 2));
        %aop_th = atan2(-M(1, 2), M(1, 1));
    else
        lan    = atan2(M(1, 3), -M(2, 3));                                       % longitude of the ascending node
        aop_th = atan2(M(3, 1), M(3, 2));                                      % argument of the periapsis (plus argument)
    end

    % Sign check
    Mch = rOrb(lan, incl, aop_th);
    Ich = M' * Mch;

    if sign(Ich(1, 1)) == -1 && sign(Ich(2, 2)) == -1
        aop_th = aop_th + pi;
        aop_th = inMPiPiInt(aop_th);
    end


end