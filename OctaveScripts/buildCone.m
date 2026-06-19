function [p1, p2, th1, ph1, th2, ph2] = buildCone(x, R, N1, N2, N3)
%
%   function [p1, p2, th1, ph1, th2, ph2] = buildCone(x, R, N1, N2, N3)
%
%   Builds a scatter of points on a spherical surface that are projected by
%   a "distant point" x
%

    if nargin < 3
        N1 = 4096;
    end

    % Distance of the point from the centre of the sphere
    d = norm(x, 2);

    % Unit vector from x to the centre of the sphere
    i1 = -x / d;

    % Builds a frame of reference at x
    i2 = zeros(size(i1));
    [~, i_mc] = min(abs(i1));
    i2(i_mc(1, 1), 1) = 1;
    i2 = i2 - i1 * dot(i1, i2);
    i2 = i2 / norm(i2, 2);
    i3 = cross(i1, i2);

    % Gets the max tangent distance of x from the sphere
    q = sqrt((d ^ 2) - (R ^ 2));

    % Gets the angle between x and the tangent to the sphere
    alpha = atan2(R, q);

    % Angles used to build the points
    theta1 = linspace(-pi, pi, N1);

    % Directions for the boundary
    ii1 = repmat(i1, 1, N1) * cos(alpha) + ...
         (repmat(i2, 1, N1).* sin(repmat(theta1, 3, 1)) + ...
          repmat(i3, 1, N1).* cos(repmat(theta1, 3, 1))) * sin(alpha);

    % Boundary points
    p1 = repmat(x, 1, N1) + q * ii1;

    % If the internal points are required
    if nargout > 1

        % Default numbers of points
        if nargin < 5
            N3 = 128;
        end
        if nargin < 4
            N2 = 128;
        end

        % Angles used to build the points
        theta2 = repmat(linspace(-pi, pi, N2), 1, N3);
        N23 = N2 * N3;
        al2 = reshape(repmat(sin(linspace(0, 0.5 * pi, N3)), N2, 1), 1, N23) * alpha;

        % Directions for the set of points
        ii2 = repmat(i1, 1, N23).* cos(repmat(al2, 3, 1)) + ...
             (repmat(i2, 1, N23).* sin(repmat(theta2, 3, 1)) + ...
              repmat(i3, 1, N23).* cos(repmat(theta2, 3, 1))).* ...
              sin(repmat(al2, 3, 1));

        % Dot product
        dprod = sum((ii2.* repmat(x, 1, N23)), 1);

        % Possible distances
        a = repmat(-dprod, 2, 1) + repmat([1; -1], 1, N23).* ...
            sqrt((repmat(dprod, 2, 1).^ 2) - ((d ^ 2) - (R ^ 2)));

        % Minimum values (first encounter with the sphere)
        [~, iam] = min(abs(a));

        % Distance vector
        qq = zeros(1, N23);
        for idx_1 = 1:N23
            qq(1, idx_1) = a(iam(1, idx_1), idx_1);
        end

        % Real values
        is_a_v   = (imag(qq) == 0);

        % Internal points
        p2 = repmat(x, 1, sum(is_a_v, 2)) + repmat(qq(1, is_a_v), 3, 1).* ii2(:, is_a_v);

        % If the polar coordinates are required for the boundary
        if nargout > 2

            % Computes the polar coordinates
            [th1, ph1, ~] = polarCoords(p1);

            % If the polar coordinates are required for the internal points
            if nargout > 4

                % Computes the polar coordinates
                [th2, ph2, ~] = polarCoords(p2);

            end

        end

    end


end