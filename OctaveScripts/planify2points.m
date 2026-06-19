function [R, xr1, xr2] = planify2points(x1, x2, pd)
%
%   function [R, xr1, xr2] = planify2points(x1, x2, pd)
%
%
%   Produces the rotation matrix that brings points x1 and x2 on the x-y
%   plane, with x1 along the x axis.
%   The y direction is determined either by x2 or, if present, by the
%   optional preferential direction pd.
%   If requested, the rotated points are provided as output.
%
%

    % x direction
    i1 = x1 / norm(x1, 2);

    % z direction
    i3 = cross(x1, x2);
    if norm(i3, 2) == 0
        i3 = [0; 0; 1];
        i3 = i3 - i1 * dot(i3, i1);
    end
    while norm(i3, 2) == 0
        i3 = rand(3, 1);
        i3 = i3 - i1 * dot(i3, i1);
    end
    i3 = i3 / norm(i3, 2);

    % y direction
    i2 = cross(i3, i1);
    i2 = i2 / norm(i2, 2);

    % Corrects for the preferential direction if requested
    if nargin > 2
        if dot(i2, pd) < 0
            i2 = -i2;
            i3 = -i3;
        end
    end

    % Rotational matrix
    R = [i1, i2, i3]';

    % Additional output if requested
    if nargout > 1
        xr1 = R * x1;
    end
    if nargout > 2
        xr2 = R * x2;
    end

end