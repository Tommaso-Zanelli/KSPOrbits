function [theta, phi, R] = polarCoords(x)
%
%   function [theta, phi] = polarCoords(x)
%
%   Computes the spherical coordinates of x.
%


    transpose = false;
    if size(x, 2) == 3
        transpose = true;
        x = x';
    end

    theta = atan2(x(2, :), x(1, :));
    if nargout > 1
    	phi = atan2(x(3, :), sqrt((x(1, :).^ 2) + (x(2, :).^ 2)));
        if nargout > 2
            R = sqrt((x(1, :).^ 2) + (x(2, :).^ 2) + (x(3, :).^ 2));
        end
	end

    if transpose
        theta = theta';
        if nargout > 1
            phi = phi';
            if nargout > 2
                R = R';
            end
        end
    end

end