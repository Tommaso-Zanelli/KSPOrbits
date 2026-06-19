function [angle] = zeroToMinusPi(angle)

    angle(angle > pi) = angle(angle > pi) - 2 * pi;

end
