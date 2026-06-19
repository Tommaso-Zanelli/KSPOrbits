function [angle] = minusPiToZero(angle)

    angle(angle < 0) = angle(angle < 0) + 2 * pi;

end
