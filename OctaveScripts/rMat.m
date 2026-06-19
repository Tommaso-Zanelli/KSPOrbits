function [R] = rMat(angle, axis)

  R = eye(3);
  v = [2, 3, 1, 2];
  q = v(1, axis:(axis + 1));
  R(q, q) = [cos(angle), -sin(angle); sin(angle), cos(angle)];

end
