function [R] = rOrb(lan, incl, aop)

  R = rMat(lan, 3) * rMat(incl, 1) * rMat(aop, 3);

end
