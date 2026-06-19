function [om] = angularVelocity(th, p, e, mu)
%
% function [om] = angularVelocity(th, p, e, mu)
%
% As implied, computes the angular velocity of an object in orbit around a
% body of gravitational constant mu.
%

	om = sqrt(mu / (p ^ 3)) * ((1 + e * cos(th)).^ 2);

end
