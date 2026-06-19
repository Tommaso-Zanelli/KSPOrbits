function [T] = argToTime(argument, p, e, mu, T_ap, ep, d)
%
% function [T] = argToTime(argument, p, e, mu, T_ap, ep, d)
%
% Computes the time given the argument and the orbital parameters
%

  if nargin < 5
    T_ap = 0;                       % the time at the apoapsis (offset)
  end
  if nargin < 6
      ep = 1e-4;                    % The size of the interval arount pi
                                    % which is to be linearized
  end
  if nargin < 7
    d = 16;                         % The negative esponent of the power of
                                    % two representing the interval around
                                    % e = 1 which is to be linearized
  end

  % Gets the sizes of the input array in order to reshape it
  nrows = size(argument, 1);
  ncolumns = size(argument, 2);
  n_el = nrows * ncolumns;

  argument = reshape(argument, n_el, 1);

  % Gets the equivalent circular angular velocity
  tc2p = sqrt((p ^ 3) / mu);


  if e < (1 - (2 ^ -d))
    % Elliptic orbit
    T = tc2p * argToTimeEllipse(argument, e, ep) + T_ap;
  elseif e >= (1 - (2 ^ -d)) && e <= 1
    % Linear combination between elliptic and parabolic orbit
    q1 = argToTimeEllipse(argument, (1 - (2 ^ (-d))), ep);
    q2 = argToTimeParabola(argument);
    q = q2 + (q1 - q2) * (1 - e) / (2 ^ -d);
    T = tc2p * q + T_ap;
  elseif e > 1 && e <= (1 + (2 ^ -d))
    % Linear combination between hyperbolic and parabolic orbit
    q1 = argToTimeHyperbola(argument, (1 + (2 ^ (-d))));
    q2 = argToTimeParabola(argument);
    q = q2 + (q1 - q2) * (e - 1) / (2 ^ -d);
    T = tc2p * q + T_ap;
  else
    % Hyperbolic orbit
    T = tc2p * argToTimeHyperbola(argument, e) + T_ap;
  end

  % Reshapes the output so that it has the same dimensions as the input
  T = reshape(T, nrows, ncolumns);

end
