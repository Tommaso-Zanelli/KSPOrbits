function [X, Y, Z] = genSphere(R, N)

  [TH, PH] = meshgrid(linspace(0, 2 * pi, N), linspace(0, pi, N));
  X = R * cos(TH).* sin(PH);
  Y = R * sin(TH).* sin(PH);
  Z = R * cos(PH);

end
