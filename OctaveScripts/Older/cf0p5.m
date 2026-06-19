function c = cf0p5(x)

    c = (4.95e-4) * (exp(12.216 * ((0.5 + 0.5 * tanh(38.2 * (abs(x) - 0.5))).^ 2)) - 1);

end