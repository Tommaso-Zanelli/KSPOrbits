function angle = inMPiPiInt(angle)

    angle = angle - 2 * pi * (floor((angle + pi) / (2 * pi)));

end