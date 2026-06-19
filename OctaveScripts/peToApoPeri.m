function [apo, peri] = peToApoPeri(p, e)

    apo  = p./ (1 + e);
    peri = p./ (1 - e);

end 