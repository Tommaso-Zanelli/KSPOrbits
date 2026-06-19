function [lan, incl, aop] = rotEquivalent(lan, incl, aop)

    i1_lan = (lan < 0);
    i2_lan = (lan >= 0);

    i1_aop = (aop < 0);
    i2_aop = (aop >= 0);

    lan(i1_lan) = lan(i1_lan) + pi;
    lan(i2_lan) = lan(i2_lan) - pi;

    incl = -incl;

    aop(i1_aop) = aop(i1_aop) + pi;
    aop(i2_aop) = aop(i2_aop) - pi;

end
