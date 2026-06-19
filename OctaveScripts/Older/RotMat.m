fclose all;
close all
format long g
more off
clear
clc

symbAvail = true;

try
    syms lan incl aop
catch exception
    symbAvail = false;
end

if symbAvail

    Ml = [cos(lan),  -sin(lan),          0
          sin(lan),   cos(lan),          0
                 0,          0,          1];

    Mi = [       1,          0,          0
                 0,  cos(incl), -sin(incl)
                 0,  sin(incl),  cos(incl)];

    Ma = [cos(aop),  -sin(aop),          0
          sin(aop),   cos(aop),          0
                 0,          0,          1];

    Ml * Mi * Ma

else

    % Tested on 11/09/2019
    fprintf('[ cos(aop)*cos(lan) - cos(incl)*sin(aop)*sin(lan), - cos(lan)*sin(aop) - cos(aop)*cos(incl)*sin(lan),  sin(incl)*sin(lan)]\n');
    fprintf('[ cos(aop)*sin(lan) + cos(incl)*cos(lan)*sin(aop),   cos(aop)*cos(incl)*cos(lan) - sin(aop)*sin(lan), -cos(lan)*sin(incl)]\n');
    fprintf('[                              sin(aop)*sin(incl),                                cos(aop)*sin(incl),           cos(incl)]\n');

end