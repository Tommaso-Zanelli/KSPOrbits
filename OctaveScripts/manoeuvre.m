function [pf, ef, lanf, inclf, aopf, T_apf] = manoeuvre(tm, dv, p0, e0, lan0, incl0, aop0, mu, T_ap0, e_bs, e_nt, n_bs, n_nt, ep, d)
%
%   [pf, ef, lanf, inclf, aopf, T_apf] = manoeuvre(tm, dv, p0, e0, lan0, incl0, aop0, T_ap0, mu, e_bs, e_nt, n_bs, n_nt, ep, d)
%
%   Applies a manoeuvre
%

    if nargin < 9
        T_ap0 = 0;                              % the time at the apoapsis (offset)
    end
	if nargin < 10
        e_bs = min([1e-3, 0.1 / e0]);           % The target residual for the bisection phase
    end
	if nargin < 11
        e_nt = 1e-9;                            % The target residual for the Newton phase
    end
	if nargin < 12
        n_bs = 250;                             % The maximum number of steps for the bisection phase
    end
	if nargin < 13
        n_nt = 2500;                            % The maximum number of steps for the Newton phase
	end
	if nargin < 14
    	ep = 1e-4;                              % The size of the interval arount pi which is to be linearized in argToTime
    end
	if nargin < 15
        d = 16;                                 % The negative esponent of the power of two representing the interval around e = 1 which is to be linearized  in argToTime
    end

    % Gets the argument of the manoeuvre
    [a0m] = timeToArg(tm, p0, e0, mu, T_ap0, e_bs, e_nt, n_bs, n_nt, ep, d);

    % Gets the position of the manoeuvre
    [xm, v0m] = paramsToPosVel(a0m, p0, e0, lan0, incl0, aop0, mu);

    % Applies the velocity (this implies it is in the absolute frame of refernce)
    vfm = v0m + dv;

    % Gets the parameters of the new orbit
    [afm, pf, ef, lanf, inclf, aopf] = posVelToParams(xm, vfm, mu);

    % Gets the absolute time to/from periapsis in the new orbit
    [tf] = argToTime(afm, pf, ef, mu, 0, ep, d);

    % Gets the time to/from periapsis of the new orbit
    T_apf = tm - tf;

end

