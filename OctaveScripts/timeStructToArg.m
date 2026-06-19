function [argument] = timeStructToArg(t, o, mu, tta)
%
% [argument] = timeStructToArg(t, o, mu, tta)
%
% Computes, iteratively, the orbital argument given the instant in time and
% the other orbital parameters. It applies a certain amout of bisections
% steps followed by Newton's method steps.
%

    % Additional arguments for "timeToArg" if unspecified
	if nargin < 4
        tta = [];
    end
    tta = check_tta(tta, o.e);

    % Computes the argument
    [argument] = timeToArg(t, o.p, o.e, mu, o.peT, ...
        tta.e_bs, tta.e_nt, tta.n_bs, tta.n_nt, tta.ep, tta.d);

end