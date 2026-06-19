function [argument] = timeToArg(t, p, e, mu, T_ap, e_bs, e_nt, n_bs, n_nt, ep, d)
%
% function [argument] = timeToArg(t, p, e, mu, T_ap, e_bs, e_nt, n_bs, n_nt, ep, d)
%
% Computes, iteratively, the orbital argument given the instant in time and
% the other orbital parameters. It applies a certain amout of bisections
% steps followed by Newton's method steps.

	if nargin < 5
        T_ap = 0;                               % the time at the apoapsis (offset)
    end
	if nargin < 6
        e_bs = min([1e-3, 0.1 / e]);            % The target residual for the bisection phase
    end
	if nargin < 7
        e_nt = 1e-9;                            % The target residual for the Newton phase
    end
	if nargin < 8
        n_bs = 250;                             % The maximum number of steps for the bisection phase
    end
	if nargin < 9
        n_nt = 2500;                            % The maximum number of steps for the Newton phase
	end
	if nargin < 10
    	ep = 1e-4;                              % The size of the interval arount pi which is to be linearized in argToTime
    end
	if nargin < 11
        d = 16;                                 % The negative esponent of the power of two representing the interval around e = 1 which is to be linearized  in argToTime
    end

    % Gets the sizes of the input array in order to reshape it
	nrows = size(t, 1);
	ncolumns = size(t, 2);
	n = nrows * ncolumns;

	t = reshape(t, n, 1);
	t = t - T_ap;

    % If the orbit is ellptic, corrects the time counting the number of
    % periods occurred
	nn = zeros(size(t));
    jumpSearch = false;
	if e < 1
        T = 2 * pi * sqrt(((p / (1 - (e ^ 2))) ^ 3) / mu);
        nn = floor((t + 0.5 * T) / T);
        t = t - T * nn;

        % If at the extremes, jumps the search to avoid NaN errors
        if abs((t / T) - 0.5) < eps(1)
            th = pi;
            jumpSearch = true;
        elseif abs((t / T) + 0.5) < eps(1)
            th = -pi;
            jumpSearch = true;
        end

    end


    if ~jumpSearch

        % Computes the limists of the search interval: the whole circumference
        % for an elliptic orbit, the valid arguments for an hyperbolic one
        if e <= 1
            b_0 = pi;
        else
            b_0 = acos((-1) / e);
        end

        % Sets up the search vectors
        a_0 =  -b_0;
        a = a_0 * ones(size(t));
        b = b_0 * ones(size(t));

        mv = 1e308;
        n_ib = 0;

        % Performs bisection
        while (mv > e_bs && n_ib < n_bs)
            c = [a, 0.5  * (a + b), b];
            er = argToTime(c, p, e, mu, 0, ep, d) - repmat(t, 1, 3);
            mv = max(abs(er(:, 2)));

            b(er(:, 1).* er(:, 2) > 0) = b(er(:, 1).* er(:, 2) > 0);
            a(er(:, 1).* er(:, 2) >= 0) = c(er(:, 1).* er(:, 2) >= 0, 2);

            b(er(:, 2).* er(:, 3) >= 0) = c(er(:, 2).* er(:, 3) >= 0, 2);
            a(er(:, 2).* er(:, 3) > 0) = a(er(:, 2).* er(:, 3) > 0);

            n_ib = n_ib + 1;
        end

        n_in = 0;
        th = c(:, 2);

        % Performs Newton
        while (mv > e_nt && n_in < n_nt)

            tt = t - argToTime(th, p, e, mu, 0, ep, d);
            mv = max(abs(tt));
            th = th + tt.* angularVelocity(th, p, e, mu);

            n_in = n_in + 1;

        end

    end

    argument = th;

    % Adds the rounds corresponding tho the periods already passed
	argument = argument + 2 * pi * nn;

    % Reshapes the output so that it has the same dimensions as the input
	argument = reshape(argument, nrows, ncolumns);

end