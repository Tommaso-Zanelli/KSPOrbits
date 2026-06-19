% Additional arguments for "timeToArg" if unspecified
function [tta] = check_tta(tta, e)

    % Adds non-existing fields if not present
    if ~isfield(tta, 'e_bs')
        tta.e_bs = [];
    end
    if ~isfield(tta, 'e_nt')
        tta.e_nt = [];
    end
    if ~isfield(tta, 'n_bs')
        tta.n_bs = [];
    end
    if ~isfield(tta, 'n_nt')
        tta.n_nt = [];
    end
    if ~isfield(tta, 'ep')
        tta.ep = [];
    end
    if ~isfield(tta, 'd')
        tta.d = [];
    end

    % Fills missing arguments for "timeToArg"
    if numel(tta.e_bs) == 0
        tta.e_bs = min([1e-3, 0.1 / e]);            % The target residual for the bisection phase
    end
	if numel(tta.e_nt) == 0
        tta.e_nt = 1e-9;                            % The target residual for the Newton phase
    end
	if numel(tta.n_bs) == 0
        tta.n_bs = 250;                             % The maximum number of steps for the bisection phase
    end
	if numel(tta.n_nt) == 0
        tta.n_nt = 2500;                            % The maximum number of steps for the Newton phase
	end
	if numel(tta.ep) == 0
    	tta.ep = 1e-4;                              % The size of the interval arount pi which is to be linearized in argToTime
    end
	if numel(tta.d) == 0
        tta.d = 16;                                 % The negative esponent of the power of two representing the interval around e = 1 which is to be linearized  in argToTime
    end

end