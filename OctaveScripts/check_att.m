% Additional arguments for "argToTime"
function [att] = check_att(att)

    if ~isfield(att, 'ep')
        att.ep = [];
    end
    if ~isfield(att, 'd')
        att.d = [];
    end
    if numel(att.ep) == 0
    	att.ep = 1e-4;              % The size of the interval arount pi
                                    % which is to be linearized
    end
    if numel(att.d) == 0
        att.d = 16;                 % The negative esponent of the power of
                                    % two representing the interval around
                                    % e = 1 which is to be linearized
    end

end