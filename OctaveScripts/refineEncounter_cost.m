function [cc, o1, o2, t_enc] = refineEncounter_cost(a, ref, oA, entry, margCost, att)
%
%   function [cc, o1, o2, t_enc] = refineEncounter_cost(a, ref, oA, entry)
%
%   Data in ref:
%   [pe1        ap1     e1      lan1        incl1	aop1    R_1
%    pe1_s      pe1_y	pe1_d	pe1_hh      pe1_mm	pe1_ss  mu_1
%    now1_s     now1_y	now1_d	now1_hh     now1_mm	now1_ss R_2     ]
%   [enc_s      enc_y	enc_d	enc_hh      enc_mm	enc_ss  mu_2
%    nowE_s     nowE_y	nowE_d	nowE_hh     nowE_mm	nowE_ss rSOI_2  ]
%   [pe2        ap2     e2      lan2        incl2	aop2    0
%    pe2_s      pe2_y	pe2_d	pe2_hh      pe2_mm	pe2_ss  0
%    now2_s     now2_y	now2_d	now2_hh     now2_mm	now2_ss 0       ]
%
%   Data in a:
%   [d_pe1;     d_ap1;	d_lan1; d_incl1;	d_aop1;	d_peT1; ...
%	 d_pe2;     d_ap2;	d_lan2; d_incl2;	d_aop2;	d_peT2          ]
%
%

    % Additional arguments for "argToTime"
    if nargin < 6
    	att = [];
    end
    att = check_att(att);

    % By default margin costs are not added
    if nargin < 5
        margCost = false;
    end

    % If unspecified, entry is assumed
    if nargin < 4
        entry = true;
    else
        if numel(entry) < 1
            entry = true;
        end
    end

    % Builds the orbit around the attractor
    peap2_p     = 10.^ (floor(log10(abs(ref(6, 1:2)))) - 5);     % Periapsis and apoapsis precision
    ref(6, 1:2)	= ref(6, 1:2) + peap2_p.* (a(7:8, 1)');     % Corrects apoapsis and periapsis for specified precision
    ref(6, 4:6)	= ref(6, 4:6) + 0.1 * (a(9:11, 1)');        % Corrects the angles for the specified precision
    [~, pet2_p] = ydhmsTimeToSeconds(ref(7, 2:6));          % Time precision 1
    [~, now2_t] = ydhmsTimeToSeconds(ref(8, 2:6));          % Time precision 2
    ref(8, 6)   = ref(8, 6) + (pet2_p + now2_t) * a(12, 1);	% Corrects time for precision

    % Converts the orbit around the attractor to a structure
    [o2] = rawOrbitToStructOrbit(ref(6:8, 1:6), ref(3, 7));

    % Builds the orbit around the main body
    peap1_p     = 10.^ (floor(log10(abs(ref(1, 1:2)))) - 5);     % Periapsis and apoapsis precision
    ref(1, 1:2)	= ref(1, 1:2) + peap1_p.* (a(1:2, 1)');     % Corrects apoapsis and periapsis for specified precision
    ref(1, 4:6)	= ref(1, 4:6) + 0.1 * (a(3:5, 1)');         % Corrects the angles for the specified precision
    [~, pet1_p] = ydhmsTimeToSeconds(ref(2, 2:6));          % Time precision 1
    [~, now1_t] = ydhmsTimeToSeconds(ref(3, 2:6));          % Time precision 2
    ref(3, 6)   = ref(3, 6) + (pet1_p + now1_t) * a(6, 1);  % Corrects time for precision

    % Converts the orbit around the attractor to a structure
    [o1] = rawOrbitToStructOrbit(ref(1:3, 1:6), ref(1, 7));

    % Encounter argument
    th2_enc     = (1 - 2 * entry) *  acos(((o2.p / ref(5, 7)) - 1) / o2.e);

    % encounter time
    t_enc       = argStructToTime(th2_enc, o2, ref(4, 7), att);

    % Encounter time check
    [encT, encT_p] = ydhmsTimeToSeconds(ref(4, 2:6));	% Time precision 1
    [encN, encN_p] = ydhmsTimeToSeconds(ref(5, 2:6));	% Time precision 2
    b = [a; ((t_enc - (ref(4, 1) *  encT) - (ref(5, 1) * encN)) / (encT_p + encN_p))];	% Margins

    % Positions
    [xx1,   vv1,	~] = structTimeToPosVelArg(t_enc, o1, ref(2, 7));
    [xx2_1,	vv2_1,  ~] = structTimeToPosVelArg(t_enc, oA, ref(2, 7));
    [xx2_2,	vv2_2,	~] = structTimeToPosVelArg(t_enc, o2, ref(4, 7));
    xx2 = xx2_1 + xx2_2;
    vv2 = vv2_1 + vv2_2;

    % Base cost
    cc = (norm((xx1 - xx2), 2) ^ 2) + (norm((vv1 - vv2), 2) ^ 2);

    % Margin costs if requested
    if margCost
        cc = cc + sum(((2 * b - 0.5 * tanh(4 * b)).^ 4), 1);
    end

end