function C = MinmusCMap(N)

    load('MinmusCMap.mat');
    R = R(size(R, 1):(-1):1, :);
    G = G(size(G, 1):(-1):1, :);
    B = B(size(B, 1):(-1):1, :);

    C = zeros(N, N, 3);

    if N == size(R, 1)
        C(:, :, 1) = R;
        C(:, :, 2) = G;
        C(:, :, 3) = B;
    else
        [x1, y1] = meshgrid(linspace(0, 1, size(R, 1)), linspace(0, 1, size(R, 1)));
        [x2, y2] = meshgrid(linspace(0, 1, N), linspace(0, 1, N));
        C(:, :, 1) = interp2(x1, y1, R, x2, y2);
        C(:, :, 2) = interp2(x1, y1, G, x2, y2);
        C(:, :, 3) = interp2(x1, y1, B, x2, y2);
    end

end 